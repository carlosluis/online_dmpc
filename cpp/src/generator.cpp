//
// Created by carlos on 19/03/19.
//
#include "iostream"
#include "generator.h"

using namespace std;
using namespace Eigen;

Generator::Generator(const Generator::Params& p) :
    _bezier(p.bezier_params),
    _model_pred(p.mpc_params.h, p.model_params),
    _model_exec(p.mpc_params.Ts, p.model_params){

    // Unpack params from struct into private members
    _h = p.mpc_params.h;
    _Ts = p.mpc_params.Ts;
    _k_hor = p.mpc_params.k_hor;
    _dim = p.bezier_params.dim;
    _l = p.bezier_params.num_segments;
    _po = p.po;
    _pf = p.pf;

    // Define two time bases:
    // 1) h_samples is the time base for the planning of trajectories
    // 2) ts_samples is a finer base for the execution in-between planning updates
    VectorXd h_samples = VectorXd::LinSpaced(_k_hor, 0, (_k_hor - 1) * _h);
    VectorXd ts_samples = VectorXd::LinSpaced(_h / _Ts, 0, _h - _Ts);

    // Energy cost function associated Hessian matrix
    _H_energy = _bezier.get_mat_energy_cost(p.mpc_params.tuning.energy_weights);

    _num_ctrl_pts = _H_energy.rows();

    // Get the matrices that sample the optimization result according to the 2 time bases
    std::vector<MatrixXd> Rho_h = _bezier.get_vec_input_sampling(h_samples);
    std::vector<MatrixXd> Rho_ts = _bezier.get_vec_input_sampling(ts_samples);

    // Get the matrices that propagate the input (i.e., optimization result) using the model
    StatePropagator Lambda_pred = _model_pred.get_lambda(_k_hor);
    StatePropagator A0_pred = _model_pred.get_A0(_k_hor);
    StatePropagator Lambda_exec = _model_exec.get_lambda(_h / _Ts);
    StatePropagator A0_exec = _model_exec.get_A0(_h / _Ts);

    // Create new matrices that map bezier curve control points to states of the robot
    _Phi_pred.pos = Lambda_pred.pos * Rho_h.at(0);
    _Phi_pred.vel = Lambda_pred.vel * Rho_h.at(0);
    _Phi_exec.pos = Lambda_exec.pos * Rho_ts.at(0);
    _Phi_exec.vel = Lambda_exec.vel * Rho_ts.at(0);

    // Build the inequality constraints regarding the physical limits of the robot
    _ineq = build_ineq_constr(p.mpc_params.limits);

    // Build Aeq matrix for the equality constraints
    _eq.A = _bezier.get_mat_eq_constr(p.bezier_params.deg_poly);

    // Set matrices to minimize goal error
    set_error_penalty_mats(p.mpc_params.tuning, _pf);

    cout << _H_f.block(0, 0, 15, 15) << endl << endl;
    cout << _H_o.block(0, 0, 15, 15) << endl << endl;
    cout << _H_r.block(0, 0, 15, 15) << endl;
}

void Generator::set_error_penalty_mats(const TuningParams& p, const MatrixXd& pf) {

    // Case 1: no collisions in the horizon - go fast
    MatrixXd S_free = MatrixXd::Zero(_dim * _k_hor, _dim * _k_hor);

    Ref<MatrixXd> block_f = S_free.block(_dim * (_k_hor - p.spd_f), _dim * (_k_hor - p.spd_f),
                                         _dim * p.spd_f, _dim * p.spd_f);

    block_f = p.s_free * MatrixXd::Identity(_dim * p.spd_f, _dim * p.spd_f);

    // Case 2: collisions in the horizon - go slower
    MatrixXd S_obs = MatrixXd::Zero(_dim * _k_hor, _dim * _k_hor);

    Ref<MatrixXd> block_o = S_obs.block(_dim * (_k_hor - p.spd_o), _dim * (_k_hor - p.spd_o),
                                         _dim * p.spd_o, _dim * p.spd_o);

    block_o = p.s_obs * MatrixXd::Identity(_dim * p.spd_o, _dim * p.spd_o);

    // Case 3: colliding scenario, change setpoints and move fast towards it

    MatrixXd S_repel = MatrixXd::Zero(_dim * _k_hor, _dim * _k_hor);

    Ref<MatrixXd> block_r = S_repel.block(0, 0, _dim * p.spd_r, _dim * p.spd_r);

    block_r = p.s_repel * MatrixXd::Identity(_dim * p.spd_r, _dim * p.spd_r);

    // Assemble matrices
    MatrixXd H_free = _Phi_pred.pos.transpose() * S_free * _Phi_pred.pos;
    MatrixXd H_obs = _Phi_pred.pos.transpose() * S_obs * _Phi_pred.pos;
    MatrixXd H_repel = _Phi_pred.pos.transpose() * S_repel * _Phi_pred.pos;

    // Get the final Hessians by multiplying by the energy consumption hessian
    _H_f = H_free + _H_energy;
    _H_o = H_obs + _H_energy;
    _H_r = H_repel + _H_energy;
}

Constraint Generator::build_ineq_constr(const PhysLimits& limits) {
    // Get position constraint
    VectorXd pos_samples = VectorXd::LinSpaced(_k_hor / 4, _h, _h + floor(_k_hor / 4) * _h * _l);
    Constraint lim_pos = _bezier.limit_derivative(0, pos_samples, limits.pmin, limits.pmax);

    // Get acceleration constraint
    VectorXd acc_samples = VectorXd::LinSpaced(_k_hor / 2, _h, (_k_hor - 1) * _h);
    Constraint lim_acc = _bezier.limit_derivative(2, acc_samples, limits.amin, limits.amax);

    // Assemble the constraints
    Constraint ineq;
    ineq.A = MatrixXd::Zero(lim_pos.A.rows() + lim_acc.A.rows(), lim_pos.A.cols());
    ineq.b = VectorXd::Zero(lim_pos.b.size() + lim_acc.b.size());

    ineq.A << lim_pos.A, lim_acc.A;
    ineq.b << lim_pos.b, lim_acc.b;

    return ineq;
}
