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
    _l = p.bezier_params.num_segments;

    // Define two time bases:
    // 1) h_samples is the time base for the planning of trajectories
    // 2) ts_samples is a finer base for the execution in-between planning updates
    VectorXd h_samples = VectorXd::LinSpaced(_k_hor, 0, (_k_hor - 1) * _h);
    VectorXd ts_samples = VectorXd::LinSpaced(_h / _Ts, 0, _h - _Ts);

    // Energy cost function associated Hessian matrix
    _H_energy = _bezier.get_mat_energy_cost(p.mpc_params.energy_weights);

    // Get the matrices that sample the optimization result according to the 2 time bases
    std::vector<MatrixXd> Rho_h = _bezier.get_vec_input_sampling(h_samples);
    std::vector<MatrixXd> Rho_ts = _bezier.get_vec_input_sampling(ts_samples);

    // Get the matrices that propagate the input (i.e., optimization result) using the model
    StatePropagator Lambda_pred = _model_pred.get_lambda(_k_hor);
    StatePropagator A0_pred = _model_pred.get_A0(_k_hor);
    StatePropagator Lambda_exec = _model_exec.get_lambda(_h / _Ts);
    StatePropagator A0_exec = _model_exec.get_A0(_h / _Ts);

    // Create new matrices that map bezier curve control points to states of the robot
//    _Phi_pred.pos = Lambda_pred.pos

    // Build the inequality constraints regarding the physical limits of the robot
    build_ineq_constr(p.mpc_params.limits, &_ineq);

    _eq.A = _bezier.get_mat_eq_constr(p.bezier_params.deg_poly);
}

void Generator::build_ineq_constr(const PhysLimits& limits, Constraint* ineq) {
    // Get position constraint
    VectorXd pos_samples = VectorXd::LinSpaced(_k_hor / 4, _h, _h + floor(_k_hor / 4) * _h * _l);
    Constraint lim_pos = _bezier.limit_derivative(0, pos_samples, limits.pmin, limits.pmax);

    // Get acceleration constraint
    VectorXd acc_samples = VectorXd::LinSpaced(_k_hor / 2, _h, (_k_hor - 1) * _h);
    Constraint lim_acc = _bezier.limit_derivative(2, acc_samples, limits.amin, limits.amax);

    // Assemble the constraints
    ineq->A = MatrixXd::Zero(lim_pos.A.rows() + lim_acc.A.rows(), lim_pos.A.cols());
    ineq->b = VectorXd::Zero(lim_pos.b.size() + lim_acc.b.size());

    ineq->A << lim_pos.A, lim_acc.A;
    ineq->b << lim_pos.b, lim_acc.b;
}
