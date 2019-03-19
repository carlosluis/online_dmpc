//
// Created by carlos on 19/03/19.
//
#include "iostream"
#include "generator.h"

using namespace std;
using namespace Eigen;

Generator::Generator(Generator::Params* p) {
    _h = p->mpc_params->h;
    _Ts = p->mpc_params->Ts;
    _k_hor = p->mpc_params->k_hor;
    _l = p->bezier_params->num_segments;

    VectorXd h_samples = VectorXd::LinSpaced(_k_hor, 0, (_k_hor - 1) * _h);
    VectorXd ts_samples = VectorXd::LinSpaced(_h / _Ts, 0, _h - _Ts);

    _bezier = new BezierCurve(p->bezier_params);
    _model_pred = new DoubleIntegrator3D(_h, p->model_params);
    _model_exec = new DoubleIntegrator3D(_Ts, p->model_params);

    _H_snap = _bezier->get_mat_energy_cost(p->mpc_params->energy_weights);

    std::vector<MatrixXd> Phi_h = _bezier->get_vec_input_sampling(h_samples);
    std::vector<MatrixXd> Phi_ts = _bezier->get_vec_input_sampling(ts_samples);

    Constraint ineq = build_ineq_constr(p->mpc_params->limits);
}

Constraint Generator::build_ineq_constr(PhysLimits *limits) {
    // Get position constraint
    VectorXd pos_samples = VectorXd::LinSpaced(_k_hor / 4, _h, _h + floor(_k_hor / 4) * _h * _l);
    Constraint lim_pos = _bezier->limit_derivative(0, pos_samples, limits->pmin, limits->pmax);

    // Get position constraint
    VectorXd acc_samples = VectorXd::LinSpaced(_k_hor / 2, _h, (_k_hor - 1) * _h);
    Constraint lim_acc = _bezier->limit_derivative(2, acc_samples, limits->amin, limits->amax);

    // Assemble the constraints
    Constraint ineq;
    ineq.A = MatrixXd::Zero(lim_pos.A.rows() + lim_acc.A.rows(), lim_pos.A.cols());
    ineq.b = VectorXd::Zero(lim_pos.b.size() + lim_acc.b.size());

    ineq.A << lim_pos.A, lim_acc.A;
    ineq.b << lim_pos.b, lim_acc.b;

    return ineq;
}
