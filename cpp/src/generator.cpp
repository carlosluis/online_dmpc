//
// Created by carlos on 19/03/19.
//
#include "iostream"
#include "generator.h"

using namespace std;
using namespace Eigen;

Generator::Generator(Generator::Params* p) {
    float h = 0.2;
    float Ts = 0.01;

    VectorXd h_samples = VectorXd::LinSpaced(16, 0, 3);
    VectorXd ts_samples = VectorXd::LinSpaced(h / Ts, 0, h - Ts);

    _bezier = new BezierCurve(p->bezier_params);
    _model_pred = new DoubleIntegrator3D(h, p->model_params);
    _model_exec = new DoubleIntegrator3D(Ts, p->model_params);

    _H_snap = _bezier->get_mat_energy_cost(p->mpc_params->energy_weights);

    std::vector<MatrixXd> Phi_h = _bezier->get_vec_input_sampling(h_samples);
    std::vector<MatrixXd> Phi_ts = _bezier->get_vec_input_sampling(ts_samples);

    cout << ts_samples << endl;
    int hola2 = 1;

}
