//
// Created by carlos on 19/03/19.
//

#include <iostream>
#include "model.h"

using namespace Eigen;
using namespace std;

StatePropagator DoubleIntegrator::get_lambda(int K) {
    StatePropagator Lambda;
    Lambda.pos = MatrixXd::Zero(_dim * K, _dim * K);
    Lambda.vel = MatrixXd::Zero(_dim * K, _dim * K);
    MatrixXd prev_row = MatrixXd::Zero(2 * _dim, _dim * K);
    MatrixXd new_row = MatrixXd::Zero(2 * _dim, _dim * K);
    MatrixXd add_b = MatrixXd::Zero(2 * _dim, _dim * K);

    int idx = 0;

    for (int k = 0; k < K; ++k) {
        add_b << MatrixXd::Zero(_model_B.rows(), _model_B.cols() * (k)), _model_B,
                 MatrixXd::Zero(_model_B.rows(), _model_B.cols() * (K - k - 1));
        new_row = _model_A * prev_row + add_b;
        Lambda.pos.middleRows(idx, 3) = new_row.middleRows(0, 3);
        Lambda.vel.middleRows(idx, 3) = new_row.middleRows(3, 3);
        prev_row = new_row;
        idx += 3;
    }

    return Lambda;
}

StatePropagator DoubleIntegrator::get_A0(int K) {
    StatePropagator A0;
    A0.pos = MatrixXd::Zero(_dim * K, 2 * _dim);
    A0.vel = MatrixXd::Zero(_dim * K, 2 * _dim);
    MatrixXd new_row = MatrixXd::Zero(2 * _dim, 2 * _dim);
    MatrixXd prev_row = MatrixXd::Identity(2 * _dim, 2 * _dim);

    int idx = 0;

    for (int  k = 0; k < K; ++k) {
        new_row = _model_A * prev_row;
        A0.pos.middleRows(idx, 3) = new_row.middleRows(0, 3);
        A0.vel.middleRows(idx, 3) = new_row.middleRows(3, 3);
        prev_row = new_row;
        idx += 3;
    }

    return A0;
}


DoubleIntegrator3D::DoubleIntegrator3D(const float& ts, const DoubleIntegrator3D::Params& p) {

    _dim = 3;
    float omega_xy = 1 / p.tau_xy;
    float omega_z = 1 / p.tau_z;

    // Create the 3D double integrator model
    _model_A = MatrixXd::Zero(6, 6);
    _model_B = MatrixXd::Zero(6, 3);
    _model_A << 1, 0, 0, ts, 0, 0,
                0, 1, 0, 0, ts, 0,
                0, 0, 1, 0, 0, ts,
                -ts * pow(omega_xy, 2), 0, 0, 1 - (2 * omega_xy * ts * p.zeta_xy), 0, 0,
                0, -ts * pow(omega_xy, 2), 0, 0, 1 - (2 * omega_xy * ts * p.zeta_xy), 0,
                0, 0, -ts * pow(omega_z, 2), 0, 0, 1 - (2 * omega_z * ts * p.zeta_z);

    _model_B << MatrixXd::Zero(3, 3),
                ts * pow(omega_xy, 2), 0, 0,
                0, ts * pow(omega_xy, 2), 0,
                0, 0, ts * pow(omega_z, 2);
}