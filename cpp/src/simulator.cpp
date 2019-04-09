//
// Created by carlos on 09/04/19.
//
#include "iostream"
#include "simulator.h"

using namespace Eigen;
using namespace std;

Simulator::Simulator(const Simulator::Params& p) :
    _generator(p.generator_params),
    _pos_std(p.position_noise_std),
    _vel_std(p.velocity_noise_std),
    _sim_model(p.generator_params.mpc_params.Ts, p.generator_params.model_params){

}

MatrixXd Simulator::generateRandomPoints(int N, const Vector3d &pmin,
                                         const Vector3d &pmax, float rmin) {
    MatrixXd pts = MatrixXd::Zero(3, N);
    Vector3d candidate = MatrixXd::Zero(3, 1);
    VectorXd dist;
    bool pass = false;

    // Generate first point
    pts.col(0) = pmin.array()
                 + (pmax - pmin).array() *
                   ((MatrixXd::Random(3, 1).array() + 1) / 2);

    for (int n = 1; n < N; ++n) {
        while (!pass) {
            // Candidate picked randomly within workspace boundaries
            candidate = pmin.array()
                        + (pmax - pmin).array() *
                          ((MatrixXd::Random(3, 1).array() + 1) / 2);

            // Calculate distance to every previous pts calculated
            dist = ((((pts.leftCols(n)).colwise()
                      -
                      candidate).array().square()).colwise().sum()).array().sqrt();

            // If the candidate is sufficiently separated from previous pts,
            // then we add it to the Matrix of valid pts
            for (int k = 0; k < n; ++k) {
                pass = dist[k] > rmin;
                if (!pass)
                    break;
            }
            if (pass)
                pts.col(n) = candidate.array();
        }
        pass = false;
    }
    return pts;
}