#include <iostream>
#include <Eigen/Dense>
#include "generator.h"

using namespace std;
using namespace Eigen;

int main() {
	cout << "Hello world!" << endl;

    // Bezier curve params
    int d = 5;
    int num_segments = 3;
    int deg_poly = 3;
    int dim = 3;
    float t_segment = 1.0;
    BezierCurve::Params bezier_params = {d, num_segments, dim, deg_poly, t_segment};

    // Model params
    float zeta_xy = 0.6502;
    float tau_xy = 0.3815;
    float zeta_z = 0.9103;
    float tau_z = 0.3;
    DoubleIntegrator3D::Params model_params = {zeta_xy, tau_xy, zeta_z, tau_z};

    // MPC params

    int T = 20; // sim duration
    int k_hor = 16; // horizon length
    float h = 0.2;
    float ts = 0.01;

    // penalty tuning
    TuningParams tune;
    tune.s_free = 100;
    tune.s_obs = 100;
    tune.s_repel = 1000;
    tune.spd_f = 3;
    tune.spd_o = 1;
    tune.spd_r = 10;
    tune.lin_coll = -pow(10, 5);
    tune.quad_coll = pow(10, 0);
    VectorXd cr = VectorXd::Zero(d + 1);
    cr(2) = .008;
    tune.energy_weights = cr;

    // Physical limits for inequality constraint
    PhysLimits limits;
    limits.pmin = (Eigen::Vector3d() << -1.5, -1.5, 0.2).finished();
    limits.pmax = (Eigen::Vector3d() << 1.5, 1.5, 2.2).finished();
    limits.amin = (Eigen::Vector3d() << -1.0, -1.0, -1.0).finished();
    limits.amax = (Eigen::Vector3d() << 1.0, 1.0, 1.0).finished();

    // Number of agents to control
    int N = 2;

    // Collision constraint for an agent
    EllipseParams ellipse_params;
    ellipse_params.order = 2;
    ellipse_params.rmin = 0.30;
    ellipse_params.c = (Eigen::Vector3d() << 1.0, 1.0, 2.0).finished();

    std::vector<EllipseParams> ellipse_vec(N, ellipse_params);

    MpcParams mpc_params = {h, ts, k_hor, tune, limits};

    // Generate a standard test for 1 vehicle moving 1.0 meters

    MatrixXd po = MatrixXd::Zero(3, N);
    Vector3d po1 = (Eigen::Vector3d() << -1.0, -1.0, 1.0).finished();
    Vector3d po2 = (Eigen::Vector3d() << -1.0, -0.8, 1.0).finished();
    po << po1, po2;

    MatrixXd pf = MatrixXd::Zero(3, N);
    Vector3d pf1 = (Eigen::Vector3d() << 1.0, 1.0, 1.0).finished();
    Vector3d pf2 = (Eigen::Vector3d() << -1.0, -1.0, 1.0).finished();
    pf << pf1, pf2;

    // Testing the Generator class
    Generator::Params p = {bezier_params, model_params, ellipse_vec, mpc_params, po, pf};
    Generator gen(p);

    State3D agent1 = {.pos = po1, .vel = 0.001*VectorXd::Ones(dim)};
    State3D agent2 = {.pos = po2, .vel = 0.001*VectorXd::Ones(dim)};
    std::vector<State3D> curr_states{agent1, agent2};
    gen.get_next_inputs(curr_states);

	return 0;
}

MatrixXd gen_rand_pts(const int &N,
                      const Vector3d &pmin,
                      const Vector3d &pmax,
                      const float &rmin) {
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