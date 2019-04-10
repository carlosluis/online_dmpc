#include <iostream>
#include <Eigen/Dense>
#include "simulator.h"

using namespace std;
using namespace Eigen;
using namespace std::chrono;

Simulator::Params getSimulationParams();

int main() {
	cout << "Hello world!" << endl;

    Simulator::Params sim_params = getSimulationParams();
    Simulator sim(sim_params);

    int T = 20; // simulation duration
//    sim.run(T);

//    State3D agent1 = {.pos = po1, .vel = 0.001*VectorXd::Ones(dim)};
//    State3D agent2 = {.pos = po2, .vel = 0.001*VectorXd::Ones(dim)};
//    State3D agent3 = {.pos = po3, .vel = 0.001*VectorXd::Ones(dim)};
//    State3D agent4 = {.pos = po4, .vel = 0.001*VectorXd::Ones(dim)};
//    std::vector<State3D> curr_states{agent1, agent2, agent3, agent4};
//
//    high_resolution_clock::time_point t1 = high_resolution_clock::now();
//
//    gen.getNextInputs(curr_states);
//
//    high_resolution_clock::time_point t2 = high_resolution_clock::now();
//    auto duration = duration_cast<microseconds>( t2 - t1 ).count();
//    cout << "Next inputs computed in = "
//         << duration/1000.0 << "ms" << endl << endl;

	return 0;
}

Simulator::Params getSimulationParams() {
    // In the future this function should be replaced by an appropriate read of a config file

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
    int N = 4;

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
    Vector3d po2 = (Eigen::Vector3d() << 1.0, 1.0, 1.0).finished();
    Vector3d po3 = (Eigen::Vector3d() << -1.0, 1.0, 1.0).finished();
    Vector3d po4 = (Eigen::Vector3d() << 1.0, -1.0, 1.0).finished();
    po << po1, po2, po3, po4;

    MatrixXd pf = MatrixXd::Zero(3, N);
    Vector3d pf1 = (Eigen::Vector3d() << 1.0, 1.0, 1.0).finished();
    Vector3d pf2 = (Eigen::Vector3d() << -1.0, -1.0, 1.0).finished();
    Vector3d pf3 = (Eigen::Vector3d() << 1.0, -1.0, 1.0).finished();
    Vector3d pf4 = (Eigen::Vector3d() << -1.0, 1.0, 1.0).finished();
    pf << pf1, pf2, pf3, pf4;

    // Testing the Generator class
    Generator::Params p = {bezier_params, model_params, ellipse_vec, mpc_params, po, pf};

//    Generator gen(p);

    // Noise levels from VICON system
    float std_position = 0.00228682;
    float std_vel = 0.0109302;

    Simulator::Params psim = {p, std_position, std_vel};

    return psim;
}

