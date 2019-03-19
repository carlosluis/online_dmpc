#include <iostream>
#include <Eigen/Dense>
#include "generator.h"

using namespace std;
using namespace Eigen;

int main() {
	cout << "Hello world!" << endl;

    BezierCurve::Params bezier_params = {5,3,3,1.0};
    DoubleIntegrator3D::Params model_params = {0.6502, 0.3815, 0.9103, 0.3};

    VectorXd t_samples = VectorXd::LinSpaced(16, 0, 3);

    //MPC tuning parameters
    VectorXd cr = VectorXd::Zero(5+1);
    cr(2) = .008;

    PhysLimits limits;
    limits.pmin = (Eigen::Vector3d() << -2.5, -2.5, 0.2).finished();
    limits.pmax = (Eigen::Vector3d() << 2.5, 2.5, 2.2).finished();
    limits.amin = (Eigen::Vector3d() << -2.0, -2.0, -2.0).finished();
    limits.amax = (Eigen::Vector3d() << 2.0, 2.0, 2.0).finished();


    MpcParams mpc_params = {0.2, 0.01, 16, &limits, cr};

    // Testing the Generator class
    Generator::Params p = {&bezier_params, &model_params, &mpc_params};
    Generator gen(&p);

	return 0;
}
