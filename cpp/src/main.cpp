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
    MpcParams mpc_params = {cr};

    // Testing the Generator class
    Generator::Params p = {&bezier_params, &model_params, &mpc_params};
    Generator gen(&p);

	return 0;
}
