#include <iostream>
#include <Eigen/Dense>
#include "bezier.h"

using namespace std;
using namespace Eigen;

int main() {
	cout << "Hello world!" << endl;
	BezierCurve test(5,3,3);
    test.set_ctrl_pts(VectorXd::Zero(54));
//    cout << test.bernstein_to_power(5) << endl;

    VectorXd t_samples = VectorXd::LinSpaced(8, 0.2, 3);

    cout << t_samples << endl << endl;

    MatrixXd Gamma = test.get_mat_sample_poly(1.0, t_samples, 5-2);

	return 0;
}
