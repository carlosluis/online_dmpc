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
    VectorXd cr = VectorXd::Zero(5+1);
    cr(2) = .008;

    VectorXd t_samples = VectorXd::LinSpaced(8, 0.2, 3);

    cout << t_samples << endl << endl;

	return 0;
}
