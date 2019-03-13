#include <iostream>
#include <Eigen/Dense>
#include "bezier.h"

using namespace std;
using namespace Eigen;

int main() {
	cout << "Hello world!" << endl;
	BezierCurve test(5,3,3);
    test.set_ctrl_pts(VectorXd::Zero(54));
    cout << nchoosek(8,2) << endl;
    cout << test.bernstein_to_power(5) << endl;

	return 0;
}
