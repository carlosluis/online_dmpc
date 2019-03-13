#include <iostream>
#include <Eigen/Dense>
#include "bezier.h"

int main() {
	std:: cout << "Hello world! \n";
	Eigen::Matrix<double, 3, 3> A;
	A.fill(10);
	std::cout << A << std::endl;
	BezierCurve test(0,0,0);
	return 0;
}
