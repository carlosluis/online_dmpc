#ifndef BEZIER_H
#define BEZIER_H

#include <Eigen/Dense>
using namespace Eigen;

class BezierCurve {
public:
	BezierCurve(int deg, int l, int dim);
	~BezierCurve(){};
	// Public members
	VectorXd ctrl_pts;

private:
	Eigen::VectorXd get_sample(float t);
	int hola;
};

#endif
