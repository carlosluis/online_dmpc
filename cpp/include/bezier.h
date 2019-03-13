#ifndef BEZIER_H
#define BEZIER_H

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <unsupported/Eigen/KroneckerProduct>

class BezierCurve {
public:
	BezierCurve(int deg, int l, int dim);
	~BezierCurve(){};

	// Public methods
    void set_ctrl_pts(Eigen::VectorXd x);

    // Public variables
    Eigen::MatrixXd T_pts_poly;

    Eigen::MatrixXd bernstein_to_power(int deg);

private:
    // Methods

	Eigen::VectorXd get_sample(float t);
    Eigen::MatrixXd increase_matrix_dim(Eigen::MatrixXd matrix, int dim);

    // Variables
    int _deg;
    int _l;
    int _dim;
    int _num_ctrl_pts;
    Eigen::VectorXd _ctrl_pts;
};

// Helper function to shift an eigen vector n positions to the right, zero-padding
static Eigen::VectorXd shift_vector(Eigen::VectorXd vec, int n){
    int size = vec.size();
    Eigen::VectorXd shifted_vec = Eigen::VectorXd::Zero(size);
    shifted_vec.segment(n, size-n) = vec.segment(0, size-n);
    return shifted_vec;
};

// Helper function to calculate nchoosek
static int nchoosek(int n, int k)
{
    int res = 1;

    // Since C(n, k) = C(n, n-k)
    if ( k > n - k )
        k = n - k;

    // Calculate value of [n * (n-1) *---* (n-k+1)] / [k * (k-1) *----* 1]
    for (int i = 0; i < k; ++i)
    {
        res *= (n - i);
        res /= (i + 1);
    }

    return res;
};

#endif
