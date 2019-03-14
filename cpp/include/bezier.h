#ifndef BEZIER_H
#define BEZIER_H

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <unsupported/Eigen/KroneckerProduct>

typedef Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,
                   Eigen::Dynamic, Eigen::RowMajor>> MapMatrixXd;

class BezierCurve {
public:
	BezierCurve(int deg, int l, int dim);
	~BezierCurve(){};

	// Public methods
    void set_ctrl_pts(Eigen::VectorXd x);

    // Public variables
    Eigen::MatrixXd T_pts_poly;

    Eigen::MatrixXd bernstein_to_power(int deg);
    void set_input_samples(float segment_duration, Eigen::VectorXd t_samples);
    Eigen::MatrixXd get_mat_sample_poly(float segment_duration, Eigen::VectorXd t_samples,
                                        int deg);

    Eigen::VectorXd get_sample(float t);

private:
    // Methods

    Eigen::MatrixXd increase_matrix_dim(Eigen::MatrixXd matrix, int dim);

    Eigen::MatrixXd build_sample_mat_poly(std::vector<Eigen::MatrixXd> Tau,
                                          int num_samples, int deg);


    // Variables
    int _deg;
    int _num_segments;
    int _dim;
    int _num_ctrl_pts;
    Eigen::VectorXd _ctrl_pts;
    Eigen::MatrixXd _Gamma;
    Eigen::MatrixXd _Beta;
    Eigen::MatrixXd _GammaBeta;
};

// Helper function to shift an Eigen vector n positions to the right, zero-padding
static Eigen::VectorXd shift_vector(Eigen::VectorXd vec, int n) {
    int size = vec.size();
    Eigen::VectorXd shifted_vec = Eigen::VectorXd::Zero(size);
    shifted_vec.segment(n, size-n) = vec.segment(0, size-n);
    return shifted_vec;
};

// Helper function to calculate nchoosek
static int nchoosek(int n, int k) {
    int res = 1;

    // Since C(n, k) = C(n, n-k)
    if ( k > n - k )
        k = n - k;

    // Calculate value of [n * (n-1) *---* (n-k+1)] / [k * (k-1) *----* 1]
    for (int i = 0; i < k; ++i) {
        res *= (n - i);
        res /= (i + 1);
    }

    return res;
};

#endif
