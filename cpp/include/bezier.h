#ifndef BEZIER_H
#define BEZIER_H

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <unsupported/Eigen/KroneckerProduct>

typedef Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,
                   Eigen::Dynamic, Eigen::RowMajor>> MapMatrixXd;

struct Constraint {
    Eigen::MatrixXd A;
    Eigen::VectorXd b;
};

class BezierCurve {
public:
    struct Params{
        int deg, num_segments, dim;
        float t_segment;
    };

	BezierCurve(BezierCurve::Params* p);
	~BezierCurve(){};

    // Public variables
    Eigen::MatrixXd T_pts_poly;

	// Public methods
    void set_ctrl_pts(Eigen::VectorXd x);

    std::vector<Eigen::MatrixXd> get_vec_input_sampling(Eigen::VectorXd t_samples);

    Eigen::MatrixXd get_mat_energy_cost(Eigen::VectorXd weights);

    std::vector<Eigen::MatrixXd> derivate_ctrl_pts();

    Constraint limit_derivative(int degree, Eigen::VectorXd t_samples,
                                Eigen::VectorXd min, Eigen::VectorXd max);

    Eigen::MatrixXd get_mat_eq_constr (int deg_poly);

private:
    // Methods
    Eigen::MatrixXd bernstein_to_power(int deg);
    Eigen::MatrixXd increase_matrix_dim(Eigen::MatrixXd matrix, int dim);

    Eigen::MatrixXd get_mat_sumsqrd_der(Eigen::VectorXd weights);

    Eigen::MatrixXd get_mat_poly_sampling(Eigen::VectorXd t_samples, int deg);

    Eigen::MatrixXd build_mat_poly_sampling(std::vector<Eigen::MatrixXd> Tau,
                                          int num_samples, int deg);

    Eigen::MatrixXd augmented_form(Eigen::MatrixXd mat);

    Eigen::MatrixXd get_mat_input_sampling(Eigen::VectorXd t_samples, int r);



    // Variables
    int _deg;
    int _num_segments;
    double _t_segment;
    int _dim;
    int _num_ctrl_pts;
    Eigen::VectorXd _ctrl_pts;
    Eigen::MatrixXd _Gamma;
    Eigen::MatrixXd _Beta;
    Eigen::MatrixXd _GammaBeta;
    std::vector<Eigen::MatrixXd> _T_ctrl_pts;
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
