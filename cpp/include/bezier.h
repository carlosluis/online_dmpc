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

struct InequalityConstraint {
    Eigen::MatrixXd A_full;
    Eigen::VectorXd b_full;
    Eigen::MatrixXd A;
    Eigen::VectorXd lower_bound;
    Eigen::VectorXd upper_bound;
};

class BezierCurve {
public:
    struct Params{
        const int& deg, num_segments, dim, deg_poly;
        const float& t_segment;
    };

	BezierCurve(const BezierCurve::Params& p);
	~BezierCurve(){};

	// Public methods

    std::vector<Eigen::MatrixXd> getMatrixInputDerivativeSampling(const Eigen::VectorXd& t_samples);

    Eigen::MatrixXd getMatrixEnergyCost(const Eigen::VectorXd& weights);

    InequalityConstraint limitDerivative(int degree, const Eigen::VectorXd& t_samples,
                                         const Eigen::VectorXd& min, const Eigen::VectorXd& max);

    Eigen::MatrixXd getMatrixEqualityConstraint (int deg_poly);

private:
    // Methods
    Eigen::MatrixXd bernsteinToPowerBasis(int deg);
    Eigen::MatrixXd increaseMatrixDimension(const Eigen::MatrixXd& matrix, int dim);

    Eigen::MatrixXd getMatrixSumSqrdDerivatives(const Eigen::VectorXd& weights);

    Eigen::MatrixXd getMatrixPolynomeSampling(const Eigen::VectorXd& t_samples, int deg);

    Eigen::MatrixXd buildMatrixPolynomeSampling(const std::vector<Eigen::MatrixXd>& Tau,
                                                int num_samples, int deg);

    Eigen::MatrixXd augmentMatrixForm(const Eigen::MatrixXd& mat);

    Eigen::MatrixXd getMatrixInputSampling(const Eigen::VectorXd& t_samples, int r);

    std::vector<Eigen::MatrixXd> derivateControlPoints();

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
static Eigen::VectorXd shiftVector(Eigen::VectorXd vec, int n) {
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
