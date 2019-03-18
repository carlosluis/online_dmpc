#include <iostream>
#include "bezier.h"

using namespace std;
using namespace Eigen;

BezierCurve::BezierCurve(int deg, int num_segments, float t_segment, int dim){
	cout << "Creating a Bezier curve with " << num_segments << " segments of degree "
         << deg << " in " << dim << " dimensions" << endl;

    // Assign to private variables
    _deg = deg;
    _num_segments = num_segments;
    _t_segment = t_segment;
    _dim = dim;

    cout << "Initializing the vector of control points as zero" << endl;
    _num_ctrl_pts = dim*num_segments*(deg+1);
    _ctrl_pts = VectorXd::Zero(_num_ctrl_pts);

    // Build the matrices defined by the bezier curve parameters
    _Beta = bernstein_to_power(_deg);

    // Obtain matrices that compute the n-th derivative of position ctrl points
    _T_ctrl_pts = derivate_ctrl_pts();
}

MatrixXd BezierCurve::bernstein_to_power(int deg) {
    MatrixXd delta = MatrixXd::Zero(deg + 1, deg + 1);
    for (int k = 0; k < deg + 1; ++k){
        for (int i = k; i <= deg; ++i){
            delta(i, k) = pow(-1, i - k) * nchoosek(deg, i) * nchoosek(i, k);
        }
    }

    return augmented_form(delta);
}

MatrixXd BezierCurve::increase_matrix_dim(Eigen::MatrixXd matrix, int dim) {
    int nrows = (int) matrix.rows();
    int ncols = (int) matrix.cols();

    MatrixXd matrix_zeropad = MatrixXd::Zero(dim,ncols);
    MatrixXd matrix_reshaped = MatrixXd::Zero(nrows, ncols*dim);
    VectorXd vec = VectorXd::Zero(ncols*dim);

    MatrixXd matrix_increased = MatrixXd::Zero(nrows*dim, ncols*dim);

    for (int i = 0; i < nrows; ++i){
        matrix_zeropad << matrix.row(i),
                MatrixXd::Zero(dim-1,ncols);

        for (int k = 0; k < ncols; ++k){
            vec.segment(k*dim, dim) = matrix_zeropad.col(k);
        }
        matrix_reshaped.row(i) = vec;
    }

    MatrixXd matrix_kron = kroneckerProduct(matrix_reshaped, VectorXd::Ones(dim));

    // Shifting operation
    for (int k =0; k < dim*nrows; ++k){
        // remainder tells us the amount of spaces we need to shift the vector to the right
        int remainder = k % dim;
        if (remainder != 0){
            matrix_increased.row(k) = shift_vector(matrix_kron.row(k),remainder);
        }
        else matrix_increased.row(k) = matrix_kron.row(k);
    }

    return matrix_increased;
}

void BezierCurve::set_ctrl_pts(Eigen::VectorXd x) {
    // check if the dimension of the ctrl pts matches in dimensions
    assert(x.size() == _num_ctrl_pts);
    _ctrl_pts = x;
}

MatrixXd BezierCurve::get_mat_sample_poly(Eigen::VectorXd t_samples, int deg) {
    float T_final = _num_segments * _t_segment;
    float t;
    int segment_id;
    int prev_segment_id;
    int samples = 0;
    bool force_update = false;
    std::vector<double> tmp;
    std::vector<MatrixXd> Tau(_num_segments);

    for (int k = 0; k < t_samples.size(); ++k) {

        if (abs(t_samples(k) - T_final) < 1e-4){
            t = _t_segment;
            force_update = true;

            for (int n = 0; n < deg + 1; ++n)
                tmp.push_back(pow(t / _t_segment, n));
            samples++;
        }

        else {
            segment_id = floor(t_samples(k) / _t_segment);
            t = fmod(t_samples(k), _t_segment);
        }

        bool new_segment = segment_id != prev_segment_id;

        if (new_segment && k != 0 || force_update){
            Tau.at(prev_segment_id) = MapMatrixXd(tmp.data(),samples, deg+1);
            tmp.clear();
            samples = 0;
        }

        for (int n = 0; n < deg + 1; ++n)
            tmp.push_back(pow(t / _t_segment, n));

        prev_segment_id = segment_id;
        samples++;
    }

    return build_mat_sample_poly(Tau, t_samples.size(), deg);
}

MatrixXd BezierCurve::build_mat_sample_poly(std::vector<Eigen::MatrixXd> Tau,
                                            int num_samples, int deg) {
    int N = deg + 1;
    MatrixXd T_sample_poly = MatrixXd::Zero(3*num_samples, 3*N*_num_segments);
    std::vector<MatrixXd> Tau_d(_num_segments);

    int nrows;
    int ncols;
    int curr_row = 0;

    for (int i = 0; i < _num_segments; ++i){
        if (Tau.at(i).cols() > 0){
            Tau_d.at(i) = increase_matrix_dim(Tau.at(i), _dim);
            nrows = Tau_d.at(i).rows();
            ncols = Tau_d.at(i).cols();
            T_sample_poly.block(curr_row, i*ncols, nrows, ncols) = Tau_d.at(i);
            curr_row += nrows;
        }
    }

    return T_sample_poly;
}

MatrixXd BezierCurve::get_mat_input_sampling(Eigen::VectorXd t_samples) {
    // based on the t_samples, create the internal matrix that will sample the optimized vector
    _Gamma = get_mat_sample_poly(t_samples, _deg);
    _GammaBeta = _Gamma*_Beta;
    return _GammaBeta;
}

VectorXd BezierCurve::get_input_sequence(Eigen::VectorXd x){
    return _GammaBeta * x;
}

MatrixXd BezierCurve::get_mat_sumsqrd_der(Eigen::VectorXd weights){
    float mult;
    MatrixXd Q_sum = MatrixXd::Zero(_deg + 1, _deg + 1);
    MatrixXd Q = MatrixXd::Zero(_deg + 1, _deg + 1);
    MatrixXd Qd;
    float T = _t_segment;

    for (int r = 0; r <= _deg; ++r){
        for (int i = 0; i <= _deg; ++i){
            for (int k = 0; k <= _deg; ++k){
                if (i >= r && k >= r){
                    mult = 1.0;
                    for (int m = 0; m <= r - 1; ++m){
                        mult = mult * (i - m) * (k - m);
                    }
                    Q(i, k) = weights(r) * 2 * mult * pow(T, (i + k - 2*r + 1)) / (i + k - 2*r + 1);
                }
                else
                    Q(i, k) = 0;
            }
        }

        Q_sum += Q;
        Q = MatrixXd::Zero(_deg + 1, _deg + 1);
    }

    return augmented_form(Q);
}

MatrixXd BezierCurve::get_mat_energy_cost(Eigen::VectorXd weights) {
    MatrixXd Alpha = get_mat_sumsqrd_der(weights);
    MatrixXd H_energy = _Beta.transpose() * Alpha * _Beta;
    return  H_energy;
}

std::vector<MatrixXd> BezierCurve::derivate_ctrl_pts() {
    // Construct matrices that transform position ctrl pts to ctrl pts of the n-th derivative.

    std::vector<MatrixXd> T_ctrl_pts;
    std::vector<MatrixXd> der_mats;
    RowVector2d rows;
    rows << -1, 1;

    for (int c = _deg; c > 0; --c){
        MatrixXd aux = MatrixXd::Zero(c, c + 1);
        RowVector2d insert_rows = c*rows;

        for (int l = 0; l < c; ++l)
            aux.block(l, l, 1, 2) = insert_rows;

        der_mats.push_back(aux);
    }

    for (int k = 0; k < _deg; ++k){
        MatrixXd aux = der_mats.at(0);
        for (int k_aux = _deg - 1; k_aux > k; --k_aux)
            aux = der_mats.at(_deg - k_aux) * aux;

        cout << aux << endl;
        T_ctrl_pts.push_back(aux);
    }
    // Reverse the order for intuitive ordering - first element corresponds to velocity ctrl pts
    std::reverse(T_ctrl_pts.begin(), T_ctrl_pts.end());
    return T_ctrl_pts;
}

Constraint BezierCurve::limit_derivative(int n, Eigen::VectorXd t_samples,
                                         VectorXd min, VectorXd max) {

    // Get the appropriate matrix that transforms ctrl points to ctrl points of the n-th derivative
    MatrixXd M = _T_ctrl_pts.at(n - 1);

    // Expand matrix to 3D
    MatrixXd Sigma_n = augmented_form(M);

    MatrixXd Beta_n = bernstein_to_power(_deg - n);
    MatrixXd Tau_n = get_mat_sample_poly(t_samples, _deg - n);

    // Build 'A' matrix for the constraint
    MatrixXd A_in = Tau_n * Beta_n * Sigma_n;
    MatrixXd A_in_t(2 * A_in.rows(), A_in.cols());
    A_in_t << A_in, -A_in;

    // Build 'b' vector for the constraint
    VectorXd b_max = max.replicate(t_samples.size(), 1);
    VectorXd b_min = min.replicate(t_samples.size(), 1);
    VectorXd b_in(2 * b_max.size());
    b_in << b_max, b_min;

    // Assemble constraint struct
    Constraint constr = {A_in, b_in};
    return constr;
}

MatrixXd BezierCurve::augmented_form(Eigen::MatrixXd mat) {
    MatrixXd mat_d = increase_matrix_dim(mat, _dim);
    return kroneckerProduct(MatrixXd::Identity(_num_segments, _num_segments), mat_d);
}





