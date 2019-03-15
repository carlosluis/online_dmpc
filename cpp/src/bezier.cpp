#include <iostream>
#include "bezier.h"

using namespace std;
using namespace Eigen;

BezierCurve::BezierCurve(int deg, int num_segments, int dim){
	cout << "Creating a Bezier curve with " << num_segments << " segments of degree "
         << deg << " in " << dim << " dimensions" << endl;

    // Assign to private variables
    _deg = deg;
    _num_segments = num_segments;
    _dim = dim;

    cout << "Initializing the vector of control points as zero" << endl;
    _num_ctrl_pts = dim*num_segments*(deg+1);
    _ctrl_pts = VectorXd::Zero(_num_ctrl_pts);

    // Build the matrices defined by the bezier curve parameters
    _Beta = bernstein_to_power(_deg);
}

MatrixXd BezierCurve::bernstein_to_power(int deg) {
    MatrixXd delta = MatrixXd::Zero(deg + 1, deg + 1);
    for (int k = 0; k < deg + 1; ++k){
        for (int i = k; i <= deg; ++i){
            delta(i, k) = pow(-1, i - k)*nchoosek(deg, i)*nchoosek(i, k);
        }
    }
    return increase_matrix_dim(delta,_dim);;
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

MatrixXd BezierCurve::get_mat_sample_poly(float segment_duration,
                                          Eigen::VectorXd t_samples,
                                          int deg) {
    float T_final = _num_segments * segment_duration;
    float t;
    int segment_id;
    int prev_segment_id;
    int samples = 0;
    bool force_update = false;
    std::vector<double> tmp;
    std::vector<MatrixXd> Tau(_num_segments);

    for (int k = 0; k < t_samples.size(); ++k) {

        if (abs(t_samples(k) - T_final) < 1e-4){
            t = segment_duration;
            force_update = true;

            for (int n = 0; n < deg + 1; ++n)
                tmp.push_back(pow(t / segment_duration, n));
            samples++;
        }

        else {
            segment_id = floor(t_samples(k) / segment_duration);
            t = fmod(t_samples(k), segment_duration);
        }

        bool new_segment = segment_id != prev_segment_id;

        if (new_segment && k != 0 || force_update){
            Tau.at(prev_segment_id) = MapMatrixXd(tmp.data(),samples, deg+1);
            tmp.clear();
            samples = 0;
        }

        for (int n = 0; n < deg + 1; ++n)
            tmp.push_back(pow(t / segment_duration, n));

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

MatrixXd BezierCurve::get_mat_input_sampling(float segment_duration,
                                             Eigen::VectorXd t_samples) {
    // based on the t_samples, create the internal matrix that will sample the optimized vector
    _Gamma = get_mat_sample_poly(segment_duration, t_samples, _deg);
    _GammaBeta = _Gamma*_Beta;
    return _GammaBeta;
}

VectorXd BezierCurve::get_input_sequence(Eigen::VectorXd x){
    return _GammaBeta * x;
}

MatrixXd BezierCurve::get_mat_sumsqrd_der(float T, Eigen::VectorXd weights){
    float mult;
    MatrixXd Q_sum = MatrixXd::Zero(_deg + 1, _deg + 1);
    MatrixXd Q = MatrixXd::Zero(_deg + 1, _deg + 1);
    MatrixXd Qd;

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

    Qd = increase_matrix_dim(Q_sum, _dim);
    return kroneckerProduct(MatrixXd::Identity(_num_segments, _num_segments), Qd);
}

MatrixXd BezierCurve::get_mat_energy_cost(float segment_duration, Eigen::VectorXd weights) {
    MatrixXd Alpha = get_mat_sumsqrd_der(segment_duration, weights);
    MatrixXd H_energy = _Beta.transpose() * Alpha * _Beta;
    return  H_energy;
}





