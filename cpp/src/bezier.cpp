#include <iostream>
#include "bezier.h"

using namespace std;
using namespace Eigen;

BezierCurve::BezierCurve(int deg, int l, int dim){
	cout << "Creating a Bezier curve with " << l << " segments of degree "
         << deg << " in " << dim << " dimensions" << endl;

    // Assign to private variables
    _deg = deg;
    _l = l;
    _dim = dim;

    cout << "Initializing the vector of control points as zero" << endl;
    _num_ctrl_pts = dim*l*(deg+1);
    _ctrl_pts = VectorXd::Zero(_num_ctrl_pts);
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
    MatrixXd aux = MatrixXd::Zero(dim,ncols);
    VectorXd aux2 = VectorXd::Zero(ncols*dim);
    MatrixXd aug_array = MatrixXd::Zero(nrows, ncols*dim);
    MatrixXd aug_array2 = MatrixXd::Zero(nrows*dim, ncols*dim);
    MatrixXd increased_matrix = MatrixXd::Zero(nrows*dim, ncols*dim);

    for (int i = 0; i < nrows; ++i){
        aux << matrix.row(i),
                MatrixXd::Zero(dim-1,ncols);

        for (int k = 0; k < ncols; ++k){
            aux2.segment(k*dim, dim) = aux.col(k);
        }
        aug_array.row(i) = aux2;
    }

    aug_array2 = kroneckerProduct(aug_array, VectorXd::Ones(dim));
    cout << aug_array2 << endl << endl;

    // Shifting operation
    for (int k =0; k < dim*ncols; ++k){
        // remainder tells us the amount of spaces we need to shift the vector to the right
        int remainder = k % dim;
        if (remainder != 0){
            increased_matrix.row(k) = shift_vector(aug_array2.row(k),remainder);
        }
        else increased_matrix.row(k) = aug_array2.row(k);
    }

    cout << increased_matrix << endl;
    return increased_matrix;
}

void BezierCurve::set_ctrl_pts(Eigen::VectorXd x) {
    // check if the dimension of the ctrl pts matches in dimensions
    assert(x.size() == _num_ctrl_pts);
    _ctrl_pts = x;
}

Eigen::VectorXd BezierCurve::get_sample(float hola){

}



