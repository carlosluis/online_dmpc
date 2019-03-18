#include <iostream>
#include <Eigen/Dense>
#include "bezier.h"

using namespace std;
using namespace Eigen;

int main() {
	cout << "Hello world!" << endl;
	BezierCurve test(5,3, 1.0,3);
    test.set_ctrl_pts(VectorXd::Zero(54));
//    cout << test.bernstein_to_power(5) << endl;
    VectorXd cr = VectorXd::Zero(5+1);
    cr(2) = .008;

    VectorXd t_samples = VectorXd::LinSpaced(16, 0, 3);

    cout << t_samples << endl << endl;

    Vector3d alim = 3*VectorXd::Ones(3);

    Constraint hola = test.limit_derivative(2, t_samples, -alim, alim);

    MatrixXd hola2 = test.get_mat_eq_constr(3);

    std::vector<MatrixXd> hola3 = test.get_vec_input_sampling(t_samples);

//    for (int k = 0; k < hola3.size(); ++k){
//        cout << hola3.at(k).block(0,0,15,15) << endl << endl;
//    }


	return 0;
}
