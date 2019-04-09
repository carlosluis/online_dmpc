//
// Created by carlos on 06/04/19.
//

#include "solver.h"

using namespace Eigen;
using namespace std;

bool Quadprog::solveQP(const QuadraticProblem& problem) {
    int num_inequalities = problem.bin_full.size();
    int num_equalities = problem.beq.size();
    int num_vars = problem.f.size();
    MatrixXd H = 2 * problem.H;
    MatrixXd H_sym = 0.5 * (H + H.transpose());

    QuadProgDense qp(num_vars, num_equalities, num_inequalities);
    qp.solve(H_sym, problem.f, problem.Aeq, problem.beq, problem.Ain_full, problem.bin_full, true);
    VectorXd hola = qp.result();
    cout << "Result: " << qp.fail() << endl;
    std::cout << hola << std::endl;
    return static_cast<bool>(qp.fail());
}

bool OOQP::solveQP(const QuadraticProblem& problem) {
    int num_vars = problem.f.size();

    bool status = ooqpei::OoqpEigenInterface::solve(2 * problem.H.sparseView(), problem.f,
                                                    problem.Aeq.sparseView(), problem.beq,
                                                    problem.Ain.sparseView(), problem.bin_lower,
                                                    problem.bin_upper,  _solution);
    return status;
}