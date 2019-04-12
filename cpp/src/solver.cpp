//
// Created by carlos on 06/04/19.
//

#include "solver.h"

using namespace Eigen;
using namespace std;

typedef Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,
        Eigen::Dynamic, Eigen::ColMajor>> MapMatrixCol;

typedef Matrix<double, Dynamic, Dynamic, RowMajor> MatrixXdRow;

bool Quadprog::solveQP(const QuadraticProblem& problem) {
    // DOESN'T WORK PROPERLY!
    int num_inequalities = problem.bin_full.size();
    int num_equalities = problem.beq.size();
    int num_vars = problem.f.size();
    MatrixXd H = 2 * problem.H;
    MatrixXd H_sym = 0.5 * (H + H.transpose());

    QuadProgDense qp(num_vars, num_equalities, num_inequalities);
    qp.solve(H_sym, problem.f, problem.Aeq, problem.beq, problem.Ain_full, problem.bin_full, false);
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

bool QpOASES::solveQP(const QuadraticProblem &problem) {
    int num_inequalities = problem.bin_lower.size();
    int num_equalities = problem.beq.size();
    int num_vars = problem.f.size();
    int num_constraints = num_inequalities + num_equalities;
    _solution = VectorXd::Zero(num_vars);

    // Since qpOASES only takes inequality constraints with lower and upper bounds, we will
    // transform the equality constraints to be beq <= Aeq <= beq
    VectorXd bin_lower = VectorXd::Zero(num_constraints);
    bin_lower << problem.bin_lower, problem.beq;
    VectorXd bin_upper = VectorXd::Zero(num_constraints);
    bin_upper << problem.bin_upper, problem.beq;

    MatrixXd Ain = MatrixXd::Zero(num_constraints, num_vars);
    Ain << problem.Ain, problem.Aeq;

    // Matrices Ain and H need to be in Rowmajor form for qpOASES
    MatrixXd H = 2 * problem.H;
    MatrixXdRow H_row = MapMatrixCol(H.data(), num_vars, num_vars);
    MatrixXdRow Ain_row = MapMatrixCol(Ain.data(), num_constraints, num_vars);

    // Start building the qpOASES object and problem to be solved
    qpOASES::QProblem qp(num_vars, num_constraints, qpOASES::HST_SEMIDEF);
    qpOASES::Options solver_options;
    solver_options.setToMPC();
    solver_options.printLevel = qpOASES::PL_LOW;
    qp.setOptions(solver_options);
    qpOASES::int_t nWSR = 1000;
    qpOASES::real_t cpu_time = 0.01;
    qpOASES::returnValue status = qp.init(H_row.data(), problem.f.data(), Ain_row.data(),
                                          NULL, NULL, bin_lower.data(), bin_upper.data(),
                                          nWSR, &cpu_time);

    qpOASES::int_t exit_flag = qpOASES::getSimpleStatus(status);
    bool solvedOK = exit_flag == 0;

    if (solvedOK)
        qp.getPrimalSolution(_solution.data());

    return solvedOK;
}