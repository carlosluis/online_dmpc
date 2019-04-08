//
// Created by carlos on 06/04/19.
//

#ifndef ONLINE_PLANNING_SOLVER_H
#define ONLINE_PLANNING_SOLVER_H

#include "iostream"
#include "bezier.h"
#include "model.h"
#include "eigen-quadprog/src/QuadProg.h"
#include "ooqp_eigen_interface/OoqpEigenInterface.hpp"

struct QuadraticProblem {
    Eigen::MatrixXd H, Aeq, Ain_full, Ain;
    Eigen::VectorXd f, beq, bin_full, bin_lower, bin_upper;
};

class BaseSolver {
public:
    BaseSolver(){};
    virtual ~BaseSolver(){};

    virtual bool solveQP(const QuadraticProblem& problem) = 0;
    virtual Eigen::VectorXd getSolution() = 0;
};

class Quadprog : public BaseSolver {
public:
    Quadprog(){};
    ~Quadprog(){};

    bool solveQP(const QuadraticProblem& problem);
    Eigen::VectorXd getSolution(){return _solution;};

private:
    Eigen::VectorXd _solution;
};

class OOQP : public BaseSolver {
public:
    OOQP(){};
    ~OOQP(){};

    bool solveQP(const QuadraticProblem& problem);
    Eigen::VectorXd getSolution(){return _solution;};

private:
    Eigen::VectorXd _solution;
};

#endif //ONLINE_PLANNING_SOLVER_H
