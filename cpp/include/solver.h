//
// Created by carlos on 06/04/19.
//

#ifndef ONLINE_PLANNING_SOLVER_H
#define ONLINE_PLANNING_SOLVER_H

#include "iostream"
#include "bezier.h"
#include "model.h"

struct QuadraticProblem {
    Eigen::MatrixXd H, Aeq, Ain;
    Eigen::VectorXd f, beq, bin;
};



#endif //ONLINE_PLANNING_SOLVER_H
