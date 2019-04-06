//
// Created by carlos on 05/04/19.
//

#ifndef ONLINE_PLANNING_AVOIDANCE_H
#define ONLINE_PLANNING_AVOIDANCE_H

#include "iostream"
#include "bezier.h"
#include "model.h"

struct EllipseParams {
    int order;
    float rmin;
    Eigen::VectorXd c;
};

struct Ellipse {
    int order;
    float rmin;
    Eigen::MatrixXd E1, E2;
};

struct CollisionConstraint {
    Constraint constraint;
    Eigen::VectorXd distance;
};

class BaseAvoider {
public:
    BaseAvoider(){};
    ~BaseAvoider(){};

    virtual Constraint get_coll_constr(const State3D& state, int agent_id) = 0;
};

class OndemandAvoider : public BaseAvoider {
public:
    OndemandAvoider (const std::vector<Eigen::MatrixXd>& horizon, const Eigen::MatrixXd& Phi_pos,
                     const Eigen::MatrixXd& A0_pos, const std::vector<EllipseParams>& p);

    Constraint get_coll_constr(const State3D& state, int agent_id);


private:

    // Vector with ellipse parameters for each agent
    std::vector<Ellipse> _ellipse;

    // Create references to matrices and value on the Generator object
    // The Avoider won't change any of these values
    const std::vector<Eigen::MatrixXd>& _horizon;
    const Eigen::MatrixXd& _Phi_pos;
    const Eigen::MatrixXd& _A0_pos;
    int _k_hor;
    int _N;
    int _dim;

    // Private methods
    std::vector<bool> check_collisions(int agent_id, int k);
    CollisionConstraint build_coll_constr(int agent_id, int k, const std::vector<bool>& neighbours,
                                          const State3D& state);
};

#endif //ONLINE_PLANNING_AVOIDANCE_H
