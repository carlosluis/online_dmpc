//
// Created by carlos on 19/03/19.
//

#ifndef ONLINE_PLANNING_GENERATOR_H
#define ONLINE_PLANNING_GENERATOR_H

#include "bezier.h"
#include "model.h"

struct PhysLimits {
    Eigen::VectorXd pmax;
    Eigen::VectorXd pmin;
    Eigen::VectorXd amax;
    Eigen::VectorXd amin;
};

struct MpcParams {
    float h, Ts;
    int k_hor;
    const PhysLimits& limits;
    Eigen::VectorXd energy_weights;
};

class Generator {
public:
    struct Params {
        const BezierCurve::Params& bezier_params;
        const DoubleIntegrator3D::Params& model_params;
        const MpcParams& mpc_params;
    };

    Generator(const Generator::Params& p);
    ~Generator(){};

private:
    float _h;
    float _Ts;
    float _k_hor;
    int _l;

    BezierCurve _bezier;
    DoubleIntegrator3D _model_pred;
    DoubleIntegrator3D _model_exec;
    Eigen::MatrixXd _H_energy;
    Constraint _ineq;
    Constraint _eq;
    StatePropagator _Phi_pred;
    StatePropagator _Phi_exec;

    // Methods
    void build_ineq_constr(const PhysLimits& lim, Constraint* ineq);
};

#endif //ONLINE_PLANNING_GENERATOR_H
