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
    PhysLimits* limits;
    Eigen::VectorXd energy_weights;
};

class Generator {
public:
    struct Params {
        BezierCurve::Params* bezier_params;
        DoubleIntegrator3D::Params* model_params;
        MpcParams* mpc_params;
    };

    Generator(Generator::Params* p);
    ~Generator(){};

private:
    float _h;
    float _Ts;
    float _k_hor;
    int _l;

    BezierCurve* _bezier;
    DoubleIntegrator* _model_pred;
    DoubleIntegrator* _model_exec;
    Eigen::MatrixXd _H_snap;

    // Methods
    Constraint build_ineq_constr(PhysLimits* lim);
};

#endif //ONLINE_PLANNING_GENERATOR_H
