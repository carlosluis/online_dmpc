//
// Created by carlos on 19/03/19.
//

#ifndef ONLINE_PLANNING_GENERATOR_H
#define ONLINE_PLANNING_GENERATOR_H

#include "bezier.h"
#include "model.h"

struct MpcParams {
//    float h, Ts;
//    int k_hor;
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
    BezierCurve* _bezier;
    DoubleIntegrator* _model_pred;
    DoubleIntegrator* _model_exec;
    Eigen::MatrixXd _H_snap;
};

#endif //ONLINE_PLANNING_GENERATOR_H
