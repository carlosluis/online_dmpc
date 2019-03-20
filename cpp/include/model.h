//
// Created by carlos on 19/03/19.
//

#ifndef ONLINE_PLANNING_MODEL_H
#define ONLINE_PLANNING_MODEL_H

#include <Eigen/Dense>

// Define a struct for lifted matrices in position and velocity separately
struct StatePropagator {
    Eigen::MatrixXd pos;
    Eigen::MatrixXd vel;
};

class DoubleIntegrator {
public:
    DoubleIntegrator(){};
    ~DoubleIntegrator(){};

    StatePropagator get_lambda(int K);
    StatePropagator get_A0(int K);

protected:
    Eigen::MatrixXd _model_A;
    Eigen::MatrixXd _model_B;
    int _dim;
};

class DoubleIntegrator3D : public DoubleIntegrator {
public:
    // Define parameters to create a 3D double integrator model
    // We assume identical dynamics for X and Y, based on quadcopter experiments
    struct Params {
        float zeta_xy,  tau_xy, zeta_z, tau_z;
    };
    DoubleIntegrator3D (const float& time_step, const DoubleIntegrator3D::Params& p);
    ~DoubleIntegrator3D(){};
};

class DoubleIntegrator2D : public DoubleIntegrator {
public:
    // To be used for ground robots or quadcopters flying in the plane
    struct Params {
        float zeta_x,  tau_x, zeta_y, tau_y;
    };

    DoubleIntegrator2D (const float& time_step, const DoubleIntegrator2D::Params& p);
    ~DoubleIntegrator2D(){};
};

#endif //ONLINE_PLANNING_MODEL_H
