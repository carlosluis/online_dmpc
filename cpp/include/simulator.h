//
// Created by carlos on 09/04/19.
//

#ifndef ONLINE_PLANNING_SIMULATOR_H
#define ONLINE_PLANNING_SIMULATOR_H

#include "generator.h"

class Simulator {
public:
    struct Params {
        const Generator::Params& generator_params;
        const float position_noise_std;
        const float velocity_noise_std;
    };

    explicit Simulator(const Simulator::Params& p);
    ~Simulator(){};

    void run(int duration);

private:

    // Members
    Generator _generator;
    DoubleIntegrator3D _sim_model;
    float _pos_std;
    float _vel_std;

    // Methods
    State3D addRandomNoise(const State3D& states);
    Eigen::MatrixXd generateRandomPoints(int N, const Eigen::Vector3d& pmin,
                                         const Eigen::Vector3d& pmax, float rmin);
};

#endif //ONLINE_PLANNING_SIMULATOR_H
