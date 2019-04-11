//
// Created by carlos on 09/04/19.
//

#ifndef ONLINE_PLANNING_SIMULATOR_H
#define ONLINE_PLANNING_SIMULATOR_H

#include "generator.h"
#include "json.hpp"
#include <random>
#include <fstream>

class Simulator {
public:
    struct Params {
        Generator::Params generator_params;
        float position_noise_std;
        float velocity_noise_std;
    };

//    explicit Simulator(const Simulator::Params& p);
    explicit Simulator(std::ifstream& config_file);
    ~Simulator(){};


    void run(int duration);
    void saveDataToFile(char const* pathAndName);

private:

    // Members
    BezierCurve::Params _bezier_params;
    DoubleIntegrator3D::Params _model_params;
    MpcParams _mpc_params;
    std::unique_ptr<Generator> _generator;
    std::unique_ptr<DoubleIntegrator3D> _sim_model;
    std::vector<Ellipse> _ellipses;
    std::vector<Eigen::MatrixXd> _inputs;
    std::vector<State3D> _current_states;
    std::vector<Eigen::MatrixXd> _trajectories;
    float _pos_std;
    float _vel_std;
    Eigen::MatrixXd _po;
    Eigen::MatrixXd _pf;
    Eigen::VectorXd _pmin;
    Eigen::VectorXd _pmax;
    int _Ncmd;
    int _N;
    float _h;
    float _Ts;

    // Methods
    Generator::Params parseJSON(std::ifstream& config_file);
    State3D addRandomNoise(const State3D& states);
    bool collisionCheck(const std::vector<Eigen::MatrixXd>& trajectories);
    bool goalCheck(const std::vector<State3D>& states);
    Eigen::MatrixXd generateRandomPoints(int N, const Eigen::Vector3d &pmin,
                                         const Eigen::Vector3d &pmax, float rmin);
};

#endif //ONLINE_PLANNING_SIMULATOR_H
