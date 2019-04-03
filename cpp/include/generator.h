//
// Created by carlos on 19/03/19.
//

#ifndef ONLINE_PLANNING_GENERATOR_H
#define ONLINE_PLANNING_GENERATOR_H

#include "bezier.h"
#include "model.h"
#include <thread>

struct PhysLimits {
    Eigen::VectorXd pmax;
    Eigen::VectorXd pmin;
    Eigen::VectorXd amax;
    Eigen::VectorXd amin;
};

struct TuningParams {
    float s_free, s_obs, s_repel;
    int spd_f, spd_o, spd_r;
    double lin_coll, quad_coll;
    Eigen::VectorXd energy_weights;
};

struct CollisionParams {
    int order;
    float rmin;
    float c;
};

struct MpcParams {
    float h, Ts;
    int k_hor;
    const TuningParams& tuning;
    const PhysLimits& limits;
    const CollisionParams coll;
};

class Generator {
public:
    struct Params {
        const BezierCurve::Params& bezier_params;
        const DoubleIntegrator3D::Params& model_params;
        const MpcParams& mpc_params;
        const Eigen::MatrixXd& po;
        const Eigen::MatrixXd& pf;
    };

    Generator(const Generator::Params& p);
    ~Generator(){};

    // Public methods
    std::vector<Eigen::MatrixXd> get_next_inputs(const std::vector<State3D>& curr_states);

private:
    float _h;
    float _Ts;
    float _k_hor;
    int _l;
    int _num_ctrl_pts;
    int _dim;
    int _Ncmd;
    int _N;
    int _d;

    Eigen::MatrixXd _po;
    Eigen::MatrixXd _pf;

    BezierCurve _bezier;
    DoubleIntegrator3D _model_pred;
    DoubleIntegrator3D _model_exec;
    Eigen::MatrixXd _H_energy;
    Constraint _ineq;
    Constraint _eq;

    // Variables related to multithreading and clustering
    const int _max_clusters;
    std::vector<std::thread> _t;
    std::vector<std::vector<int>> _cluster;

    // State propagators using the model
    StatePropagator _Lambda_pred;
    StatePropagator _Lambda_exec;
    StatePropagator _A0_pred;
    StatePropagator _A0_exec;
    StatePropagator _Phi_pred;
    StatePropagator _Phi_exec;

    // Matrices and vectors to minimize goal error
    Eigen::MatrixXd _H_f;
    Eigen::MatrixXd _H_o;
    Eigen::MatrixXd _H_r;
    Eigen::MatrixXd _Hlin_f;
    Eigen::MatrixXd _Hlin_o;
    Eigen::MatrixXd _Hlin_r;
    std::vector<Eigen::RowVectorXd> _fpf_free;
    std::vector<Eigen::RowVectorXd> _fpf_obs;

    // Horizon variables, one for solving and one for updating
    std::vector<Eigen::MatrixXd> _newhorizon;
    std::vector<Eigen::MatrixXd> _oldhorizon;
    std::vector<Eigen::MatrixXd> _x0_ref;

    // Methods
    Constraint build_ineq_constr(const PhysLimits& lim);
    void set_error_penalty_mats(const TuningParams& p, const Eigen::MatrixXd& pf);
    void init_clusters();
    void init_generator();

};

#endif //ONLINE_PLANNING_GENERATOR_H
