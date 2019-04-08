//
// Created by carlos on 19/03/19.
//

#ifndef ONLINE_PLANNING_GENERATOR_H
#define ONLINE_PLANNING_GENERATOR_H

#include "bezier.h"
#include "model.h"
#include "avoidance.h"
#include "solver.h"
#include <thread>

struct PhysLimits {
    Eigen::VectorXd pmax, pmin, amax, amin;
};

struct TuningParams {
    float s_free, s_obs, s_repel;
    int spd_f, spd_o, spd_r;
    double lin_coll, quad_coll;
    Eigen::VectorXd energy_weights;
};

struct MpcParams {
    const float& h, Ts;
    const int& k_hor;
    const TuningParams& tuning;
    const PhysLimits& limits;
};

class Generator {
public:
    struct Params {
        const BezierCurve::Params& bezier_params;
        const DoubleIntegrator3D::Params& model_params;
        const std::vector<EllipseParams>& ellipse;
        const MpcParams& mpc_params;
        const Eigen::MatrixXd& po, pf;
    };

    Generator(const Generator::Params& p);
    ~Generator(){};

    // Public methods
//    void get_next_inputs();
    void get_next_inputs(const std::vector<State3D>& curr_states);
//    std::vector<Eigen::MatrixXd> get_next_inputs(const std::vector<State3D>& curr_states);


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
    int _deg_poly;

    double _lin_coll;
    double _quad_coll;

    Eigen::MatrixXd _po;
    Eigen::MatrixXd _pf;

    float _max_cost;
    float _min_cost;

    BezierCurve _bezier;
    DoubleIntegrator3D _model_pred;
    DoubleIntegrator3D _model_exec;
    Eigen::MatrixXd _H_energy;
    InequalityConstraint _ineq;
    Constraint _eq;

    // Variables related to multithreading and clustering
    const int _max_clusters;
    std::vector<std::thread> _t;
    std::vector<std::vector<int>> _cluster;

    // Matrix that samples the input for t in [t0, t0 + h], with sampling time ts
    Eigen::MatrixXd _Rho_position;

    // State propagators using the model
    StatePropagator _Lambda_pred;
    StatePropagator _Lambda_exec;
    StatePropagator _A0_pred;
    StatePropagator _A0_exec;
    StatePropagator _Phi_pred;
    StatePropagator _Phi_exec;

    // Collision avoidance module
    std::unique_ptr<BaseAvoider> _avoider;

    // QP solver moduler

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
    void init_generator();
    InequalityConstraint build_ineq_constr(const PhysLimits& lim);
    void set_error_penalty_mats(const TuningParams& p, const Eigen::MatrixXd& pf);
    void init_clusters();
    Eigen::MatrixXd get_init_ref(const State3D& state, const Eigen::MatrixXd& ref);
    void solve_cluster(const std::vector<State3D>& curr_states,
                       const std::vector<int>& agents);

    QuadraticProblem build_qp(const Constraint& collision, const State3D& state,
                              int agent_id);

};

#endif //ONLINE_PLANNING_GENERATOR_H
