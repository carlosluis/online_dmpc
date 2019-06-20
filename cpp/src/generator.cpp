//
// Created by carlos on 19/03/19.
//
#include "iostream"
#include "generator.h"

using namespace std;
using namespace Eigen;
using namespace std::chrono;

Generator::Generator(const Generator::Params& p) :
    _bezier(p.bezier_params),
    _model_pred(p.mpc_params.h, p.model_params),
    _model_exec(p.mpc_params.Ts, p.model_params),
    _max_clusters(1),
    _max_cost(0.08),
    _min_cost(-0.01)
{
    // Unpack params from struct into private members
    _h = p.mpc_params.h;
    _Ts = p.mpc_params.Ts;
    _k_hor = p.mpc_params.k_hor;
    _dim = p.bezier_params.dim;
    _l = p.bezier_params.num_segments;
    _d = p.bezier_params.deg;
    _deg_poly = p.bezier_params.deg_poly;
    _po = p.po;
    _pf = p.pf;
    _solver_name = p.solver_name;
    _N = _po.cols();
    _Ncmd = _pf.cols();
    _fpf_free.reserve(_Ncmd);
    _fpf_obs.reserve(_Ncmd);
    _lin_coll = p.mpc_params.tuning.lin_coll;
    _quad_coll = p.mpc_params.tuning.quad_coll;

    // Define two time bases:
    // 1) h_samples is the time base for the planning of trajectories
    // 2) ts_samples is a finer base for the execution in-between planning updates
    VectorXd ts_samples = VectorXd::LinSpaced(_h / _Ts, 0, _h - _Ts);
    VectorXd h_samples = VectorXd::LinSpaced(_k_hor, 0, (_k_hor - 1) * _h);

    // Energy cost function associated Hessian matrix
    _H_energy = _bezier.getMatrixEnergyCost(p.mpc_params.tuning.energy_weights);

    _num_ctrl_pts = _H_energy.rows();

    // Get the matrices that sample the optimization result according to the 2 time bases
    _Rho_h = _bezier.getMatrixInputDerivativeSampling(h_samples);
    std::vector<MatrixXd> Rho_ts = _bezier.getMatrixInputDerivativeSampling(ts_samples);
    _Rho_pos_ts = Rho_ts[0];

    // Get the matrices that propagate the input (i.e., optimization result) using the model
    _Lambda_pred = _model_pred.get_lambda(_k_hor);
    _A0_pred = _model_pred.get_A0(_k_hor);
    _Lambda_exec = _model_exec.get_lambda(_h / _Ts);
    _A0_exec = _model_exec.get_A0(_h / _Ts);

    // Create new matrices that map bezier curve control points to states of the robot
    _Phi_pred.pos = _Lambda_pred.pos * _Rho_h.at(0);
    _Phi_pred.vel = _Lambda_pred.vel * _Rho_h.at(0);
    _Phi_exec.pos = _Lambda_exec.pos * Rho_ts.at(0);
    _Phi_exec.vel = _Lambda_exec.vel * Rho_ts.at(0);

    // Build the inequality constraints regarding the physical limits of the robot
    _ineq = buildInequalityConstraint(p.mpc_params.limits);

    // Build Aeq matrix for the equality constraints
    _eq.A = _bezier.getMatrixEqualityConstraint(p.bezier_params.deg_poly);

    // Set matrices to minimize goal error
    setErrorPenaltyMatrices(p.mpc_params.tuning, _pf);

    // Create threads and clusters to solve in parallel
    initClusters();

    // Initialize several variables used to generate the trajectories
    initGenerator();

    // Create the avoider object for collision avoidance
    _avoider = std::make_unique<OndemandAvoider>(_oldhorizon, _Phi_pred.pos,
                                                 _A0_pred.pos, p.ellipse);
}

void Generator::setErrorPenaltyMatrices(const TuningParams &p, const Eigen::MatrixXd &pf){

    // Case 1: no collisions in the horizon - go fast
    MatrixXd S_free = MatrixXd::Zero(_dim * _k_hor, _dim * _k_hor);

    Ref<MatrixXd> block_f = S_free.block(_dim * (_k_hor - p.spd_f), _dim * (_k_hor - p.spd_f),
                                         _dim * p.spd_f, _dim * p.spd_f);

    block_f = p.s_free * MatrixXd::Identity(_dim * p.spd_f, _dim * p.spd_f);

    // Case 2: collisions in the horizon - go slower
    MatrixXd S_obs = MatrixXd::Zero(_dim * _k_hor, _dim * _k_hor);

    Ref<MatrixXd> block_o = S_obs.block(_dim * (_k_hor - p.spd_o), _dim * (_k_hor - p.spd_o),
                                         _dim * p.spd_o, _dim * p.spd_o);

    block_o = p.s_obs * MatrixXd::Identity(_dim * p.spd_o, _dim * p.spd_o);

    // Case 3: colliding scenario, change setpoints and move fast towards it

    MatrixXd S_repel = MatrixXd::Zero(_dim * _k_hor, _dim * _k_hor);

    Ref<MatrixXd> block_r = S_repel.block(0, 0, _dim * p.spd_r, _dim * p.spd_r);

    block_r = p.s_repel * MatrixXd::Identity(_dim * p.spd_r, _dim * p.spd_r);

    // Assemble matrices
    MatrixXd H_free = _Phi_pred.pos.transpose() * S_free * _Phi_pred.pos;
    MatrixXd H_obs = _Phi_pred.pos.transpose() * S_obs * _Phi_pred.pos;
    MatrixXd H_repel = _Phi_pred.pos.transpose() * S_repel * _Phi_pred.pos;

    // Get the final Hessians by multiplying by the energy consumption hessian
    _H_f = H_free + _H_energy;
    _H_o = H_obs + _H_energy;
    _H_r = H_repel + _H_energy;

    // Build the matrices that multiply the initial condition to obtain the linear term of the cost
    _Hlin_f = _A0_pred.pos.transpose() * S_free * _Phi_pred.pos;
    _Hlin_o = _A0_pred.pos.transpose() * S_obs * _Phi_pred.pos;
    _Hlin_r = _A0_pred.pos.transpose() * S_repel * _Phi_pred.pos;

    // Build the constant part of the linear term of the cost function
    for (int i = 0; i < _Ncmd; i++) {
        VectorXd pfi = _pf.col(i);
        _fpf_free.push_back(pfi.replicate(_k_hor, 1).transpose() * S_free * _Phi_pred.pos);
        _fpf_obs.push_back(pfi.replicate(_k_hor, 1).transpose() * S_obs * _Phi_pred.pos);
    }

}

InequalityConstraint Generator::buildInequalityConstraint(const PhysicalLimits &limits){
    // Get position constraint
    VectorXd pos_samples = VectorXd::LinSpaced((_k_hor)/3, 1.0, floor(_k_hor / 3) * _h * _l );
//    cout << pos_samples << endl;
    InequalityConstraint lim_pos = _bezier.limitDerivative(0, pos_samples, limits.pmin, limits.pmax);

    InequalityConstraint lim_pos2 = _bezier.limitControlPoints(0, limits.pmin, limits.pmax);

    // Get acceleration constraint
    VectorXd acc_samples = VectorXd::LinSpaced((_k_hor)/2, _h, (_k_hor - 1) * _h);
    InequalityConstraint lim_acc = _bezier.limitDerivative(2, acc_samples, limits.amin, limits.amax);

    // Assemble the constraints
    InequalityConstraint ineq;
    ineq.A_full = MatrixXd::Zero(lim_pos.A_full.rows() + lim_acc.A_full.rows(), lim_pos.A_full.cols());
    ineq.b_full = VectorXd::Zero(lim_pos.b_full.size() + lim_acc.b_full.size());
    ineq.A_full << lim_acc.A_full, lim_pos.A_full;
    ineq.b_full << lim_acc.b_full, lim_pos.b_full;

    ineq.A = MatrixXd::Zero(lim_pos.A.rows() + lim_acc.A.rows(), lim_pos.A.cols());
    ineq.lower_bound = VectorXd::Zero(lim_pos.lower_bound.size() + lim_acc.lower_bound.size());
    ineq.upper_bound = VectorXd::Zero(lim_pos.upper_bound.size() + lim_acc.upper_bound.size());

    ineq.A << lim_acc.A, lim_pos.A;
    ineq.lower_bound << lim_acc.lower_bound, lim_pos.lower_bound;
    ineq.upper_bound << lim_acc.upper_bound, lim_pos.upper_bound;

    return ineq;
}

void Generator::initClusters() {
    int n_clusters = min(_Ncmd, _max_clusters);
    _t.resize(n_clusters);
    _cluster.resize(n_clusters);
    int agents_per_cluster = _Ncmd / n_clusters;
    int residue = _Ncmd % n_clusters;
    int cluster_capacity;
    int N_index = 0;

    for (int i = 0; i < n_clusters; i++) {
        if (residue != 0) {
            cluster_capacity = agents_per_cluster + 1;
            residue--;
        }
        else cluster_capacity = agents_per_cluster;

        int curr_index = N_index;
        for (int j = curr_index; j < curr_index + cluster_capacity; j++) {
            _cluster[i].push_back(j);
            N_index = j + 1;
        }
    }
}

void Generator::initGenerator() {
    _x0_ref.reserve(_N);
    _next_inputs.reserve(_N);
    _newhorizon.reserve(_N);
    _oldhorizon.reserve(_N);
    MatrixXd pos_aux = MatrixXd::Zero(_dim, _k_hor);
    MatrixXd init_ref = MatrixXd::Zero(_dim, _d + 1);

    for (int i = 0; i < _N; i++) {
        VectorXd poi = _po.col(i);
        VectorXd voi = 0.001 * VectorXd::Ones(_dim);
        init_ref << poi, voi, MatrixXd::Zero(_dim, _d - 1);
        _x0_ref.push_back(init_ref);
        _oldhorizon.push_back(poi.replicate(1, _k_hor));
        _next_inputs.push_back(poi.replicate(1, _h / _Ts));
    }
    _newhorizon = _oldhorizon;
}

std::vector<MatrixXd> Generator::getNextInputs(const vector<State3D>& current_states) {

    // Launch all the threads to get the next set of inputs
    for (int i = 0; i < _cluster.size(); i++)
        _t[i] = std::thread(&Generator::solveCluster, this, ref(current_states), _cluster[i]);

    // Wait for all agents to solve the optimization
    for (int i = 0; i < _cluster.size(); ++i)
        _t[i].join();

    // Update the old horizon with the newhorizon for the next optimization
    _oldhorizon = _newhorizon;

    return _next_inputs;
}

void Generator::solveCluster(const std::vector<State3D> &current_states,
                              const std::vector<int> &agents) {

    VectorXd err_pos, cost, denominator;
    for (int i = agents.front(); i <= agents.back(); i++) {
        // Pick initial condition for the reference
        _x0_ref[i] = getInitialReference(current_states[i], _x0_ref[i]);

        // Build collision constraint
//        high_resolution_clock::time_point t1 = high_resolution_clock::now();
        Constraint collision = _avoider->getCollisionConstraint(current_states[i], i);
//        high_resolution_clock::time_point t2 = high_resolution_clock::now();
//        auto duration = duration_cast<microseconds>( t2 - t1 ).count();
//        cout << "Time building constraint = "
//             << duration/1000.0 << "ms" << endl << endl;

        // We need to increase the size of other matrices if collision constraints are included
//        t1 = high_resolution_clock::now();
        QuadraticProblem problem = buildQP(collision, current_states[i], i);
//        t2 = high_resolution_clock::now();
//        duration = duration_cast<microseconds>( t2 - t1 ).count();
//        cout << "Time assembling QP = "
//             << duration/1000.0 << "ms" << endl << endl;

        // Create a new solver pointer of base class
        unique_ptr<BaseSolver> solver;
        if (_solver_name == kOoqp)
            solver = make_unique<OOQP>();
        else if (_solver_name == kQpoases)
            solver = make_unique<QpOASES>();

        // Solve QP and get solution vector x
//        t1 = high_resolution_clock::now();
        bool solvedOK = solver->solveQP(problem);
//        t2 = high_resolution_clock::now();
//        duration = duration_cast<microseconds>( t2 - t1 ).count();
//        cout << "Time solving QP = "
//             << duration/1000.0 << "ms" << endl << endl;

        if (solvedOK){
            // Extract solution
            VectorXd solution = solver->getSolution();

            // Strip away the slack variables from the solution
            VectorXd u = solution.head(_num_ctrl_pts);

            // Update horizon for the agent
//            t1 = high_resolution_clock::now();
            _newhorizon[i] = updateHorizon(u, current_states[i]);
//            t2 = high_resolution_clock::now();
//            duration = duration_cast<microseconds>( t2 - t1 ).count();
//            cout << "Time updating horizon = "
//                 << duration/1000.0 << "ms" << endl << endl;

            // Update the initial condition for the next optimization
//            t1 = high_resolution_clock::now();
            _x0_ref[i] = updateInitialReference(u);
//            t2 = high_resolution_clock::now();
//            duration = duration_cast<microseconds>( t2 - t1 ).count();
//            cout << "Time updating initial point of reference = "
//                 << duration/1000.0 << "ms" << endl << endl;

            // Update the next input sequence for agent i
//            t1 = high_resolution_clock::now();
            _next_inputs[i] = vec2mat(_Rho_pos_ts * u, _dim);
//            t2 = high_resolution_clock::now();
//            duration = duration_cast<microseconds>( t2 - t1 ).count();
//            cout << "Time updating next input sequence = "
//                 << duration/1000.0 << "ms" << endl << endl;
        }
        else {
            // QP failed - repeat previous solution
            cout << "QP failed" << endl;
        }

    }
}

QuadraticProblem Generator::buildQP(const Constraint& collision, const State3D& state,
                                     int agent_id) {
    // Determine if a collision happened, based on the size of the constraint
    int num_neigh = collision.b.size() / 2;
    int num_vars = _ineq.A_full.cols() + num_neigh;
    int num_ineq = _ineq.A_full.rows() + 2 * num_neigh;
    int num_eq = _eq.A.rows();

    MatrixXd Ain_full = MatrixXd::Zero(num_ineq, num_vars);
    MatrixXd Ain = MatrixXd::Zero(_ineq.A.rows() + 2 * num_neigh, num_vars);
    MatrixXd Aeq = MatrixXd::Zero(num_eq, num_vars);
    MatrixXd H = MatrixXd::Zero(num_vars, num_vars);
    VectorXd bin_full = VectorXd::Zero(num_ineq);
    VectorXd bin_lower = VectorXd::Zero(_ineq.A.rows() + 2 * num_neigh);
    VectorXd bin_upper = VectorXd::Zero(_ineq.A.rows() + 2 * num_neigh);
    VectorXd beq = VectorXd::Zero(num_eq);
    VectorXd f = VectorXd::Zero(num_vars);
    VectorXd x0 = VectorXd::Zero(2 * _dim);
    x0 << state.pos, state.vel;

    // Build the complete inequality constraint
    Ain_full << _ineq.A_full, MatrixXd::Zero(_ineq.A_full.rows(), num_neigh),
                collision.A;

    Ain << _ineq.A, MatrixXd::Zero(_ineq.A.rows(), num_neigh),
            collision.A;

    bin_lower << _ineq.lower_bound, -50*VectorXd::Ones(num_neigh), -50*VectorXd::Ones(num_neigh);
    bin_upper << _ineq.upper_bound, collision.b;

    bin_full << _ineq.b_full, collision.b;
    Aeq << _eq.A, MatrixXd::Zero(num_eq, num_neigh);

    // build the constant vector for the equality constraint
    for (int i = 0; i < _deg_poly + 1; i++)
        beq.segment(_dim * i * _l, _dim) = _x0_ref[agent_id].col(i);

    // Augment matrices in case of collisions to account for the additional slack variables
    if (collision.b.size() > 0) {
        MatrixXd H_eps = _quad_coll * MatrixXd::Identity(num_neigh, num_neigh);
        VectorXd f_eps = _lin_coll * VectorXd::Ones(num_neigh);
        H << _H_o, MatrixXd::Zero(_H_o.rows(), num_neigh),
             MatrixXd::Zero(num_neigh, _H_o.cols()), H_eps;
        f << -2 * (_fpf_obs[agent_id] - x0.transpose() * _Hlin_o).transpose(), f_eps;
    }
    else {
        H = _H_f;
        f << -2 * (_fpf_free[agent_id] - x0.transpose() * _Hlin_f).transpose();
    }

    // Assemble the problem and return it
    QuadraticProblem problem = {H, Aeq, Ain_full, Ain, f, beq, bin_full, bin_lower, bin_upper};
    return problem;
}

MatrixXd Generator::getInitialReference(const State3D &state, const MatrixXd& ref) {
    VectorXd err_pos, cost, denominator;
    err_pos = state.pos - ref.col(0);
    denominator = -(state.vel.array() + 0.01 * state.vel.array().sign());
    cost = err_pos.array().pow(5) / denominator.array();

    if ((cost.array() > _max_cost).any() ||  (cost.array() < _min_cost).any()) {
        MatrixXd new_ref = MatrixXd::Zero(_dim, _d + 1);
        new_ref << state.pos, state.vel, MatrixXd::Zero(_dim, _d - 1);
        return new_ref;
    }
    else return ref;
}

MatrixXd Generator::updateHorizon(const Eigen::VectorXd &u, const State3D& states) {

    VectorXd initial_states = VectorXd::Zero(2 * _dim);
    initial_states << states.pos, states.vel;
    return vec2mat(_Phi_pred.pos * u + _A0_pred.pos * initial_states, _dim);
}

MatrixXd Generator::updateInitialReference(const Eigen::VectorXd &u) {
    MatrixXd init_ref = MatrixXd::Zero(_dim, _d + 1);
    for (int r = 0; r < _d + 1; r++)
        init_ref.col(r) =  _Rho_h[r].middleRows(_dim, _dim) * u;

    return init_ref;
}