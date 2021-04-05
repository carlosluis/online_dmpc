//
// Created by carlos on 09/04/19.
//
#include "iostream"
#include "simulator.h"

using namespace Eigen;
using namespace std;
using namespace std::chrono;
using json = nlohmann::json;

Simulator::Simulator(std::ifstream& config_file) {
    Generator::Params p = parseJSON(config_file);
    _generator = std::make_unique<Generator>(p);
    _sim_model = std::make_unique<DoubleIntegrator3D>(p.mpc_params.Ts, p.model_params);

    _h = p.mpc_params.h;
    _Ts = p.mpc_params.Ts;
    _po = p.po;
    _pf = p.pf;
    _pmin = p.mpc_params.limits.pmin;
    _pmax = p.mpc_params.limits.pmax;
    _current_states.reserve(_Ncmd);
    _trajectories.reserve(_Ncmd);

    _ellipses = _generator->getEllipses();

    for (int i = 0; i < _Ncmd; i++) {
        State3D agent_i = {_po.col(i), VectorXd::Zero(_po.rows())};
        _current_states.push_back(addRandomNoise(agent_i));
    }
}

Generator::Params Simulator::parseJSON(std::ifstream& config_file) {
    json j;
    config_file >> j;

    // Parse the data into variables and create the generator params object
    _N = j["N"];
    _Ncmd = j["Ncmd"];

    std::string solver_name =  j["solver"];
    Solver qp_solver;

    if (!solver_name.compare("qpoases"))
        qp_solver = kQpoases;
    else
        throw std::invalid_argument("Invalid solver '" + solver_name + " '");

    // Bezier curve parameters
    _bezier_params = {j["d"], j["num_segments"], j["dim"], j["deg_poly"], j["t_segment"]};

    // Dynamic model used in MPC
    _model_params = {j["zeta_xy"], j["tau_xy"], j["zeta_z"], j["tau_z"]};

    // Tuning parameters for MPC
    VectorXd energy_weights = VectorXd::Zero(static_cast<int>(j["d"]) + 1);
    energy_weights(2) = j["acc_cost"];

    TuningParams tune = {j["s_free"], j["s_obs"], j["s_repel"], j["spd_f"], j["spd_o"],
                         j["spd_r"], j["lin_coll"], j["quad_coll"], energy_weights};

    // Physical limits for inequality constraint
    vector<double> pmin_vec = j["pmin"];
    vector<double> pmax_vec = j["pmax"];
    vector<double> amin_vec = j["amin"];
    vector<double> amax_vec = j["amax"];

    Vector3d pmin(pmin_vec.data());
    Vector3d pmax(pmax_vec.data());
    Vector3d amin(amin_vec.data());
    Vector3d amax(amax_vec.data());

    PhysicalLimits limits = {pmax, pmin, amax, amin};

    _mpc_params = {j["h"], j["ts"], j["k_hor"], tune, limits};

    // Ellipsoidal collision constraint
    EllipseParams ellipse_params;
    ellipse_params.order = j["order"];
    ellipse_params.rmin = j["rmin"];
    ellipse_params.c = (Eigen::Vector3d() << 1.0, 1.0, j["height_scaling"]).finished();
    std::vector<EllipseParams> ellipse_vec(_Ncmd, ellipse_params);

    EllipseParams ellipse_params_obs;
    ellipse_params_obs.order = j["order_obs"];
    ellipse_params_obs.rmin = j["rmin_obs"];
    ellipse_params_obs.c = (Eigen::Vector3d() << 1.0, 1.0, j["height_scaling_obs"]).finished();
    std::vector<EllipseParams> ellipse_vec_obs(_N - _Ncmd, ellipse_params_obs);

    // Combine the two
    ellipse_vec.insert(ellipse_vec.end(), ellipse_vec_obs.begin(), ellipse_vec_obs.end());

    _pos_std = j["std_position"];
    _vel_std = j["std_velocity"];

    // Generate the test
    std::string test_type =  j["test"];
    if (!test_type.compare("default")) {
        // get the values from the json file for a fixed test
        vector<vector<double>> po_vec = j["po"];
        _po = MatrixXd::Zero(j["dim"], _N);
        for (int i = 0; i < _N; i++) {
            Vector3d poi(po_vec[i].data());
            _po.col(i) = poi;
        }

        vector<vector<double>> pf_vec = j["pf"];
        _pf = MatrixXd::Zero(j["dim"], _Ncmd);
        for (int i = 0; i < _Ncmd; i++) {
            Vector3d pfi(pf_vec[i].data());
            _pf.col(i) = pfi;
        }

    }
    else if (!test_type.compare("random")) {
        // Generate a random test
        _po = generateRandomPoints(_N, limits.pmin.array() + 0.3,
                                   limits.pmax.array() - 0.3, ellipse_params.rmin + 0.2);

        _pf = generateRandomPoints(_Ncmd, limits.pmin.array() + 0.3,
                                   limits.pmax.array() - 0.3, ellipse_params.rmin + 0.2);
    }
    else throw std::invalid_argument("Invalid test type '" + test_type + " '");

    Generator::Params p = {_bezier_params, _model_params, ellipse_vec,
                           _mpc_params, _po, _pf, qp_solver};
    return p;
}

void Simulator::run(int duration) {

    // Amount of time steps to simulate, with time step duration  = _Ts
    auto K = static_cast<int>(floor(duration / _Ts));

    // Allocate the memory to record the trajectory of each agent
    for (int i = 0; i < _Ncmd; i++) {
        _trajectories.push_back(MatrixXd::Zero(_po.rows(), K));
    }

    // Calculate the amount of time steps to wait before running again the optimizer
    auto max_count = static_cast<int>(_h / _Ts);
    int count = max_count;
    high_resolution_clock::time_point t1, t2;

    // Loop through all the time steps for the duration specified
    for (int k = 0; k < K; k++) {
        if (count == max_count) {
            // Get next series of inputs
            t1 = high_resolution_clock::now();
            _inputs = _generator->getNextInputs(_current_states);
            t2 = high_resolution_clock::now();
            auto mpc_duration = duration_cast<microseconds>( t2 - t1 ).count();
            cout << "Solving frequency = " << 1000000.0 / mpc_duration << " Hz" << endl;
            count = 0;
        }

        // Apply inputs to all agents to get the next states
        for (int i = 0; i < _Ncmd; i++) {
            State3D sim_states = _sim_model->applyInput(_current_states[i], _inputs[i].col(count));
            _current_states[i] = addRandomNoise(sim_states);
            _trajectories[i].col(k) = _current_states[i].pos;
        }
        count++;
    }

    // After simulation finished, run a check on the solution
    collisionCheck(_trajectories);
    goalCheck(_current_states);
}


bool Simulator::collisionCheck(const std::vector<Eigen::MatrixXd> &trajectories) {
    float rmin_check = 0.15;
    int order = 2;
    VectorXd c_check = (Eigen::Vector3d() << 1.0, 1.0, 3.0).finished();
    MatrixXd E_check = c_check.asDiagonal();
    MatrixXd E1_check = E_check.inverse();

    MatrixXd differ;
    VectorXd dist;
    bool violation = false;
    double min_dist;
    int pos;

    for (int i = 0; i < _Ncmd; i++) {
        for (int j = i + 1; j < _Ncmd; j++) {
            if (i != j) {
                differ = E1_check * (trajectories[i] - trajectories[j]);
                dist = pow(((differ.array().pow(order)).colwise().sum()),1.0 / order);
                min_dist = dist.minCoeff(&pos);

                if (min_dist < rmin_check) {
                    violation = true;
                    cout << "Collision constraint violation: ";
                    cout << "Vehicles " << i << " and " << j;
                    cout << " will be " << min_dist << "m";
                    cout << " apart @ t = " << (float)pos * _Ts << "s" << endl;
                }
            }
        }
    }

    if (!violation)
        cout << "No collisions found!" << endl;
    return violation;
}

bool Simulator::goalCheck(const std::vector<State3D> &states) {
    Vector3d diff;
    double dist;
    bool reached_goal = true;
    float goal_tolerance = 0.1;  // 10 centimeters tolerance

    for (int i = 0; i < _Ncmd; i++) {
        diff = states[i].pos - _pf.col(i);
        dist = pow(((diff.array().pow(2)).sum()),1.0 / 2);
        if (dist > goal_tolerance){
            cout << "Vehicle " << i << " did not reached its goal by " << dist << " m" << endl;
            reached_goal = false;
        }
    }

    if (reached_goal)
        cout << "All the vehicles reached their goals!" << endl;

    return reached_goal;
}

State3D Simulator::addRandomNoise(const State3D &states) {
    std::random_device rd;
    std::mt19937 gen(rd());
    VectorXd state_vector = VectorXd::Zero(6);
    state_vector << states.pos, states.vel;
    double sample = 0.0;
    std::normal_distribution<double> distribution_position(0.0, _pos_std);
    std::normal_distribution<double> distribution_velocity(0.0, _vel_std);
    for (int i = 0; i < 6; i++) {
        if (i < 3)
            sample = distribution_position(gen);
        else
            sample = distribution_velocity(gen);

        state_vector[i] += sample;
    }

    State3D result = {state_vector.segment(0, 3), state_vector.segment(3, 3)};
    return result;
}

void Simulator::saveDataToFile(char const *pathAndName) {

    ofstream file(pathAndName, ios::out | ios::trunc);
    if(file)  // succeeded at opening the file
    {
        // instructions
        cout << "Writing solution to text file..." << endl;

        // write a few simulation parameters used in the reading end
        file << _N << " " <<  _Ncmd << " " << _pmin.transpose() << " " << _pmax.transpose() << endl;
        file << _po << endl;
        file << _pf << endl;

        for(int i=0; i < _Ncmd; ++i)
            file << _trajectories[i] << endl;

        file.close();  // close the file after finished
    }

    else
    {
        cerr << "Error while trying to open file" << endl;
    }

}

MatrixXd Simulator::generateRandomPoints(int N, const Vector3d &pmin,
                                         const Vector3d &pmax, float rmin) {
    MatrixXd pts = MatrixXd::Zero(3, N);
    Vector3d candidate = MatrixXd::Zero(3, 1);
    VectorXd dist;
    bool pass = false;

    // Generate first point
    pts.col(0) = pmin.array()
                 + (pmax - pmin).array() *
                   ((MatrixXd::Random(3, 1).array() + 1) / 2);

    for (int n = 1; n < N; ++n) {
        while (!pass) {
            // Candidate picked randomly within workspace boundaries
            candidate = pmin.array()
                        + (pmax - pmin).array() *
                          ((MatrixXd::Random(3, 1).array() + 1) / 2);

            // Calculate distance to every previous pts calculated
            dist = ((((pts.leftCols(n)).colwise()
                      -
                      candidate).array().square()).colwise().sum()).array().sqrt();

            // If the candidate is sufficiently separated from previous pts,
            // then we add it to the Matrix of valid pts
            for (int k = 0; k < n; ++k) {
                pass = dist[k] > rmin;
                if (!pass)
                    break;
            }
            if (pass)
                pts.col(n) = candidate.array();
        }
        pass = false;
    }
    return pts;
}