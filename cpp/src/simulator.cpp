//
// Created by carlos on 09/04/19.
//
#include "iostream"
#include "simulator.h"

using namespace Eigen;
using namespace std;
using namespace std::chrono;

Simulator::Simulator(const Simulator::Params& p) :
    _generator(p.generator_params),
    _sim_model(p.generator_params.mpc_params.Ts, p.generator_params.model_params),
    _pos_std(p.position_noise_std),
    _vel_std(p.velocity_noise_std)
{
    // Initialize the current states with the initial locations
    _h = p.generator_params.mpc_params.h;
    _Ts = p.generator_params.mpc_params.Ts;
    _po = p.generator_params.po;
    _pf = p.generator_params.pf;
    _pmin = p.generator_params.mpc_params.limits.pmin;
    _pmax = p.generator_params.mpc_params.limits.pmax;
    _N = _po.cols();
    _Ncmd = _pf.cols();
    _current_states.reserve(_Ncmd);
    _trajectories.reserve(_Ncmd);

    _ellipses = _generator.getEllipses();

    for (int i = 0; i < _Ncmd; i++) {
        State3D agent_i = {_po.col(i), VectorXd::Zero(_po.rows())};
        _current_states.push_back(addRandomNoise(agent_i));
    }

    cout << "I successfully build a new Simulator" << endl;
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
            _inputs = _generator.getNextInputs(_current_states);
            t2 = high_resolution_clock::now();
            auto mpc_duration = duration_cast<microseconds>( t2 - t1 ).count();
            cout << "Solving frequency = " << 1000000.0 / mpc_duration << " Hz" << endl;
            count = 0;
        }

        // Apply inputs to all agents to get the next states
        for (int i = 0; i < _Ncmd; i++) {
            State3D sim_states = _sim_model.applyInput(_current_states[i], _inputs[i].col(count));
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