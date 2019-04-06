//
// Created by carlos on 05/04/19.
//

#include <avoidance.h>

using namespace std;
using namespace Eigen;

OndemandAvoider::OndemandAvoider (const std::vector<Eigen::MatrixXd>& horizon,
                                  const Eigen::MatrixXd& Phi_pos,
                                  const Eigen::MatrixXd& A0_pos,
                                  const std::vector<EllipseParams>& p) :
    _horizon(horizon),
    _Phi_pos(Phi_pos),
    _A0_pos(A0_pos)
{
    Ellipse tmp;
    MatrixXd E;
    for (int i = 0; i < p.size(); i++) {
        tmp.order = p[i].order;
        tmp.rmin = p[i].rmin;
        E = p[i].c.asDiagonal();
        tmp.E1 = E.inverse();
        tmp.E2 = tmp.E1.array().pow(2);
        _ellipse.push_back(tmp);
    }

    _k_hor = _horizon[0].cols();
    _dim = _horizon[0].rows();
    _N = _horizon.size();
}

vector<bool> OndemandAvoider::check_collisions(int agent_id, int k) {
    vector<bool> violation(_N);
    vector<bool> neighbours(_N);
    VectorXd pi = _horizon[agent_id].col(k);
    VectorXd pj = VectorXd::Zero(_dim);
    VectorXd diff = VectorXd::Zero(_dim);
    double dist;
    int order;

    for (int j = 0; j < _N; j++) {
        if (j != agent_id) {
            order = _ellipse[j].order;
            pj = _horizon[j].col(k);
            diff = _ellipse[j].E1 * (pi - pj);
            dist = pow(((diff.array().pow(order)).sum()), 1.0 / order);
            violation[j] = dist < _ellipse[j].rmin;
            neighbours[j] = dist < _ellipse[j].rmin * 2;
        }
        else {
            violation[j] = false;
            neighbours[j] = false;
        }
    }

    if (std::find(violation.begin(), violation.end(), true) == violation.end())
        neighbours.clear();

    return neighbours;

}

CollisionConstraint OndemandAvoider::build_coll_constr(int agent_id, int k,
                                                       const std::vector<bool> &neighbours,
                                                       const State3D &state) {
    int num_neighbours = std::count(neighbours.begin(), neighbours.end(), true);
    int num_variables = _Phi_pos.cols();
    MatrixXd Ain = MatrixXd::Zero(num_neighbours, num_variables);
    VectorXd bin = VectorXd::Zero(num_neighbours);
    MatrixXd diff_row = MatrixXd::Zero(1, _dim * _k_hor);
    VectorXd pi = _horizon[agent_id].col(k);
    VectorXd pj = VectorXd::Zero(_dim);
    VectorXd diff = VectorXd::Zero(_dim);
    double dist, r;
    int order;

    VectorXd initial_states = VectorXd::Zero(2 * _dim);
    VectorXd distance = VectorXd::Zero(num_neighbours);
    initial_states << state.pos, state.vel;
    int idx = 0;

    for (int j = 0; j < _N; j++) {
        if (neighbours[j]) {
            order = _ellipse[j].order;
            pj = _horizon[j].col(k);
            diff = _ellipse[j].E1 * (pi - pj);
            dist = pow(((diff.array().pow(order)).sum()), 1.0 / order);
            diff = (_ellipse[j].E2 * (pi - pj)).array().pow(order - 1);

            r = pow(dist, order - 1) * (_ellipse[j].rmin - dist) + diff.transpose() * pi
                    - diff.transpose() * _A0_pos.middleRows(_dim * k, _dim) * initial_states;

            distance[idx] = pow(dist, order - 1);

            diff_row << MatrixXd::Zero(1, _dim * k),
                        diff.transpose(),
                        MatrixXd::Zero(1, _dim * (_k_hor - k - 1));

            Ain.row(idx) = -diff_row * _Phi_pos;
            bin[idx] = -r;
            idx++;
        }
    }

    // Assemble constraint and return it
    Constraint constraint = {Ain, bin};
    CollisionConstraint collision = {constraint, distance};
    return collision;
}

Constraint OndemandAvoider::get_coll_constr(const State3D &state, int agent_id) {
    // iterate through each time step
    vector<bool> neighbours;
    CollisionConstraint collision;
    Constraint soft_constraint;
    for (int k = 0; k < _k_hor; k++) {
        neighbours = check_collisions(agent_id, k);
        bool violation = neighbours.size() > 0;

        if (violation) {
            collision = build_coll_constr(agent_id, k, neighbours, state);

            int num_neighbours = collision.constraint.b.size();
            int num_variables = _Phi_pos.cols();
            MatrixXd Ain_soft = MatrixXd::Zero(2 * num_neighbours, num_variables + num_neighbours);
            VectorXd bin_soft = VectorXd::Zero(2 * num_neighbours);

            MatrixXd diag_distance = collision.distance.asDiagonal();

            Ain_soft << collision.constraint.A, diag_distance,
                        MatrixXd::Zero(num_neighbours, num_variables),
                        MatrixXd::Identity(num_neighbours, num_neighbours);

            bin_soft << collision.constraint.b, VectorXd::Zero(num_neighbours);
            soft_constraint = {Ain_soft, bin_soft};
            break;
        }
    }
    return soft_constraint;
}
