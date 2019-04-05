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
}

Constraint OndemandAvoider::get_coll_constr(const State3D &state, const int &agent_id) {

}
