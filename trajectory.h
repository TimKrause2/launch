#include "caams.hpp"
#include "gsim.h"

class Trajectory{
    int Npoints;
    Eigen::MatrixXd positions;
    Eigen::MatrixXd velocities;
public:
    Trajectory(int Npoints);
    void SetPoint(int i,Body *body);
    void GetPoint(int i,Body *body);
    void Draw(int Npoints);

};
