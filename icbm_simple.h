#include "psopt.h"
#include "orbit.h"

#define N_ELLIP_POINTS 100

struct ICBMSimpleData
{
    double mass; // mass of the secondary body
    double thrust; // acceleration thrust
    double mu; // GM of the primary body
    double ring_mass;
    double ring_radius;
    double G;
    int N_ellip_points;
    MatrixXd k;
    MatrixXd ellint_1;
    MatrixXd ellint_2;
    double T1; // burn time guess
    Eigen::Vector3d r_vehicle0;
    Eigen::Vector3d v_vehicle0;
    Eigen::Vector3d r_boost;
    Eigen::Vector3d v_boost;
    Eigen::Matrix3d A_g2a; // transfrom global to analytic
    OrbitalParams params;
};

bool icbm_simple_launch(
        ICBMSimpleData *data,
        Eigen::MatrixXd &controls, // phase 1 controls
        Eigen::MatrixXd &time); // phase 1 time
