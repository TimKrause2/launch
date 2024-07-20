#include "psopt.h"

#define N_ELLIP_POINTS 100

struct ICBMData
{
    double mass; // mass of the secondary body
    double thrust; // acceleration thrust
    double mu; // GM of the primary body
    double mag_omega; // magnitude of the omega vector
    double ring_mass;
    double ring_radius;
    double G;
    int N_ellip_points;
    MatrixXd k;
    MatrixXd ellint_1;
    MatrixXd ellint_2;
    double T1; // burn time guess
    double T2; // coast time guess
    Eigen::Vector3d omega;
    Eigen::Vector3d r_vehicle0;
    Eigen::Vector3d v_vehicle0;
    Eigen::Vector3d r_target0;
    Eigen::Vector3d v_target0;
    Eigen::Vector3d r_boost;
    Eigen::Vector3d v_boost;
    Eigen::Vector3d ecc_v; // eccentricity vector
    Eigen::Vector3d n_h; // angular momentum vector (noramlzied)
    double a; // semi-major length
    double e; // eccentricity
    double i; // inclination
    double Om; // right ascension of the ascending node
    double om; // argument of perigee

};

bool icbm_launch(
        ICBMData *data,
        Eigen::MatrixXd &controls, // phase 1 controls
        Eigen::MatrixXd &time); // phase 1 time
