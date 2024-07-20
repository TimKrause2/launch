#include "psopt.h"
#include "orbit.h"

#define N_ELLIP_POINTS 100

struct DelayedLaunchData
{
    double mass;
    double thrust; // thrust acceleration
    double mu; // GM of the primary body
    double mag_omega; // magnitude of the angular velocity vector
    Eigen::Vector3d omega; // primary body angular velocity vector in celestial frame
    double ring_mass; // spheroid approximation ring mass
    double ring_radius; // ring radius
    double G; // gravitational constant
    double T1; // phase 1 time guess and return value
    double T2; // phase 2 time; fixed
    double T3; // phase 3 time guess and return value
    Eigen::Vector3d r1i; // phase 1 initial position
    Eigen::Vector3d v1i; // phase 1 initial velocity
    OrbitalParams params_final; // desired orbital parameters and guess for nu
    int N_ellip_points;
    MatrixXd k;
    MatrixXd ellint_1;
    MatrixXd ellint_2;
};

void delayed_launch(
        DelayedLaunchData *data,
        Eigen::MatrixXd  &controls, // phase 3 controls; thrust vector
        Eigen::MatrixXd  &time); // phase 3 time
