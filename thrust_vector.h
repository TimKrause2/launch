#include "psopt.h"
#include "orbit.h"

struct ThrustVectorData
{
    double mass;
    double thrust;
    double mu;
    Eigen::Vector3d sp;
    Eigen::Matrix3d Jp;
};

void thrust_vector_simple(
        Eigen::Vector3d const &r0,
        Eigen::Vector3d const &v0,
        OrbitalParams const &oe_final,
        ThrustVectorData *data,
        Eigen::MatrixXd control,
        Eigen::MatrixXd time);
