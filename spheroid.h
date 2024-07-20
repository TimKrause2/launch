#include "gsim.h"

class Spheroid : public Body
{
public:
    double ring_mass;
    double ring_radius;
    Eigen::Vector3d omega_p;
    Spheroid(Spheroid &s);
    Spheroid(
            double mass,
            double radius,
            Eigen::Vector3d &position,
            Eigen::Vector3d &velocity,
            Eigen::Vector4d &p,
            Eigen::Vector4d &pdot,
            Eigen::Matrix3d &Jp,
            char* tex_image,
            double ring_mass,
            double ring_radius);

    void Update(double dt);

    Eigen::Vector3d rk_acceleration(
            Eigen::Vector3d const &rk_position,
            Eigen::Vector3d const &rk_velocity);

    void SetMass(
            double g_ref,
            double r_equ,
            double one_over_f);

};

Eigen::Vector3d g_accel_ring(
        double ring_mass,
        double ring_radius,
        Eigen::Vector3d rp);
