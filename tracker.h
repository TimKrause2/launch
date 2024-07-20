#include "gsim.h"

class Tracker : public Body
{
    Eigen::Vector3d vTrack;
    double theta_x;
    double theta_y;
    double torque_b_z;
public:
    Tracker(double mass,
            double radius,
            Eigen::Vector3d &position,
            Eigen::Vector3d &velocity,
            Eigen::Vector4d &p,
            Eigen::Vector4d &pdot,
            Eigen::Matrix3d &Jp,
            char* tex_image);
    void Prepare(void);
    void ForceAndTorque(
            double dt,
            Eigen::Vector3d &force,   // force in global coordinate space
            Eigen::Vector3d &torque); // torque in body coordinate space
    void setVTrack(Eigen::Vector3d const &vTrack);
};
