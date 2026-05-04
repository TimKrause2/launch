#include <iostream>
#include "gsim.h"

double mass_earth = 5.97219e24;
double radius_earth = 6.371e6;
double mag_omega_earth = (2*M_PI/86164.0905);

double body_rot_energy(Body* body)
{
    Eigen::Vector3d omega_p = caams::omega_p_p_dot(
                body->p, body->pdot);
    return 0.5*omega_p.dot(body->Jp*omega_p);
}





int main(void)
{
    Eigen::Vector3d pos_earth = Eigen::Vector3d::Zero();
    Eigen::Vector3d vel_earth;
    vel_earth << 1.0, 0.0, 0.0;
    Eigen::Vector4d p_identity;
    p_identity << 1.0, 0.0, 0.0, 0.0;
    Eigen::Vector4d p_earth = p_identity;
    Eigen::Vector3d omega_p_earth(0.0, 0.0, mag_omega_earth);
    Eigen::Vector4d pdot_earth = caams::p_dot_omega_p(p_earth, omega_p_earth);
    Eigen::Matrix3d Jp_earth = caams::J_p_sphere(mass_earth,radius_earth);

    Body *body = new Body(
                mass_earth,
                radius_earth,
                pos_earth,
                vel_earth,
                p_earth,
                pdot_earth,
                Jp_earth,
                NULL);

    System system;

    system.AddBody(body);


    double Krot0 = body_rot_energy(body);

    std::cout.precision(17);

    for(;;){
        system.rkIntegrate(1.0/60.0);
        double Krot = body_rot_energy(body);
        std::cout << "ratio:" << Krot/Krot0 << std::endl;
        std::cout << "position:\n" << body->m_position << std::endl;
        std::cout << "velocity:\n" << body->m_velocity << std::endl;
    }

    return 0;
}
