#include "spheroid.h"
#include <boost/math/special_functions/ellint_1.hpp>
#include <boost/math/special_functions/ellint_2.hpp>
#include <iostream>

Spheroid::Spheroid(Spheroid &s)
    :
      Body(
          s.m_mass,
          s.m_radius,
          s.m_position,
          s.m_velocity,
          s.p,
          s.pdot,
          s.Jp,
          NULL),
      ring_mass(s.ring_mass),
      ring_radius(s.ring_radius),
      omega_p(s.omega_p)
{}

Spheroid::Spheroid(
        double mass,
        double radius,
        Eigen::Vector3d &position,
        Eigen::Vector3d &velocity,
        Eigen::Vector4d &p,
        Eigen::Vector4d &pdot,
        Eigen::Matrix3d &Jp,
        char* tex_image,
        double ring_mass,
        double ring_radius):
    Body(
        mass,
        radius,
        position,
        velocity,
        p,
        pdot,
        Jp,
        tex_image),
    ring_mass(ring_mass),
    ring_radius(ring_radius)
{
    omega_p = caams::omega_p_p_dot(p,pdot);
}

void Spheroid::Update(double dt)
{
    pdot = caams::p_dot_omega_p(p, omega_p);
}

Eigen::Vector3d Spheroid::rk_acceleration(
        Eigen::Vector3d const &rk_position,
        Eigen::Vector3d const &rk_velocity)
{
    // calculate the acceleration due to the point mass
    Eigen::Vector3d a_point = Body::rk_acceleration(rk_position, rk_velocity);

    // calculate the acceleration due to the ring
    Eigen::Vector3d r = rk_position - m_rk_position;
    Eigen::Matrix3d A_basis = caams::Ap(rk_p);
    Eigen::Vector3d rp = A_basis.transpose()*r;
    Eigen::Vector3d a_ring_p = g_accel_ring(
                ring_mass,
                ring_radius,
                rp);
    Eigen::Vector3d a_ring = A_basis*a_ring_p;

    return a_point + a_ring;
}

void Spheroid::SetMass(
        double g_ref,
        double r_equ,
        double one_over_f)
{
    double f = 1/one_over_f;
    double r_pole = r_equ*(1-f);
    double r_avg = (r_pole + 2*r_equ)/3.0;
    Eigen::Matrix2d A;
    A(0,0) = G_gravity/r_pole/r_pole;
    A(0,1) = G_gravity*r_pole/pow(r_avg*r_avg + r_pole*r_pole, 3.0/2.0);
    A(1,0) = G_gravity/r_equ/r_equ;
    Eigen::Vector3d r_equ_v(r_equ, 0.0, 0.0);
    Eigen::Vector3d a_ring_v = g_accel_ring(1.0, r_avg, r_equ_v);
    A(1,1) = a_ring_v.norm();
    Eigen::Vector2d y(g_ref,g_ref);
    Eigen::Vector2d x = A.inverse()*y;

    std::cout << "M_point:" << x(0) << std::endl;
    std::cout << "M_ring :" << x(1) << std::endl;

    m_mass = x(0);
    ring_mass = x(1);
    m_radius = r_equ;
    ring_radius = r_avg;
}


Eigen::Vector3d g_accel_ring(
        double ring_mass,
        double ring_radius,
        Eigen::Vector3d rp)
{
    Eigen::Vector3d rs = rp/ring_radius;
    double rho = sqrt(rs(0)*rs(0) + rs(1)*rs(1));
    double z = rs(2);
    double length = 2*M_PI*ring_radius;
    double lambda = ring_mass/length;
    double A = -2*G_gravity*lambda/ring_radius;
    double B = (1+rho)*(1+rho) + z*z;
    double m = 4*rho/B;
    double k = sqrt(m);
    double C = sqrt(B)*((rho-1)*(rho-1)+z*z);
    double K_m = boost::math::ellint_1(k);
    double E_m = boost::math::ellint_2(k);
    double F_rho;
    if(rho==0.0){
        F_rho = 0.0;
    }else{
        F_rho = A/rho/C*((rho*rho-1-z*z)*E_m +
                         ((1-rho)*(1-rho)+z*z)*K_m);
    }
    double F_z = A*2*z/C*E_m;
    Eigen::Vector3d n_rho;
    if(rho==0){
        n_rho = Eigen::Vector3d::Zero();
    }else{
        n_rho << rs.head<2>(), 0.0;
        n_rho.normalize();
    }
    Eigen::Vector3d n_z(0,0,1);
    Eigen::Vector3d a_ring_p = F_rho*n_rho + F_z*n_z;
    return a_ring_p;
}

