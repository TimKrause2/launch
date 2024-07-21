#include "orbit.h"
#include <math.h>
#include <iostream>
#include <limits>
#include <cfloat>
#include <GL/gl.h>

#define G_gravity 6.67384e-11
#define E_MIN  1.0e-20

void ModEquOrbitalParams::FromOrbitalParams(OrbitalParams &c)
{
    p = c.a*(1.0-c.e*c.e);
    f = c.e*cos(c.omega + c.Omega);
    g = c.e*sin(c.omega + c.Omega);
    h = tan(c.i/2.0)*cos(c.Omega);
    k = tan(c.i/2.0)*sin(c.Omega);
    L = c.Omega + c.omega + c.nu;
}

void ModEquOrbitalParams::ToOrbitalParams(OrbitalParams &c)
{
    double f2pg2 = f*f + g*g;
    c.a = p/(1.0 - f2pg2);
    c.e = sqrt(f2pg2);
    double h2pk2 = h*h + k*k;
    c.i = 2.0*atan(sqrt(h2pk2));//atan2(2.0*sqrt(h2pk2), 1 - h2pk2);
    c.Omega = atan2(k, h);
    c.omega = atan2(g*h - f*k, f*h + g*k);
    c.nu = L - (c.Omega + c.omega);
    if(c.nu>M_PI){
        c.nu -= 2.0*M_PI;
    }
    if(c.nu<=-M_PI){
        c.nu += 2.0*M_PI;
    }
    if(c.omega>M_PI){
        c.omega -= 2.0*M_PI;
    }
    if(c.omega<=-M_PI){
        c.omega += 2.0*M_PI;
    }
}

Eigen::VectorXd ModEquOrbitalParams::ToVector(void)
{
    Eigen::VectorXd v(6,1);
    v << p, f, g, h, k, L;
    return v;
}

void ModEquOrbitalParams::FromVector(Eigen::VectorXd v)
{
    p = v(0);
    f = v(1);
    g = v(2);
    h = v(3);
    k = v(4);
    L = v(5);
}

void orbital_plane_elements(
        Eigen::Vector3d const &h,
        Eigen::Vector3d const &e,
        OrbitalParams &params)
{
    Eigen::Vector3d n_h = h.normalized();
    Eigen::Vector3d n_e = e.normalized();
    Eigen::Vector3d n_an;
    params.i = acos(n_h(2));
    if(params.i==0.0){
        n_an = Eigen::Vector3d(1,0,0);
        params.Omega = 0.0;
    }else{
        Eigen::Vector3d z_axis(0,0,1);
        n_an = z_axis.cross(n_h).normalized();
        params.Omega = acos(n_an(0));
        if(n_an(1)<0.0)
            params.Omega*=-1.0;
    }
    params.omega = acos(n_an.dot(n_e));
    if(n_e(2)<0.0)
        params.omega*=-1.0;
}

void orbital_elements(
        Eigen::Vector3d r,
        Eigen::Vector3d v,
        double mass_pri,
        double mass_sec,
        OrbitalParams &params)
{
    double mag_r = r.norm();
    double mag_v = v.norm();

    // calculate the specific anglular momentum
    Eigen::Vector3d h = r.cross(v);
    Eigen::Vector3d h_norm = h.normalized();

    // calculate the standard gravitational parameter
    double mu = G_gravity*(mass_pri+mass_sec);

    // calculate the semi-major axis length
    double E = mag_v*mag_v/2.0 - mu/mag_r;
    params.a = -mu/2.0/E;

    // compute the eccentricity vector
    //caams::matrix e_vec(caams::SS(v)*h*(1.0/mu)-r*(1.0/mag_r));
    Eigen::Vector3d e_vec = v.cross(h)/mu - r.normalized();
    params.e = e_vec.norm();

    orbital_plane_elements(h, e_vec, params);

    params.nu = acos(e_vec.dot(r)/mag_r/params.e);
    if(r.dot(v)<0.0)
        params.nu*=-1.0;


    // // compute the node vector
    // //caams::matrix node_vec(caams::SS(caams::matrix(3,1,0.0,0.0,1.0))*h);
    // Eigen::Vector3d z_axis;
    // z_axis << 0.0, 0.0, 1.0;
    // Eigen::Vector3d node_vec = z_axis.cross(h);

    // double mag_node_vec = node_vec.norm();
    // if(mag_node_vec == 0.0){
    //     // orbit is not inclined
    //     // choose Omega arbitrarily
    //     params.Omega = 0.0;
    //     if(h_norm(2)>0.0){
    //         params.i=0.0;
    //     }else{
    //         params.i=M_PI;
    //     }
    //     params.omega = atan2(e_vec(1),e_vec(0));
    //     params.nu = atan2(r(1),r(0)) - params.omega;
    // }else{
    //     // inclined
    //     params.Omega = atan2(node_vec(1),node_vec(0));

    //     params.i = acos(h_norm(2));

    //     // calculate the y-axis for the orbit
    //     //caams::matrix node_x(node_vec*(1.0/mag_node_vec));
    //     //caams::matrix node_y(caams::SS(h_norm)*node_x);
    //     Eigen::Vector3d node_x = node_vec.normalized();
    //     Eigen::Vector3d node_y = h_norm.cross(node_x);

    //     //caams::matrix r_dot_node_x(~r*node_x);
    //     //caams::matrix r_dot_node_y(~r*node_y);
    //     double r_dot_node_x = r.dot(node_x);
    //     double r_dot_node_y = r.dot(node_y);
    //     params.nu = atan2(r_dot_node_y,r_dot_node_x);

    //     double e_dot_node_x = e_vec.dot(node_x);
    //     double e_dot_node_y = e_vec.dot(node_y);
    //     params.omega = atan2(e_dot_node_y,e_dot_node_x);
    //     params.nu -= params.omega;
    // }

    // if(params.nu>M_PI)params.nu-=2.0*M_PI;
    // if(params.nu<-M_PI)params.nu+=2.0*M_PI;
}

Eigen::Vector3d eccentricity_vector(
    Eigen::Vector3d r, // location of satellite to primary body
    Eigen::Vector3d v, // velocity of satellite relative to primary
    double mass_pri,
    double mass_sec)
{
    Eigen::Vector3d h = r.cross(v);

    // calculate the standard gravitational parameter
    double mu = G_gravity*(mass_pri+mass_sec);

    // compute the eccentricity vector
    Eigen::Vector3d e_vec = v.cross(h)/mu - r.normalized();
    return e_vec;
}

void orbit_state(
    OrbitalParams &params,
    double mass_pri,
    double mass_sec,
    Eigen::Vector3d &r,
    Eigen::Vector3d &v)
{
    double mu=G_gravity*(mass_pri+mass_sec);

    // positive x-axis is in the same direction as the eccentricity vector
    // calculate the position and velocity in the
    double cos_nu = cos(params.nu);
    double sin_nu = sin(params.nu);
    double mag_r = params.a*(1-params.e*params.e)/(1+params.e*cos_nu);
    Eigen::Vector3d r_orbit;
    r_orbit << mag_r*cos_nu, mag_r*sin_nu, 0.0;
    // the tangent to the curve
    double A = params.a*(1-params.e*params.e);
    Eigen::Vector3d v_norm(-sin_nu, params.e + cos_nu, 0.0);
    double mag_v = sqrt(mu/A);
    Eigen::Vector3d v_orbit = v_norm*mag_v;

    // transform into the universal coordinate system
    Eigen::Matrix3d A_basis = orbit_basis(params);

    r = A_basis*r_orbit;
    v = A_basis*v_orbit;
}

void orbit_initialize(
        OrbitalParams &params,
        double mass_pri, // mass of the primary in kilograms
        double mass_sec, // mass of the secondary in kilograms
        Eigen::Vector3d &pos_pri, // position of the primary in meters
        Eigen::Vector3d &pos_sec, // position of the secondary in meters
        Eigen::Vector3d &vel_pri, // velocity of the primary in meters/second
        Eigen::Vector3d &vel_sec) // velocity of the secondary in meters/second
{
    Eigen::Vector3d r;
    Eigen::Vector3d v;
    orbit_state( params, mass_pri, mass_sec, r, v);

//    double mu=G_gravity*(mass_pri+mass_sec);

//    // positive x-axis is in the same direction as the eccentricity vector
//    // calculate the position and velocity in the
//    double cos_nu = cos(params.nu);
//    double sin_nu = sin(params.nu);
//    double mag_r = params.a*(1-params.e*params.e)/(1+params.e*cos_nu);
//    //caams::matrix r(3,1,mag_r*cos_nu,mag_r*sin_nu,0.0);
//    Eigen::Vector3d r;
//    r << mag_r*cos_nu, mag_r*sin_nu, 0.0;
//    // the tangent to the curve
//    double den = 1 + params.e*cos_nu;
//    double A = params.a*(1-params.e*params.e);
//    //caams::matrix dr(3,1,
//    //                 A*params.e*cos_nu*sin_nu/den/den-A*sin_nu/den,
//    //                 A*params.e*sin_nu*sin_nu/den/den+A*cos_nu/den,
//    //                 0.0);
//    Eigen::Vector3d dr;
//    dr << A*params.e*cos_nu*sin_nu/den/den-A*sin_nu/den,
//          A*params.e*sin_nu*sin_nu/den/den+A*cos_nu/den,
//          0.0;
//    Eigen::Vector3d v_norm = dr.normalized();
//    double mag_v = sqrt(mu*(2.0/mag_r - 1.0/params.a));
//    Eigen::Vector3d v = v_norm*mag_v;

    Eigen::Vector3d r_bc = r*(mass_sec/(mass_pri+mass_sec));
    Eigen::Vector3d v_bc = v*(mass_sec/(mass_pri+mass_sec));

    pos_pri = -r_bc;
    vel_pri = -v_bc;
    pos_sec = r-r_bc;
    vel_sec = v-v_bc;

}

void polar_orbit_initialize(
        OrbitalParams &params,
        double mass_pri,
        double mass_sec,
        double &theta,
        double &theta_dot,
        double &r,
        double &r_dot)
{
    double cos_nu = cos(params.nu);
    double sin_nu = sin(params.nu);
    theta = params.nu;
    r = params.a*(1.0-params.e*params.e)/(1.0+params.e*cos_nu);
    double mu=G_gravity*(mass_pri+mass_sec);
    // the tangent to the curve
    double A = params.a*(1-params.e*params.e);
    Eigen::Vector2d v_norm(-sin_nu, params.e + cos_nu);
    double mag_v = sqrt(mu/A);
    Eigen::Vector2d v = v_norm*mag_v;
    Eigen::Vector2d n_r;
    n_r << cos_nu, sin_nu;
    Eigen::Vector2d n_theta;
    n_theta << -sin_nu, cos_nu;
    r_dot = v.dot(n_r);
    theta_dot = v.dot(n_theta)/r;
}

void angular_limits(
        double theta1,
        double theta2,
        double &theta_init,
        double &theta_final)
{
    if(theta1>=0){
        if(theta2>=0){
            if((theta2-theta1)>0){
                theta_init = theta1;
                theta_final = theta2;
            }else{
                theta_init = theta1;
                theta_final = theta2+2*M_PI;
            }
        }else{
            theta_init = theta1;
            theta_final = theta2+2*M_PI;
        }
    }else{
        if(theta2>=0){
            theta_init = theta1;
            theta_final = theta2;
        }else{
            if((theta2-theta1)>0){
                theta_init = theta1;
                theta_final = theta2;
            }else{
                theta_init = theta1;
                theta_final = theta2+2*M_PI;
            }
        }
    }
}


/*
 * matrix format for rk4_df
 * it is 4 rows by 1 column
 * 0 theta
 * 1 theta_dot
 * 2 r
 * 3 r_dot
 */


Eigen::Vector4d rk4_polar_df(Eigen::Vector4d x, double mu){
    Eigen::Vector4d df;
    double theta = x(0);
    double theta_dot = x(1);
    double r = x(2);
    double r_dot = x(3);
//    df << x(1),
//        -2.0*x(1)*x(3)/x(2),
//        x(3),
//        x(1)*x(1)*x(2) - mu/x(2)/x(2);
    df << theta_dot,
        -2.0*theta_dot*r_dot/r,
        r_dot,
        theta_dot*theta_dot*r - mu/r/r;
    return df;
}

void rk4_polar_integrate(Eigen::Vector4d &x, double mu, double dt)
{
    Eigen::Vector4d k1 = rk4_polar_df(x,mu);
    Eigen::Vector4d k2 = rk4_polar_df(x+dt/2*k1,mu);
    Eigen::Vector4d k3 = rk4_polar_df(x+dt/2*k2,mu);
    Eigen::Vector4d k4 = rk4_polar_df(x+dt*k3,mu);
    x += dt/6*(k1+k2*2+k3*2+k4);
}

#define THETA_STEP_MAX (2*M_PI/4096)

void rk4_polar_integrate_adaptive(
        Eigen::Vector4d &x, double mu,
        double dt_total)
{
    double abs_dt_total = abs(dt_total);
    double sign_dt = (dt_total<0.0)?-1.0:1.0;
    double t = 0.0;
    while(t<abs_dt_total){
        double dt_max = abs(THETA_STEP_MAX/x(1));
        double dt_remain = abs_dt_total - t;
        double dt;
        if(dt_remain<dt_max){
            dt = dt_remain;
        }else{
            dt = dt_max;
        }
        rk4_polar_integrate(x,mu,dt*sign_dt);
        t+=dt;
    }
}

double orbit_time(
        OrbitalParams &params,
        double mass_pri,
        double mass_sec,
        double nu1,
        double nu2)
{
    double mu = G_gravity*(mass_pri+mass_sec);
    double nu_init;
    double nu_final;
    angular_limits(nu1,nu2,nu_init,nu_final);
    double theta;
    double theta_dot;
    double r;
    double r_dot;
    OrbitalParams params_init = params;
    params_init.nu = nu_init;
    polar_orbit_initialize(params_init,
                           mass_pri,mass_sec,
                           theta,theta_dot,
                           r,r_dot);
    //caams::matrix x(4,1,theta,theta_dot,r,r_dot);
    Eigen::Vector4d x;
    x << theta,theta_dot,r,r_dot;

    double T_period = period(params.a,mass_pri,mass_sec);
    double dt = T_period/64;
    double t=0.0;


    // get close to the final angle with coarse steps
    while(x(0)<nu_final){
        rk4_polar_integrate_adaptive(x,mu,dt);
        t+=dt;
    }

    // refine the guess
    do{
        double dnu = nu_final - x(0);
        dt = dnu/x(1);
        if(abs(dt)<(t*2.2e-16)){
            break;
        }
        rk4_polar_integrate_adaptive(x,mu,dt);
        t+=dt;
    }while(true);

    return t;
}

void orbit_integrate(
    OrbitalParams &params0,
    OrbitalParams &params1,
    double mass_pri,
    double mass_sec,
    double delta_t)
{
    double mu = G_gravity*(mass_pri+mass_sec);
    double theta;
    double theta_dot;
    double r;
    double r_dot;
    polar_orbit_initialize(params0,
                           mass_pri,mass_sec,
                           theta,theta_dot,
                           r,r_dot);
    Eigen::Vector4d x;
    x << theta, theta_dot, r, r_dot;

    rk4_polar_integrate_adaptive(x,mu,delta_t);

    params1 = params0;
    params1.nu = x(0);
}

Eigen::Matrix3d orbit_basis(OrbitalParams &params)
{
    // Eigen::Vector3d z_axis(0,0,1);
    // Eigen::Vector3d x_axis(1,0,0);
    // Eigen::Matrix3d A_periapsis = caams::AAA(params.omega,z_axis);
    // Eigen::Matrix3d A_inclination = caams::AAA(params.i,x_axis);
    // Eigen::Matrix3d A_ascending = caams::AAA(params.Omega,z_axis);
    // Eigen::Matrix3d A_basis = A_ascending*A_inclination*A_periapsis;
    // return A_basis;

    double cO = cos(params.Omega);
    double sO = sin(params.Omega);
    double co = cos(params.omega);
    double so = sin(params.omega);
    double ci = cos(params.i);
    double si = sin(params.i);

    Eigen::Matrix3d A;
    A(0,0) = cO*co - sO*ci*so;
    A(1,0) = cO*ci*so + sO*co;
    A(2,0) = si*so;

    A(0,1) = -cO*so - sO*ci*co;
    A(1,1) = cO*ci*co - sO*so;
    A(2,1) = si*co;

    A(0,2) = sO*si;
    A(1,2) = -cO*si;
    A(2,2) = ci;

    for(int row=0;row<3;row++){
        for(int col=0;col<3;col++){
            if(fabs(A(row,col))<1e-15)
                A(row,col)=0.0;
        }
    }

    return A;
}

double orbit_radius(
    OrbitalParams &params,
    double theta)
{
    double p = params.a*(1 - params.e*params.e);
    double r = p/(1 + params.e*cos(theta - params.omega));
    return r;
}

double d_orbit_radius(
    OrbitalParams &params,
    double theta)
{
    double p = params.a*(1 - params.e*params.e);
    double num = params.e*p*sin(theta - params.omega);
    double den = 1 + params.e*cos(theta - params.omega);
    double dr = num/den/den;
    return dr;
}

double refine_guess(
    OrbitalParams &params1,
    OrbitalParams &params2,
    double theta_g)
{
    while(1){
        double f = orbit_radius(params1, theta_g)
                   - orbit_radius(params2, theta_g);
        double df = d_orbit_radius(params1, theta_g)
                    - d_orbit_radius(params2, theta_g);

        double dtheta = f/df;
        if( fabs(dtheta) < fabs(theta_g)*DBL_EPSILON*2 ){
            return theta_g;
        }
        theta_g -= dtheta;
    }
}




#define N_DR_POINTS 1024

bool orbit_radii_intersect(
    OrbitalParams &params1,
    OrbitalParams &params2,
    double &theta_intersect1, // delta-r positive
    double &theta_intersect2) // delta-r negative
{
    if(params1.i!=params2.i || params1.Omega!=params2.Omega){
        std::cout << "orbit_radii_intersect: Orbits are not co-planar." << std::endl;
        return false;
    }

    bool i1_found = false;
    bool i2_found = false;

    double dr0 = orbit_radius(params1, 0.0)
               - orbit_radius(params2, 0.0);
    double dr1 = dr0;
    double dr2;
    double theta1=0.0;
    double theta2;
    double dtheta = 2.0*M_PI/N_DR_POINTS;

    for(int i=1;i<=N_DR_POINTS;i++){
        if(i==N_DR_POINTS){
            theta2 = 2.0*M_PI;
            dr2 = dr0;
        }else{
            theta2 = 2.0*M_PI*i/N_DR_POINTS;
            dr2 = orbit_radius(params1, theta2)
              - orbit_radius(params2, theta2);
        }
        if(dr1==0){
            if(dr2>0.0){
                // delta-r positive
                theta_intersect1 = theta1;
                i1_found = true;
            }else if(dr2<0.0){
                // delta-r negative
                theta_intersect2 = theta1;
                i2_found = true;
            }
        }else if(dr1>0.0){
            if(dr2<=0.0){
                // delta-r negative
                double alpha = dr1/(dr1-dr2);
                theta_intersect2 = theta1 + alpha*dtheta;
                i2_found = true;
            }
        }else{
            // dr1 < 0.0
            if(dr2>=0){
                // delta-r positive
                double alpha = dr1/(dr1-dr2);
                theta_intersect1 = theta1 + alpha*dtheta;
                i1_found = true;
            }
        }
        if(i1_found && i2_found){
            break;
        }
        dr1 = dr2;
        theta1 = theta2;
    }

    if(i1_found && i2_found){
        theta_intersect1 = refine_guess(params1,params2,theta_intersect1);
        theta_intersect2 = refine_guess(params1,params2,theta_intersect2);
        return true;
    }
    return false;
}

void orbit_planes_intersect(
        OrbitalParams &params1,
        OrbitalParams &params2,
        double &nu1, // true anomoly of intersection orbit 1
        double &nu2, // true anomoly of intersection orbit 2
        double &theta) // angle between the orbits
{
    // basis for orbit 1
    Eigen::Matrix3d A_basis1 = orbit_basis(params1);

    // basis for orbit 2
    Eigen::Matrix3d A_basis2 = orbit_basis(params2);

    // z-axis of basis is the normal of the orbital plane
    Eigen::Vector3d n1 = A_basis1.col(2);
    Eigen::Vector3d n2 = A_basis2.col(2);

    // cross product of the two normals gives the
    // line vector
    Eigen::Vector3d l = n1.cross(n2);

    // resolve the line vector in each basis
    // space
    Eigen::Vector3d l_1 = A_basis1.transpose()*l;
    Eigen::Vector3d l_2 = A_basis2.transpose()*l;

    nu1 = atan2(l_1(1),l_1(0));
    nu2 = atan2(l_2(1),l_2(0));
    theta = acos(n1.dot(n2));
}


double semi_major(double T, double m1, double m2)
{
    double r = T/2.0/M_PI;
    r*=r;
    r*=G_gravity*(m1+m2);
    r = pow(r,1.0/3.0);
    return r;
}

double period(double a, double m1, double m2)
{
    return 2.0*M_PI*sqrt(a*a*a/G_gravity/(m1+m2));
}

double v_circular(double r, double m1, double m2)
{
    return sqrt(G_gravity*(m1+m2)/r);
}

Eigen::Vector3d direction_from_coord(
    double latitude, double longitude)
{
    double c_lat = cos(latitude);
    double s_lat = sin(latitude);
    double c_lon = cos(longitude);
    double s_lon = sin(longitude);
    Eigen::Vector3d r;
    r << c_lat*c_lon, c_lat*s_lon, s_lat;
    return r;
}

void coord_from_direction(
    Eigen::Vector3d dir,
    double &latitude,
    double &longitude)
{
    Eigen::Vector3d n = dir.normalized();
    double theta = acos(n(2));
    latitude = M_PI/2.0 - theta;
    longitude = atan2(n(1), n(0));
}

Eigen::Matrix3d vehicle_basis(
    Eigen::Vector3d &r, // location of vehicle relative to orbited body
    Eigen::Vector3d &v) // velocity of vehicle relative to orbited body

{
    Eigen::Vector3d h = r.cross(v);
    Eigen::Vector3d n_h = h.normalized();
    Eigen::Vector3d n_r = r.normalized();
    Eigen::Vector3d n_theta = n_h.cross(n_r);
    Eigen::Matrix3d A;
    A.col(0) = n_theta;
    A.col(1) = n_h;
    A.col(2) = n_r;
    return A;
}

double semi_major_from_periapsis_and_apoapsis(
    double r_periapsis, double r_apoapsis)
{
    return (r_periapsis + r_apoapsis)/2.0;
}

double eccentricity_from_periapsis_and_apoapsis(
    double r_periapsis, double r_apoapsis)
{
    double e = (r_apoapsis - r_periapsis)/(r_apoapsis + r_periapsis);
    return fabs(e);
}

double nu_from_r( // true anomoly from a location on the orbit
        OrbitalParams params, // orbital paramaters
        Eigen::Vector3d r) // r, location, assumed to be in the plane
{
    Eigen::Matrix3d A = orbit_basis(params);
    Eigen::Vector3d r_b = A.transpose()*r;
    return atan2(r_b(1), r_b(0));
}

Eigen::Vector3d n_r_earth_frame(double latitude, double longitude)
{
    return Eigen::Vector3d(cos(latitude)*cos(longitude),
                           cos(latitude)*sin(longitude),
                           sin(latitude));
}

Eigen::Vector3d orbit_n_plane(OrbitalParams const &params)
{
    return Eigen::Vector3d(
                sin(params.Omega)*sin(params.i),
                -cos(params.Omega)*sin(params.i),
                cos(params.i));
}

Eigen::Vector3d orbit_ascending_node(OrbitalParams const &params)
{
    return Eigen::Vector3d(
                cos(params.Omega),
                sin(params.Omega),
                0.0);
}

bool orbit_ring_intersect(
        OrbitalParams const &params,
        Eigen::Vector3d const &rs,
        Eigen::Vector3d &ri1,
        Eigen::Vector3d &ri2)
{
    double mag_rs = rs.norm();
    double latitude = M_PI/2.0 - acos(rs(2)/mag_rs);

    if(params.i<fabs(latitude) || params.i>(M_PI-fabs(latitude))){
        return false;
    }

    Eigen::Vector3d n = orbit_n_plane(params);
    double Rxy = sqrt(rs(0)*rs(0) + rs(1)*rs(1));
    if(fabs(n(0))>fabs(n(1))){
        // Nx is greater
        double A = n(0)*n(0) + n(1)*n(1);
        double B = 2*n(1)*n(2)*rs(2);
        double C = n(2)*n(2)*rs(2)*rs(2) - Rxy*Rxy*n(0)*n(0);
        double sqrt_disc = sqrt(B*B - 4*A*C);
        double yc1 = (-B + sqrt_disc)/2/A;
        double yc2 = (-B - sqrt_disc)/2/A;
        double xc1 = -(n(1)*yc1 + n(2)*rs(2))/n(0);
        double xc2 = -(n(1)*yc2 + n(2)*rs(2))/n(0);
        ri1 << xc1, yc1, rs(2);
        ri2 << xc2, yc2, rs(2);
    }else{
        // Ny is greater
        double A = n(0)*n(0) + n(1)*n(1);
        double B = 2*n(0)*n(2)*rs(2);
        double C = n(2)*n(2)*rs(2)*rs(2) - Rxy*Rxy*n(1)*n(1);
        double sqrt_disc = sqrt(B*B - 4*A*C);
        double xc1 = (-B + sqrt_disc)/2/A;
        double xc2 = (-B - sqrt_disc)/2/A;
        double yc1 = -(n(0)*xc1 + n(2)*rs(2))/n(1);
        double yc2 = -(n(0)*xc2 + n(2)*rs(2))/n(1);
        ri1 << xc1, yc1, rs(2);
        ri2 << xc2, yc2, rs(2);
    }
    Eigen::Vector3d na = orbit_ascending_node(params);
    if(na.dot(ri1)<0.0){
        Eigen::Vector3d t = ri1;
        ri1 = ri2;
        ri2 = t;
    }

    return true;
}

double r_perigee_1(
        double r1,
        double ra,
        double theta,
        double alpha)
{
    return -r1*ra*(1+cos(alpha+theta))/
            (r1*(1-cos(alpha+theta))-2.0*ra);
}

double f_uneven_alpha(
        double r1,
        double r2,
        double ra,
        double theta,
        double alpha)
{
    return -4.0*r1*r2*sin(alpha)*sin(theta)
            -2.0*r1*ra*(1+cos(alpha+theta))
            +2.0*r2*ra*(1+cos(alpha-theta));
}

double df_uneven_alpha(
        double r1,
        double r2,
        double ra,
        double theta,
        double alpha)
{
    return -4.0*r1*r2*cos(alpha)*sin(theta)
            +2.0*r1*ra*sin(alpha+theta)
            -2.0*r2*ra*sin(alpha-theta);
}

double f_uneven_solve(
        double r1,
        double r2,
        double ra,
        double theta)
{
    double alpha = 0.0;
    if(r1==r2){
        return alpha;
    }

    while(1)
    {
        double y = f_uneven_alpha(r1,r2,ra,theta,alpha);
        double dy = df_uneven_alpha(r1,r2,ra,theta,alpha);
        double d_alpha = y/dy;
        if(fabs(d_alpha)<1e-12){
            break;
        }
        alpha -= d_alpha;
    }

    return alpha;
}

void ballistic_orbit_params(
        OrbitalParams &params,
        Eigen::Vector3d const &r_start,
        Eigen::Vector3d const &r_end,
        double radius_apogee)
{
    Eigen::Vector3d n0 = r_start.normalized();
    Eigen::Vector3d n1 = r_end.normalized();
    double phi = acos(n0.dot(n1));
    double theta = M_PI - phi/2;

    double r1 = r_start.norm();
    double r2 = r_end.norm();
    double alpha = f_uneven_solve(r1,r2,radius_apogee,theta);
    double radius_perigee = r_perigee_1(r1, radius_apogee, theta, alpha);
    params.a = (radius_apogee + radius_perigee)/2.0;
    params.e = (radius_apogee - radius_perigee)/(radius_apogee + radius_perigee);
    Eigen::Vector3d n_h = n0.cross(n1).normalized();
    Eigen::Vector3d n_e = (-(n0+n1)).normalized();
    Eigen::Matrix3d A_e = caams::AAA(-alpha,n_h);
    n_e = A_e*n_e;
    orbital_plane_elements(n_h,n_e,params);
    // Eigen::Vector3d z_axis(0,0,1);
    // Eigen::Vector3d n_an = z_axis.cross(n_h).normalized();
    // params.i = acos(n_h(2));
    // params.Omega = acos(n_an(0));
    // if(n_an(1)<0.0) params.Omega*=-1.0;
    // params.omega = acos(n_an.dot(n_e));
    // if(n_e(2)<0.0) params.omega*=-1.0;
    params.nu = theta;
}

Eigen::Matrix3d orbit_basis_90(
        Eigen::Matrix3d A_basis)
{
    Eigen::Matrix3d A_basis_90;
    A_basis_90.col(0) = A_basis.col(2);
    A_basis_90.col(1) =-A_basis.col(1);
    A_basis_90.col(2) = A_basis.col(0);
    return A_basis_90;
}







