#include "icbm_simple.h"
#include "caams.hpp"
#include <cmath>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <boost/math/special_functions/ellint_1.hpp>
#include <boost/math/special_functions/ellint_2.hpp>

using namespace PSOPT;

typedef Eigen::Matrix<adouble, 3, 1> aVector3d;
typedef Eigen::Matrix<adouble, 4, 1> aVector4d;
typedef Eigen::Matrix<adouble, 3, 3> aMatrix3d;
typedef Eigen::Matrix<adouble, 6, 1> aVector6d;
typedef Eigen::Matrix<adouble, 5, 1> aVector5d;
typedef Eigen::Matrix<double, 6, 1> Vector6d;

#define RAD_PER_DEG (M_PI/180.0)

static aVector3d g_accel_ring(ICBMSimpleData *md,
                              aVector3d rp);
static aMatrix3d Ap(aVector4d const & p);
static void rv2os(aVector3d &r, aVector3d &v, double mu,
                  adouble &a, aVector3d &e, aVector3d &h, aVector5d &params);


//////////////////////////////////////////////////////////////////////////
///////////////////  Define the end point (Mayer) cost function //////////
//////////////////////////////////////////////////////////////////////////

static adouble endpoint_cost(adouble* initial_states, adouble* final_states,
                      adouble* parameters,adouble& t0, adouble& tf,
                      adouble* xad, int iphase, Workspace* workspace)
{
    return tf;
}

//////////////////////////////////////////////////////////////////////////
///////////////////  Define the integrand (Lagrange) cost function  //////
//////////////////////////////////////////////////////////////////////////

static adouble integrand_cost(adouble* states, adouble* controls, adouble* parameters,
                       adouble& time, adouble* xad, int iphase, Workspace* workspace)
{
    return  0.0;
}

//////////////////////////////////////////////////////////////////////////
///////////////////  Define the DAE's ////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

static void dae(adouble* derivatives, adouble* path, adouble* states,
         adouble* controls, adouble* parameters, adouble& time,
         adouble* xad, int iphase, Workspace* workspace)
{
    ICBMSimpleData *md = (ICBMSimpleData*) workspace->problem->user_data;
    Eigen::Map<Eigen::Matrix<adouble,6,1>> x(states);
    Eigen::Map<Eigen::Matrix<adouble,6,1>> dx(derivatives);
    aVector3d r_vehicle = x.segment<3>(0);
    aVector3d v_vehicle = x.segment<3>(3);

    adouble rad = sqrt( dot( r_vehicle.data(), r_vehicle.data(), 3)  );
    adouble muoverradcubed = (md->mu)/(pow(rad,3));
    aVector3d ag = -muoverradcubed*r_vehicle;
    aVector3d ar = g_accel_ring(md, r_vehicle);

    if(iphase==1){
        Eigen::Map<Eigen::Matrix<adouble, 3, 1>> u(controls);
        aVector3d at = u*md->thrust;

        dx.segment<3>(0) = v_vehicle;
        dx.segment<3>(3) = ag + ar + at;

        path[0] = dot(u.data(), u.data(), 3);
        path[1] = r_vehicle.norm();

    }
}

////////////////////////////////////////////////////////////////////////////
///////////////////  Define the events function ////////////////////////////
////////////////////////////////////////////////////////////////////////////

static void events(adouble* e, adouble* initial_states, adouble* final_states,
            adouble* parameters,adouble& t0, adouble& tf, adouble* xad,
            int iphase, Workspace* workspace)
{
    ICBMSimpleData *md = (ICBMSimpleData*) workspace->problem->user_data;

    if(iphase==1){
        Eigen::Map<Eigen::Matrix<adouble, 11, 1>> ev(e);
        Eigen::Map<Eigen::Matrix<adouble, 6, 1>> xi(initial_states);
        Eigen::Map<Eigen::Matrix<adouble, 6, 1>> xf(final_states);
        aVector3d r = xf.segment<3>(0);
        aVector3d v = xf.segment<3>(3);
        aVector3d ecc;
        aVector3d h;
        adouble a;
        aVector5d params;
        aMatrix3d A_g2a = md->A_g2a.cast<adouble>();
        aVector3d ra = A_g2a*r;
        aVector3d va = A_g2a*v;
        rv2os(ra, va, md->mu, a, ecc, h, params);

        ev.segment<6>(0) = xi;
        ev.segment<5>(6) = params;
    }
}

///////////////////////////////////////////////////////////////////////////
///////////////////  Define the phase linkages function ///////////////////
///////////////////////////////////////////////////////////////////////////

static void linkages( adouble* linkages, adouble* xad, Workspace* workspace)
{
}

////////////////////////////////////////////////////////////////////////////
///////////////////  Define the main routine ///////////////////////////////
////////////////////////////////////////////////////////////////////////////

bool icbm_simple_launch(
        ICBMSimpleData *data,
        Eigen::MatrixXd &controls, // phase 1 controls
        Eigen::MatrixXd &time) // phase 1 time
{
    ////////////////////////////////////////////////////////////////////////////
    ///////////////////  Declare key structures ////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////

    Alg  algorithm;
    Sol  solution;
    Prob problem;

    ////////////////////////////////////////////////////////////////////////////
    ///////////////////  Register problem name  ////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////

    problem.name        		= "Ballistic Missile";
    problem.outfilename                 = "ballistic.txt";
    problem.user_data = (void*)data;

    // Initialize the elliptic integral tables
    data->N_ellip_points = N_ELLIP_POINTS;
    data->k.resize(1, N_ELLIP_POINTS);
    data->ellint_1.resize(1, N_ELLIP_POINTS);
    data->ellint_2.resize(1, N_ELLIP_POINTS);

    double alpha_min = 1e-12;
    double alpha0 = pow(alpha_min, 1.0/(N_ELLIP_POINTS-1));
    for(int i=0;i<N_ELLIP_POINTS;i++){
        double k = 1.0 - pow(alpha0, i);
        data->k(i) = k;
        data->ellint_1(i) = boost::math::ellint_1(k);
        data->ellint_2(i) = boost::math::ellint_2(k);
    }

    ////////////////////////////////////////////////////////////////////////////
    ////////////  Define problem level constants & do level 1 setup ////////////
    ////////////////////////////////////////////////////////////////////////////

    problem.nphases   			= 1;
    problem.nlinkages           = 0;

    psopt_level1_setup(problem);

    /////////////////////////////////////////////////////////////////////////////
    /////////   Define phase related information & do level 2 setup  ////////////
    /////////////////////////////////////////////////////////////////////////////

    problem.phases(1).nstates   = 6;
    problem.phases(1).ncontrols = 3;
    problem.phases(1).nevents   = 11;
    problem.phases(1).npath     = 2;
    problem.phases(1).nodes    = (RowVectorXi(1) << 25).finished();

    psopt_level2_setup(problem, algorithm);

    ////////////////////////////////////////////////////////////////////////////
    ///////////////////  Enter problem bounds information //////////////////////
    ////////////////////////////////////////////////////////////////////////////

    double r_max = data->r_vehicle0.norm()*1000.0;
    double v_max = 20e3;
    double r_min = -r_max;
    double v_min = -v_max;

    for(int iphase=1; iphase<=1; iphase++)
    {
        problem.phases(iphase).bounds.lower.states(0) = r_min;
        problem.phases(iphase).bounds.lower.states(1) = r_min;
        problem.phases(iphase).bounds.lower.states(2) = r_min;

        problem.phases(iphase).bounds.lower.states(3) = v_min;
        problem.phases(iphase).bounds.lower.states(4) = v_min;
        problem.phases(iphase).bounds.lower.states(5) = v_min;

        problem.phases(iphase).bounds.upper.states(0) = r_max;
        problem.phases(iphase).bounds.upper.states(1) = r_max;
        problem.phases(iphase).bounds.upper.states(2) = r_max;

        problem.phases(iphase).bounds.upper.states(3) = v_max;
        problem.phases(iphase).bounds.upper.states(4) = v_max;
        problem.phases(iphase).bounds.upper.states(5) = v_max;
    }

    Eigen::Matrix<double, 11, 1> events1upper,events1lower;
    events1upper.segment<3>(0) = data->r_vehicle0;
    events1upper.segment<3>(3) = data->v_vehicle0;
    events1upper.segment<5>(6) <<
            data->params.a,
            data->params.e,
            data->params.i,
            data->params.Omega,
            data->params.omega;

    events1lower.segment<3>(0) = data->r_vehicle0;
    events1lower.segment<3>(3) = data->v_vehicle0;
    events1lower.segment<5>(6) = events1upper.segment<5>(6);

    problem.phases(1).bounds.lower.events = events1lower;
    problem.phases(1).bounds.upper.events = events1upper;

    problem.phases(1).bounds.lower.controls(0) = -1;
    problem.phases(1).bounds.lower.controls(1) = -1;
    problem.phases(1).bounds.lower.controls(2) = -1;

    problem.phases(1).bounds.upper.controls(0) = 1;
    problem.phases(1).bounds.upper.controls(1) = 1;
    problem.phases(1).bounds.upper.controls(2) = 1;

    problem.phases(1).bounds.lower.path(0) = 1;
    problem.phases(1).bounds.lower.path(1) = data->r_vehicle0.norm();
    problem.phases(1).bounds.upper.path(0) = 1;
    problem.phases(1).bounds.upper.path(1) = data->r_vehicle0.norm()*5.0;

    problem.phases(1).bounds.lower.StartTime = 0.0;
    problem.phases(1).bounds.upper.StartTime = 0.0;

    problem.phases(1).bounds.lower.EndTime = data->T1 * 0.2;
    problem.phases(1).bounds.upper.EndTime = data->T1 * 5.0;


    ////////////////////////////////////////////////////////////////////////////
    ///////////////////  Define & register initial guess ///////////////////////
    ////////////////////////////////////////////////////////////////////////////

    int nnodes1 = problem.phases(1).nodes(0);

    MatrixXd states1(6, nnodes1);
    MatrixXd controls1(3, nnodes1);
    MatrixXd time1(1, nnodes1);

    states1.row(0) = linspace(data->r_vehicle0(0), data->r_boost(0), nnodes1);
    states1.row(1) = linspace(data->r_vehicle0(1), data->r_boost(1), nnodes1);
    states1.row(2) = linspace(data->r_vehicle0(2), data->r_boost(2), nnodes1);
    states1.row(3) = linspace(data->v_vehicle0(0), data->v_boost(0), nnodes1);
    states1.row(4) = linspace(data->v_vehicle0(1), data->v_boost(1), nnodes1);
    states1.row(5) = linspace(data->v_vehicle0(2), data->v_boost(2), nnodes1);

    // states1.row(0) = ones(1,nnodes1)*data->r_vehicle0(0);
    // states1.row(1) = ones(1,nnodes1)*data->r_vehicle0(1);
    // states1.row(2) = ones(1,nnodes1)*data->r_vehicle0(2);
    // states1.row(3) = ones(1,nnodes1)*data->v_vehicle0(0);
    // states1.row(4) = ones(1,nnodes1)*data->v_vehicle0(1);
    // states1.row(5) = ones(1,nnodes1)*data->v_vehicle0(2);

    // controls1.row(0) = ones(1, nnodes1);
    // controls1.row(1) = zeros(1, nnodes1);
    // controls1.row(2) = zeros(1, nnodes1);

    Eigen::Vector3d n_r = data->r_vehicle0.normalized();
    controls1.row(0) = ones(1, nnodes1)*n_r(0);
    controls1.row(1) = ones(1, nnodes1)*n_r(1);
    controls1.row(2) = ones(1, nnodes1)*n_r(2);

    time1 = linspace(0.0, data->T1, nnodes1);

    problem.phases(1).guess.states = states1;
    problem.phases(1).guess.controls = controls1;
    problem.phases(1).guess.time = time1;

    ////////////////////////////////////////////////////////////////////////////
    ///////////////////  Register problem functions  ///////////////////////////
    ////////////////////////////////////////////////////////////////////////////

    problem.integrand_cost 	= &integrand_cost;
    problem.endpoint_cost 	= &endpoint_cost;
    problem.dae 		= &dae;
    problem.events 		= &events;
    problem.linkages		= &linkages;

    ////////////////////////////////////////////////////////////////////////////
    ///////////////////  Enter algorithm options  //////////////////////////////
    ////////////////////////////////////////////////////////////////////////////

    algorithm.nlp_iter_max = 200;
    algorithm.nlp_tolerance = 1.e-6;
    algorithm.nlp_method = "IPOPT";
    algorithm.scaling = "automatic";
    algorithm.derivatives = "automatic";
    algorithm.collocation_method = "trapezoidal";
    algorithm.mesh_refinement = "manual";
    algorithm.ode_tolerance = 1.e-6;

    ////////////////////////////////////////////////////////////////////////////
    ///////////////////  Now call PSOPT to solve the problem   //////////////////
    ////////////////////////////////////////////////////////////////////////////

    std::cout << std::setprecision(15);
    std::cout << "a:" << data->params.a << std::endl;
    std::cout << "e:" << data->params.e << std::endl;
    std::cout << "i:" << data->params.i*180/M_PI << std::endl;
    std::cout << "Omega:" << data->params.Omega*180/M_PI << std::endl;
    std::cout << "omega:" << data->params.omega*180/M_PI << std::endl;


    psopt(solution, problem, algorithm);

    ////////////////////////////////////////////////////////////////////////////
    ///////////  Extract relevant variables from solution structure   //////////
    ////////////////////////////////////////////////////////////////////////////

    if(solution.error_flag || solution.nlp_return_code!=Solve_Succeeded)
        return false;

    controls = solution.get_controls_in_phase(1);
    time = solution.get_time_in_phase(1);
    data->T1 = time(time.cols()-1);

    return true;
}

static aVector3d g_accel_ring(
        ICBMSimpleData *md,
        aVector3d rp)
{
    aVector3d rs = rp/md->ring_radius;
    adouble rho = sqrt(rs(0)*rs(0) + rs(1)*rs(1));
    adouble z = rs(2);
    adouble length = 2*M_PI*md->ring_radius;
    adouble lambda = md->ring_mass/length;
    adouble A = -2*md->G*lambda/md->ring_radius;
    adouble B = (1+rho)*(1+rho) + z*z;
    adouble m = 4*rho/B;
    adouble k = sqrt(m);
    //std::cout << std::setprecision(16);
    //std::cout << "k:" << k.value() << std::endl;
    adouble C = sqrt(B)*((rho-1)*(rho-1)+z*z);
    adouble K_m;//boost::math::ellint_1(k);
    adouble E_m;//boost::math::ellint_2(k);
    smooth_linear_interpolation(&K_m, k, md->k, md->ellint_1, md->N_ellip_points);
    smooth_linear_interpolation(&E_m, k, md->k, md->ellint_2, md->N_ellip_points);
    adouble F_rho;
    // if(rho==0.0){
    //     F_rho = 0.0;
    // }else{
        F_rho = A/rho/C*((rho*rho-1-z*z)*E_m +
                         ((1-rho)*(1-rho)+z*z)*K_m);
    // }
    adouble F_z = A*2*z/C*E_m;
    aVector3d n_rho;
    // if(rho==0){
    //     n_rho = aVector3d::Zero();
    // }else{
        n_rho << rs.head<2>(), 0.0;
        n_rho.normalize();
    // }
    aVector3d n_z(0,0,1);
    aVector3d a_ring_p = F_rho*n_rho + F_z*n_z;
    return a_ring_p;
}

static aMatrix3d Ap(aVector4d const & p)
{
    adouble e0_2 = p(0)*p(0);
    adouble e1_2 = p(1)*p(1);
    adouble e2_2 = p(2)*p(2);
    adouble e3_2 = p(3)*p(3);
    adouble t_e0_e1 = 2.0*p(0)*p(1);
    adouble t_e0_e2 = 2.0*p(0)*p(2);
    adouble t_e0_e3 = 2.0*p(0)*p(3);
    adouble t_e1_e2 = 2.0*p(1)*p(2);
    adouble t_e1_e3 = 2.0*p(1)*p(3);
    adouble t_e2_e3 = 2.0*p(2)*p(3);
    aMatrix3d A;
    A(0,0) = e0_2 + e1_2 - e2_2 - e3_2;
    A(1,0) = t_e0_e3 + t_e1_e2;
    A(2,0) = t_e1_e3 - t_e0_e2;
    A(0,1) = t_e1_e2 - t_e0_e3;
    A(1,1) = e0_2 - e1_2 + e2_2 - e3_2;
    A(2,1) = t_e2_e3 + t_e0_e1;
    A(0,2) = t_e1_e3 + t_e0_e2;
    A(1,2) = t_e2_e3 - t_e0_e1;
    A(2,2) = e0_2 - e1_2 - e2_2 + e3_2;
    return A;
}

static void rv2os(aVector3d &r, aVector3d &v, double mu,
                  adouble &a, aVector3d &ev, aVector3d &h, aVector5d &params)
{
    //h = r.cross(v);
    cross(r.data(), v.data(), h.data());
    aVector3d vxh;
    cross(v.data(), h.data(), vxh.data());
    ev = vxh/mu - r.normalized();
    adouble p = h.squaredNorm()/mu;
    adouble e = ev.norm();
    aVector3d n_e = ev.normalized();
    a = p/(1-e*e);
    aVector3d n_h = h.normalized();
    aVector3d n_an;
    aVector3d z_axis(0,0,1);
    cross(z_axis.data(), n_h.data(), n_an.data());
    n_an.normalize();
    adouble i = acos(n_h(2));

    adouble Omega;
    adouble omega;


#define USE_SMOOTH_HEAVISIDE
    double a_eps = 0.001;


#ifndef USE_SMOOTH_HEAVISIDE
    Omega = acos(n_an(0));
    if(n_an(1)<0.0)
        Omega = -Omega;
    omega = acos(n_e.dot(n_an));
    if(n_e(2) < 0.0)
        omega = -omega;
#endif

#ifdef USE_SMOOTH_HEAVISIDE
    Omega =  smooth_heaviside( (n_an(1)+PSOPT_extras::GetEPS()), a_eps )*acos(n_an(0))
                 +smooth_heaviside( -(n_an(1)+PSOPT_extras::GetEPS()), a_eps )*(-acos(n_an(0)));

    omega =  smooth_heaviside( (n_e(2)), a_eps )*acos(n_e.dot(n_an))
                         +smooth_heaviside( -(n_e(2)), a_eps )*(-acos(n_e.dot(n_an)));
#endif

    params << a, e, i, Omega, omega;

    // std::cout << std::setprecision(15);
    // std::cout << "a:" << a.value() << std::endl;
    // std::cout << "e:" << e.value() << std::endl;
    // std::cout << "i:" << i.value()*180.0/M_PI << std::endl;
    // std::cout << "Om:" << Omega.value()*180/M_PI << std::endl;
    // std::cout << "om:" << omega.value()*180/M_PI << std::endl;

    // std::cout << "rv2os: a:" << a.value() << " e:" << mag_e.value() << std::endl;
    // std::cout << "h:" << h(0).value() << "," << h(1).value() << "," << h(2).value() << std::endl;
    // std::cout << "e:" << e(0).value() << "," << e(1).value() << "," << e(2).value() << std::endl;

    // double dummy = 0/0;
}



