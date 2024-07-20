#include "ballistic.h"
#include "caams.hpp"
#include "orbit.h"
#include <cmath>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <boost/math/special_functions/ellint_1.hpp>
#include <boost/math/special_functions/ellint_2.hpp>

using namespace PSOPT;

typedef Eigen::Matrix<adouble, 3, 1> aVector3d;
typedef Eigen::Matrix<adouble, 12, 1> aVector12d;
typedef Eigen::Matrix<double, 12, 1> Vector12d;

#define RAD_PER_DEG (M_PI/180.0)

static aVector3d g_accel_ring(BallisticData *md,
                              aVector3d rp);


//////////////////////////////////////////////////////////////////////////
///////////////////  Define the end point (Mayer) cost function //////////
//////////////////////////////////////////////////////////////////////////

static adouble endpoint_cost(adouble* initial_states, adouble* final_states,
                      adouble* parameters,adouble& t0, adouble& tf,
                      adouble* xad, int iphase, Workspace* workspace)
{
    return 0.0;
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
    BallisticData *md = (BallisticData*) workspace->problem->user_data;
    Eigen::Map<Eigen::Matrix<adouble,12,1>> x(states);
    Eigen::Map<Eigen::Matrix<adouble,12,1>> dx(derivatives);
    aVector3d r_vehicle = x.segment<3>(0);
    aVector3d v_vehicle = x.segment<3>(3);
    aVector3d r_target = x.segment<3>(6);
    aVector3d v_target = x.segment<3>(9);
    aVector3d omega = md->omega.cast<adouble>();
    aVector3d v_target_c;
    aVector3d a_target_c;
    cross( omega.data(), r_target.data(), v_target_c.data() );
    cross( omega.data(), v_target_c.data(), a_target_c.data() );
    dx.segment<3>(6) = v_target_c;
    dx.segment<3>(9) = a_target_c;

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

    }else if(iphase==2){
        dx.segment<3>(0) = v_vehicle;
        dx.segment<3>(3) = ag + ar;
    }
}

////////////////////////////////////////////////////////////////////////////
///////////////////  Define the events function ////////////////////////////
////////////////////////////////////////////////////////////////////////////

static void events(adouble* e, adouble* initial_states, adouble* final_states,
            adouble* parameters,adouble& t0, adouble& tf, adouble* xad,
            int iphase, Workspace* workspace)
{
    BallisticData *md = (BallisticData*) workspace->problem->user_data;

    if(iphase==1){
        Eigen::Map<Eigen::Matrix<adouble, 12, 1>> ev(e);
        Eigen::Map<Eigen::Matrix<adouble, 12, 1>> xi(initial_states);
        ev = xi;
    }else if(iphase == 2){
        Eigen::Map<Eigen::Matrix<adouble, 3, 1>> ev(e);
        Eigen::Map<Eigen::Matrix<adouble, 12, 1>> xf(final_states);
        aVector3d r_vehicle = xf.segment<3>(0);
        aVector3d r_target = xf.segment<3>(6);
        ev = r_vehicle - r_target;
    }
}

///////////////////////////////////////////////////////////////////////////
///////////////////  Define the phase linkages function ///////////////////
///////////////////////////////////////////////////////////////////////////

static void linkages( adouble* linkages, adouble* xad, Workspace* workspace)
{
    aVector12d x1f,x2i;
    get_final_states( x1f.data(), xad, 1, workspace );
    get_initial_states( x2i.data(), xad, 2, workspace );

    Eigen::Map<Eigen::Matrix<adouble, 12, 1>> l(linkages);
    l = x1f - x2i;
}

////////////////////////////////////////////////////////////////////////////
///////////////////  Define the main routine ///////////////////////////////
////////////////////////////////////////////////////////////////////////////

void ballistic_launch(
        BallisticData *data,
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

    double alpha_min = 1e-9;
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

    problem.nphases   			= 2;
    problem.nlinkages           = 12;

    psopt_level1_setup(problem);

    /////////////////////////////////////////////////////////////////////////////
    /////////   Define phase related information & do level 2 setup  ////////////
    /////////////////////////////////////////////////////////////////////////////

    problem.phases(1).nstates   = 12;
    problem.phases(1).ncontrols = 3;
    problem.phases(1).nevents   = 12;
    problem.phases(1).npath     = 1;
    problem.phases(1).nodes    = (RowVectorXi(1) << 150).finished();

    problem.phases(2).nstates   = 12;
    problem.phases(2).ncontrols = 0;
    problem.phases(2).nevents   = 3;
    problem.phases(2).npath     = 0;
    problem.phases(2).nodes    = (RowVectorXi(1) << 500).finished();

    problem.bounds.lower.linkage = Eigen::VectorXd::Zero(12);
    problem.bounds.upper.linkage = Eigen::VectorXd::Zero(12);

    psopt_level2_setup(problem, algorithm);

    ////////////////////////////////////////////////////////////////////////////
    ///////////////////  Enter problem bounds information //////////////////////
    ////////////////////////////////////////////////////////////////////////////

    double r_max = max(data->r_vehicle0.norm(), data->r_target0.norm())*10.0;
    double v_max = 20e3;
    double r_min = -r_max;
    double v_min = -v_max;

    for(int iphase=1; iphase<=2; iphase++)
    {
        problem.phases(iphase).bounds.lower.states(0) = r_min;
        problem.phases(iphase).bounds.lower.states(1) = r_min;
        problem.phases(iphase).bounds.lower.states(2) = r_min;
        problem.phases(iphase).bounds.lower.states(6) = r_min;
        problem.phases(iphase).bounds.lower.states(7) = r_min;
        problem.phases(iphase).bounds.lower.states(8) = r_min;

        problem.phases(iphase).bounds.lower.states(3) = v_min;
        problem.phases(iphase).bounds.lower.states(4) = v_min;
        problem.phases(iphase).bounds.lower.states(5) = v_min;
        problem.phases(iphase).bounds.lower.states(9) = v_min;
        problem.phases(iphase).bounds.lower.states(10) = v_min;
        problem.phases(iphase).bounds.lower.states(11) = v_min;

        problem.phases(iphase).bounds.upper.states(0) = r_max;
        problem.phases(iphase).bounds.upper.states(1) = r_max;
        problem.phases(iphase).bounds.upper.states(2) = r_max;
        problem.phases(iphase).bounds.upper.states(6) = r_max;
        problem.phases(iphase).bounds.upper.states(7) = r_max;
        problem.phases(iphase).bounds.upper.states(8) = r_max;

        problem.phases(iphase).bounds.upper.states(3) = v_max;
        problem.phases(iphase).bounds.upper.states(4) = v_max;
        problem.phases(iphase).bounds.upper.states(5) = v_max;
        problem.phases(iphase).bounds.upper.states(9) = v_max;
        problem.phases(iphase).bounds.upper.states(10) = v_max;
        problem.phases(iphase).bounds.upper.states(11) = v_max;
    }

    Vector12d events1;
    events1.segment<3>(0) = data->r_vehicle0;
    events1.segment<3>(3) = data->v_vehicle0;
    events1.segment<3>(6) = data->r_target0;
    events1.segment<3>(9) = data->v_target0;

    problem.phases(1).bounds.lower.events = events1;
    problem.phases(1).bounds.upper.events = events1;

    problem.phases(2).bounds.lower.events = Eigen::Vector3d::Zero();
    problem.phases(2).bounds.upper.events = Eigen::Vector3d::Zero();

    problem.phases(1).bounds.lower.controls(0) = -1;
    problem.phases(1).bounds.lower.controls(1) = -1;
    problem.phases(1).bounds.lower.controls(2) = -1;

    problem.phases(1).bounds.upper.controls(0) = 1;
    problem.phases(1).bounds.upper.controls(1) = 1;
    problem.phases(1).bounds.upper.controls(2) = 1;

    problem.phases(1).bounds.lower.path(0) = 1;
    problem.phases(1).bounds.upper.path(0) = 1;

    problem.phases(1).bounds.lower.StartTime = 0.0;
    problem.phases(1).bounds.upper.StartTime = 0.0;

    problem.phases(1).bounds.lower.EndTime = data->T1;
    problem.phases(1).bounds.upper.EndTime = data->T1;

    problem.phases(2).bounds.lower.StartTime = 0.0;
    problem.phases(2).bounds.upper.StartTime = 0.0;

    problem.phases(2).bounds.lower.EndTime = data->T2 * 0.2;
    problem.phases(2).bounds.upper.EndTime = data->T2 * 3.0;

    ////////////////////////////////////////////////////////////////////////////
    ///////////////////  Define & register initial guess ///////////////////////
    ////////////////////////////////////////////////////////////////////////////

    int nnodes1 = problem.phases(1).nodes(0);
    int nnodes2 = problem.phases(2).nodes(0);

    MatrixXd states1(12, nnodes1);
    MatrixXd controls1(3, nnodes1);
    MatrixXd time1(1, nnodes1);
    MatrixXd states2(12, nnodes2);
    MatrixXd time2(1, nnodes2);

    states1.row(0) = ones(1,nnodes1)*data->r_vehicle0(0);
    states1.row(1) = ones(1,nnodes1)*data->r_vehicle0(1);
    states1.row(2) = ones(1,nnodes1)*data->r_vehicle0(2);
    states1.row(3) = ones(1,nnodes1)*data->v_vehicle0(0);
    states1.row(4) = ones(1,nnodes1)*data->v_vehicle0(1);
    states1.row(5) = ones(1,nnodes1)*data->v_vehicle0(2);
    states1.row(6) = ones(1,nnodes1)*data->r_target0(0);
    states1.row(7) = ones(1,nnodes1)*data->r_target0(1);
    states1.row(8) = ones(1,nnodes1)*data->r_target0(2);
    states1.row(9) = ones(1,nnodes1)*data->v_target0(0);
    states1.row(10) = ones(1,nnodes1)*data->v_target0(1);
    states1.row(11) = ones(1,nnodes1)*data->v_target0(2);

    states2.row(0) = ones(1,nnodes2)*data->r_target0(0);
    states2.row(1) = ones(1,nnodes2)*data->r_target0(1);
    states2.row(2) = ones(1,nnodes2)*data->r_target0(2);
    states2.row(3) = ones(1,nnodes2)*data->v_target0(0);
    states2.row(4) = ones(1,nnodes2)*data->v_target0(1);
    states2.row(5) = ones(1,nnodes2)*data->v_target0(2);
    states2.row(6) = ones(1,nnodes2)*data->r_target0(0);
    states2.row(7) = ones(1,nnodes2)*data->r_target0(1);
    states2.row(8) = ones(1,nnodes2)*data->r_target0(2);
    states2.row(9) = ones(1,nnodes2)*data->v_target0(0);
    states2.row(10) = ones(1,nnodes2)*data->v_target0(1);
    states2.row(11) = ones(1,nnodes2)*data->v_target0(2);

    controls1.row(0) = ones(1, nnodes1);
    controls1.row(1) = zeros(1, nnodes1);
    controls1.row(2) = zeros(1, nnodes1);

    time1 = linspace(0.0, data->T1, nnodes1);
    time2 = linspace(0.0, data->T2, nnodes2);

    problem.phases(1).guess.states = states1;
    problem.phases(1).guess.controls = controls1;
    problem.phases(1).guess.time = time1;

    problem.phases(2).guess.states = states2;
    problem.phases(2).guess.time = time2;


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

    algorithm.nlp_iter_max = 1000;
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

    psopt(solution, problem, algorithm);

    ////////////////////////////////////////////////////////////////////////////
    ///////////  Extract relevant variables from solution structure   //////////
    ////////////////////////////////////////////////////////////////////////////

    controls = solution.get_controls_in_phase(1);
    time = solution.get_time_in_phase(1);

    time2 = solution.get_time_in_phase(2);
    data->T2 = time2(time2.cols()-1);

    states2 = solution.get_states_in_phase(2);
    Vector12d x2f = states2.col(states2.cols()-1);

    Eigen::Vector3d rf_vehicle = x2f.segment<3>(0);
    Eigen::Vector3d rf_target = x2f.segment<3>(6);
    Eigen::Vector3d rf_target_c = caams::AAA(data->mag_omega*(data->T1+data->T2), data->omega)*data->r_target0;
    double latitude;
    double long_vehicle;
    double long_target;
    double long_target_c;
    coord_from_direction(rf_vehicle, latitude, long_vehicle);
    coord_from_direction(rf_target, latitude, long_target);
    coord_from_direction(rf_target_c, latitude, long_target_c);
    std::cout << std::setprecision(15);
    std::cout << "long_vehicle :" << long_vehicle/RAD_PER_DEG << std::endl;
    std::cout << "long_target  :" << long_target/RAD_PER_DEG << std::endl;
    std::cout << "long_target_c:" << long_target_c/RAD_PER_DEG << std::endl;

}

static aVector3d g_accel_ring(
        BallisticData *md,
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
