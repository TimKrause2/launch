#include "delayed_launch.h"
#include "caams.hpp"
#include <iostream>
#include <math.h>
#include <boost/math/special_functions/ellint_1.hpp>
#include <boost/math/special_functions/ellint_2.hpp>

using namespace PSOPT;

typedef Eigen::Matrix<double, 6, 1> Vector6d;
typedef Eigen::Matrix<adouble, 3, 1> aVector3d;
typedef Eigen::Matrix<adouble, 6, 1> aVector6d;

static void rv2oe(adouble* rv, adouble* vv, double mu, adouble* oe);
static void oe2rv(Vector6d& oe, double mu, Eigen::Vector3d &ri, Eigen::Vector3d &vi);
static aVector3d g_accel_ring(DelayedLaunchData *md,
                              aVector3d rp);

//////////////////////////////////////////////////////////////////////////
///////////////////  Define the end point (Mayer) cost function //////////
//////////////////////////////////////////////////////////////////////////

static adouble endpoint_cost(adouble* initial_states, adouble* final_states,
                      adouble* parameters,adouble& t0, adouble& tf,
                      adouble* xad, int iphase, Workspace* workspace)
{
    if(iphase==3)
        return tf;
    else
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
    DelayedLaunchData *md = (DelayedLaunchData*) workspace->problem->user_data;
    Eigen::Map<Eigen::Matrix<adouble,6,1>> x(states);
    Eigen::Map<Eigen::Matrix<adouble,6,1>> dx(derivatives);
    aVector3d r = x.segment<3>(0);
    aVector3d v = x.segment<3>(3);

    if(iphase==1){
        aVector3d omega = md->omega.cast<adouble>();
        aVector3d vc;
        cross( omega.data(), r.data(), vc.data() );
        aVector3d ac;
        cross( omega.data(), vc.data(), ac.data() );

        dx.segment<3>(0) = vc;
        dx.segment<3>(3) = ac;
    }else if(iphase==2){
        adouble rad = sqrt( dot( r.data(), r.data(), 3)  );
        adouble muoverradcubed = (md->mu)/(pow(rad,3));
        aVector3d ag = -muoverradcubed*r;
        aVector3d n_r = r/rad;
        aVector3d at = n_r*md->thrust;
        aVector3d ar = g_accel_ring(md, r);
        dx.segment<3>(0) = v;
        dx.segment<3>(3) = ag + at + ar;
    }else if(iphase==3){
        Eigen::Map<Eigen::Matrix<adouble,3,1>> u(controls);
        adouble rad = sqrt( dot( r.data(), r.data(), 3)  );
        adouble muoverradcubed = (md->mu)/(pow(rad,3));
        aVector3d ag = -muoverradcubed*r;
        // thrust force in body coordinate system
        aVector3d a_thrust = u * md->thrust;
        aVector3d ar = g_accel_ring(md, r);

        dx.segment<3>(0) = v;
        dx.segment<3>(3) = ag + a_thrust + ar;

        path[0] = dot(u.data(), u.data(), 3);
    }
}

////////////////////////////////////////////////////////////////////////////
///////////////////  Define the events function ////////////////////////////
////////////////////////////////////////////////////////////////////////////

static void events(adouble* e, adouble* initial_states, adouble* final_states,
            adouble* parameters,adouble& t0, adouble& tf, adouble* xad,
            int iphase, Workspace* workspace)
{
    DelayedLaunchData *md = (DelayedLaunchData*) workspace->problem->user_data;

    if(iphase==1){
        Eigen::Map<Eigen::Matrix<adouble, 6, 1>> ev(e);
        Eigen::Map<Eigen::Matrix<adouble, 6, 1>> xi(initial_states);
        ev = xi;
    }else if(iphase == 3){
        Eigen::Map<Eigen::Matrix<adouble, 5, 1>> ev(e);
        Eigen::Map<Eigen::Matrix<adouble, 6, 1>> xf(final_states);
        aVector3d rf = xf.segment<3>(0);
        aVector3d vf = xf.segment<3>(3);
        Eigen::Matrix<adouble, 6, 1> oe;
        rv2oe(rf.data(), vf.data(), md->mu, oe.data());
        ev = oe.head<5>();
    }
}

///////////////////////////////////////////////////////////////////////////
///////////////////  Define the phase linkages function ///////////////////
///////////////////////////////////////////////////////////////////////////

static void linkages( adouble* linkages, adouble* xad, Workspace* workspace)
{
    aVector6d x1f,x2i,x2f,x3i;
    get_final_states( x1f.data(), xad, 1, workspace );
    get_final_states( x2f.data(), xad, 2, workspace );
    get_initial_states( x2i.data(), xad, 2, workspace );
    get_initial_states( x3i.data(), xad, 3, workspace );

    Eigen::Map<Eigen::Matrix<adouble, 12, 1>> l(linkages);
    l.segment<6>(0) = x1f - x2i;
    l.segment<6>(6) = x2f - x3i;
}

////////////////////////////////////////////////////////////////////////////
///////////////////  Define the main routine ///////////////////////////////
////////////////////////////////////////////////////////////////////////////

void delayed_launch(
        DelayedLaunchData *data,
        Eigen::MatrixXd  &controls, // phase 3 controls; thrust vector
        Eigen::MatrixXd  &time) // phase 3 time
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

    problem.name        		= "Delayed Launch";
    problem.outfilename                 = "dlydlaunch.txt";
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

    problem.nphases   			= 3;
    problem.nlinkages           = 12;

    psopt_level1_setup(problem);

    /////////////////////////////////////////////////////////////////////////////
    /////////   Define phase related information & do level 2 setup  ////////////
    /////////////////////////////////////////////////////////////////////////////

    problem.phases(1).nstates   = 6;
    problem.phases(1).ncontrols = 0;
    problem.phases(1).nevents   = 6;
    problem.phases(1).npath     = 0;
    problem.phases(1).nodes    = (RowVectorXi(2) << 20, 50).finished();

    problem.phases(2).nstates   = 6;
    problem.phases(2).ncontrols = 0;
    problem.phases(2).nevents   = 0;
    problem.phases(2).npath     = 0;
    problem.phases(2).nodes    = (RowVectorXi(2) << 20, 50).finished();

    problem.phases(3).nstates   = 6;
    problem.phases(3).ncontrols = 3;
    problem.phases(3).nevents   = 5;
    problem.phases(3).npath     = 1;
    problem.phases(3).nodes    = (RowVectorXi(2) << 20, 50).finished();

    problem.bounds.lower.linkage = Eigen::VectorXd::Zero(12);
    problem.bounds.upper.linkage = Eigen::VectorXd::Zero(12);

    psopt_level2_setup(problem, algorithm);

    ////////////////////////////////////////////////////////////////////////////
    ///////////////////  Enter problem bounds information //////////////////////
    ////////////////////////////////////////////////////////////////////////////

    Vector6d oe;
    oe(0) = data->params_final.a;
    oe(1) = data->params_final.e;
    oe(2) = data->params_final.i;
    oe(3) = data->params_final.Omega;
    oe(4) = data->params_final.omega;
    oe(5) = data->params_final.nu;

    Eigen::Vector3d r3f;
    Eigen::Vector3d v3f;
    oe2rv(oe, data->mu, r3f, v3f);

    aVector3d arg = r3f.cast<adouble>();
    aVector3d avg = v3f.cast<adouble>();
    adouble aoe[6];

    rv2oe(arg.data(), avg.data(), data->mu, aoe);

    std::cout << "a=" << aoe[0].value() << std::endl;
    std::cout << "e=" << aoe[1].value() << std::endl;
    std::cout << "i=" << aoe[2].value()*180/M_PI << std::endl;
    std::cout << "Om=" << aoe[3].value()*180/M_PI << std::endl;
    std::cout << "om=" << aoe[4].value()*180/M_PI << std::endl;
    std::cout << "nu=" << aoe[5].value()*180/M_PI << std::endl;

    double r_max = r3f.norm()*2;
    double v_max = v3f.norm()*2;
    double r_min = -r_max;
    double v_min = -v_max;

    for(int iphase=1; iphase<=3; iphase++)
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

    Vector6d events1;
    events1.segment<3>(0) = data->r1i;
    events1.segment<3>(3) = data->v1i;

    problem.phases(1).bounds.lower.events = events1;
    problem.phases(1).bounds.upper.events = events1;

    problem.phases(3).bounds.lower.events = oe.head<5>();
    problem.phases(3).bounds.upper.events = oe.head<5>();

    problem.phases(3).bounds.lower.controls(0) = -1;
    problem.phases(3).bounds.lower.controls(1) = -1;
    problem.phases(3).bounds.lower.controls(2) = -1;

    problem.phases(3).bounds.upper.controls(0) = 1;
    problem.phases(3).bounds.upper.controls(1) = 1;
    problem.phases(3).bounds.upper.controls(2) = 1;

    problem.phases(3).bounds.lower.path(0) = 1;
    problem.phases(3).bounds.upper.path(0) = 1;

    problem.phases(1).bounds.lower.StartTime = 0.0;
    problem.phases(1).bounds.upper.StartTime = 0.0;

    problem.phases(1).bounds.lower.EndTime = data->T1 * 0.1;
    problem.phases(1).bounds.upper.EndTime = data->T1;

    problem.phases(2).bounds.lower.StartTime = 0.0;
    problem.phases(2).bounds.upper.StartTime = 0.0;

    problem.phases(2).bounds.lower.EndTime = data->T2;
    problem.phases(2).bounds.upper.EndTime = data->T2;

    problem.phases(3).bounds.lower.StartTime = 0.0;
    problem.phases(3).bounds.upper.StartTime = 0.0;

    problem.phases(3).bounds.lower.EndTime = data->T3*0.5;
    problem.phases(3).bounds.upper.EndTime = data->T3*2.0;

    ////////////////////////////////////////////////////////////////////////////
    ///////////////////  Define & register initial guess ///////////////////////
    ////////////////////////////////////////////////////////////////////////////

    Eigen::Vector3d r2i;
    Eigen::Vector3d v2i;

    Eigen::Matrix3d Arot = caams::AAA(data->mag_omega*data->T1, data->omega);
    r2i = Arot*data->r1i;
    v2i = Arot*data->v1i;

    double ag = data->mu / r2i.squaredNorm();
    double a_net = data->thrust - ag;
    Eigen::Vector3d n_r = r2i.normalized();
    double t2 = data->T2*data->T2;
    Eigen::Vector3d dr2 = 0.5*a_net*t2*n_r;
    Eigen::Vector3d dv2 = a_net*data->T2*n_r;
    Eigen::Vector3d r3i = r2i + dr2;
    Eigen::Vector3d v3i = v2i + dv2;

    int nnodes1 = problem.phases(1).nodes(0);
    int nnodes2 = problem.phases(2).nodes(0);
    int nnodes3 = problem.phases(3).nodes(0);

    MatrixXd states1(6, nnodes1);
    MatrixXd time1(1, nnodes1);
    MatrixXd states2(6, nnodes2);
    MatrixXd time2(1, nnodes2);
    MatrixXd states3(6, nnodes3);
    MatrixXd time3(1, nnodes3);
    MatrixXd controls3(3, nnodes3);

    states1.row(0) = linspace(data->r1i(0), r2i(0), nnodes1);
    states1.row(1) = linspace(data->r1i(1), r2i(1), nnodes1);
    states1.row(2) = linspace(data->r1i(2), r2i(2), nnodes1);
    states1.row(3) = linspace(data->v1i(0), v2i(0), nnodes1);
    states1.row(4) = linspace(data->v1i(1), v2i(1), nnodes1);
    states1.row(5) = linspace(data->v1i(2), v2i(2), nnodes1);

    states2.row(0) = linspace(r2i(0), r3i(0), nnodes2);
    states2.row(1) = linspace(r2i(1), r3i(1), nnodes2);
    states2.row(2) = linspace(r2i(2), r3i(2), nnodes2);
    states2.row(3) = linspace(v2i(0), v3i(0), nnodes2);
    states2.row(4) = linspace(v2i(1), v3i(1), nnodes2);
    states2.row(5) = linspace(v2i(2), v3i(2), nnodes2);

    states3.row(0) = linspace(r3i(0), r3f(0), nnodes3);
    states3.row(1) = linspace(r3i(1), r3f(1), nnodes3);
    states3.row(2) = linspace(r3i(2), r3f(2), nnodes3);
    states3.row(3) = linspace(v3i(0), v3f(0), nnodes3);
    states3.row(4) = linspace(v3i(1), v3f(1), nnodes3);
    states3.row(5) = linspace(v3i(2), v3f(2), nnodes3);

    time1 = linspace(0.0, data->T1, nnodes1);
    time2 = linspace(0.0, data->T2, nnodes2);
    time3 = linspace(0.0, data->T3, nnodes3);

    controls3.row(0) = ones(1, nnodes3);
    controls3.row(1) = zeros(1, nnodes3);
    controls3.row(2) = zeros(1, nnodes3);

    problem.phases(1).guess.states = states1;
    problem.phases(1).guess.time = time1;

    problem.phases(2).guess.states = states2;
    problem.phases(2).guess.time = time2;

    problem.phases(3).guess.states = states3;
    problem.phases(3).guess.time = time3;
    problem.phases(3).guess.controls = controls3;

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
    algorithm.derivatives = "numerical";
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

    controls = solution.get_controls_in_phase(3);
    time = solution.get_time_in_phase(3);
    time1 = solution.get_time_in_phase(1);
    data->T1 = time1(time1.cols()-1);
    data->T3 = time(time.cols()-1);

    states1 = solution.get_states_in_phase(1);

    MatrixXd r, mag_r;

    r = states1.block(0,0,3,states1.cols());
    double mag_r0 = data->r1i.norm();
    mag_r = sum_columns(elemProduct(r,r)).cwiseSqrt() - mag_r0*ones(1,states1.cols());

    Save(mag_r, "mag_r.dat");

}





//
// define the auxilary functions
//

static void rv2oe(adouble* rv, adouble* vv, double mu, adouble* oe)
{
       int j;

       adouble K[3]; K[0] = 0.0; K[1]=0.0; K[2]=1.0;

       adouble hv[3];
       cross(rv,vv, hv);
       // for(int i=0;i<3;i++)
       //     std::cout << "hv[" << i << "]=" << hv[i].value() << std::endl;

       adouble nv[3];
       cross(K, hv, nv);
       // for(int i=0;i<3;i++)
       //     std::cout << "nv[" << i << "]=" << nv[i].value() << std::endl;

       adouble n  = sqrt(  dot(nv,nv,3) );

       adouble h2 = dot(hv,hv,3);

       adouble v2 = dot(vv,vv,3);

       adouble r         = sqrt(dot(rv,rv,3));

       adouble ev[3];
       for(j=0;j<3;j++){
           ev[j] = 1/mu *( (v2-mu/r)*rv[j] - dot(rv,vv,3)*vv[j] );
           //std::cout << "ev[" << j << "]:" << ev[j].value() << std::endl;
        }

       adouble p         = h2/mu;

    adouble e  = sqrt(dot(ev,ev,3));		// eccentricity
    adouble a  = p/(1-e*e);	  			// semimajor axis
    adouble i  = acos(hv[2]/sqrt(h2));		// inclination



//#define USE_SMOOTH_HEAVISIDE
//        double a_eps = 0.1;

#ifndef USE_SMOOTH_HEAVISIDE
    adouble Om = acos(nv[0]/n);			// RAAN
    if ( nv[1] < -PSOPT_extras::GetEPS() ){		// fix quadrant
        Om = 2*pi-Om;
    }
#endif

#ifdef USE_SMOOTH_HEAVISIDE

        adouble Om =  smooth_heaviside( (nv[1]+PSOPT_extras::GetEPS()), a_eps )*acos(nv[0]/n)
                     +smooth_heaviside( -(nv[1]+PSOPT_extras::GetEPS()), a_eps )*(2*pi-acos(nv[0]/n));
#endif

#ifndef USE_SMOOTH_HEAVISIDE
    adouble om = acos(dot(nv,ev,3)/n/e);		// arg of periapsis
    if ( ev[2] < 0 ) {				// fix quadrant
        om = 2*pi-om;
    }
#endif

#ifdef USE_SMOOTH_HEAVISIDE
        adouble om =  smooth_heaviside( (ev[2]), a_eps )*acos(dot(nv,ev,3)/n/e)
                     +smooth_heaviside( -(ev[2]), a_eps )*(2*pi-acos(dot(nv,ev,3)/n/e));
#endif

#ifndef USE_SMOOTH_HEAVISIDE
    adouble nu = acos(dot(ev,rv,3)/e/r);		// true anomaly
    if ( dot(rv,vv,3) < 0 ) {			// fix quadrant
        nu = 2*pi-nu;
    }
#endif

#ifdef USE_SMOOTH_HEAVISIDE
       adouble nu =  smooth_heaviside( dot(rv,vv,3), a_eps )*acos(dot(ev,rv,3)/e/r)
                     +smooth_heaviside( -dot(rv,vv,3), a_eps )*(2*pi-acos(dot(ev,rv,3)/e/r));
#endif

        oe[0] = a;
        oe[1] = e;
        oe[2] = i;
        oe[3] = Om;
        oe[4] = om;
        oe[5] = nu;

        return;
}

static void oe2rv(Vector6d& oe, double mu, Eigen::Vector3d &ri, Eigen::Vector3d &vi)
{
    double a=oe(0), e=oe(1), i=oe(2), Om=oe(3), om=oe(4), nu=oe(5);
    double p = a*(1-e*e);
    double r = p/(1+e*cos(nu));
    MatrixXd rv(3,1);
        rv(0) = r*cos(nu);
        rv(1) = r*sin(nu);
        rv(2) = 0.0;

    MatrixXd vv(3,1);

        vv(0) = -sin(nu);
        vv(1) = e+cos(nu);
        vv(2) = 0.0;
        vv    *= sqrt(mu/p);

    double cO = cos(Om),  sO = sin(Om);
    double co = cos(om),  so = sin(om);
    double ci = cos(i),   si = sin(i);

    MatrixXd R(3,3);
   R(0,0)=  cO*co-sO*so*ci;  R(0,1)=  -cO*so-sO*co*ci; R(0,2)=  sO*si;
    R(1,0)=	 sO*co+cO*so*ci; R(1,1)=  -sO*so+cO*co*ci; R(1,2)=-cO*si;
    R(2,0)=	  so*si;         R(2,1)=    co*si;         R(2,2)=  ci;

    ri = R*rv;
    vi = R*vv;

        return;
}

static aVector3d g_accel_ring(
        DelayedLaunchData *md,
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
    spline_interpolation(&K_m, k, md->k, md->ellint_1, md->N_ellip_points);
    spline_interpolation(&E_m, k, md->k, md->ellint_2, md->N_ellip_points);
    adouble F_rho;
    if(rho==0.0){
        F_rho = 0.0;
    }else{
        F_rho = A/rho/C*((rho*rho-1-z*z)*E_m +
                         ((1-rho)*(1-rho)+z*z)*K_m);
    }
    adouble F_z = A*2*z/C*E_m;
    aVector3d n_rho;
    if(rho==0){
        n_rho = aVector3d::Zero();
    }else{
        n_rho << rs.head<2>(), 0.0;
        n_rho.normalize();
    }
    aVector3d n_z(0,0,1);
    aVector3d a_ring_p = F_rho*n_rho + F_z*n_z;
    return a_ring_p;
}
