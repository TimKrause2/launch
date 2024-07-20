#include "thrust_vector.h"
#include <iostream>
#include <math.h>

using namespace PSOPT;

typedef Eigen::Matrix<adouble, 3, 1> aVector3d;
typedef Eigen::Matrix<adouble, 4, 1> aVector4d;
typedef Eigen::Matrix<adouble, 3, 3> aMatrix3d;
typedef Eigen::Matrix<adouble, 3, 4> aMatrix3x4d;
typedef Eigen::Matrix<adouble, 4, 4> aMatrix4d;
typedef Eigen::Matrix<double, 6, 1>  Vector6d;



void rv2oe(adouble* rv, adouble* vv, double mu, adouble* oe);
void oe2rv(Vector6d& oe, double mu, Eigen::Vector3d &ri, Eigen::Vector3d &vi);

aMatrix3d Ap(aVector4d const &p);
aMatrix3x4d L(aVector4d const &p);

#define RAD_PER_DEG (M_PI/180.0)

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
    ThrustVectorData *md = (ThrustVectorData*) workspace->problem->user_data;
    Eigen::Map<Eigen::Matrix<adouble,6,1>> x(states);
    Eigen::Map<Eigen::Matrix<adouble,6,1>> dx(derivatives);
    Eigen::Map<Eigen::Matrix<adouble,3,1>> u(controls);
    aVector3d r = x.segment<3>(0);
    aVector3d v = x.segment<3>(3);

    adouble rad = sqrt( dot( r.data(), r.data(), 3)  );

    adouble muoverradcubed = (md->mu)/(pow(rad,3));
    aVector3d ag = -muoverradcubed*r;

    // thrust force in body coordinate system
    aVector3d a_thrust = u * md->thrust;

    dx.segment<3>(0) = v;
    dx.segment<3>(3) = ag + a_thrust;

    path[0] = dot(u.data(), u.data(), 3);

}

////////////////////////////////////////////////////////////////////////////
///////////////////  Define the events function ////////////////////////////
////////////////////////////////////////////////////////////////////////////

static void events(adouble* e, adouble* initial_states, adouble* final_states,
            adouble* parameters,adouble& t0, adouble& tf, adouble* xad,
            int iphase, Workspace* workspace)
{
    ThrustVectorData *md = (ThrustVectorData*) workspace->problem->user_data;
    Eigen::Map<Eigen::Matrix<adouble, 11, 1>> ev(e);
    Eigen::Map<Eigen::Matrix<adouble, 6, 1>> xi(initial_states);
    Eigen::Map<Eigen::Matrix<adouble, 6, 1>> xf(final_states);
    aVector3d rf = xf.segment<3>(0);
    aVector3d vf = xf.segment<3>(3);
    Eigen::Matrix<adouble, 6, 1> oe;
    rv2oe(rf.data(), vf.data(), md->mu, oe.data());

    ev.segment<6>(0) = xi;
    ev.segment<5>(6) = oe.head<5>();

    // for(int i=0;i<11;i++){
    //     std::cout << "e[" << i << "]:" << ev[i].value() << std::endl;
    // }


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

void thrust_vector_simple(
        Eigen::Vector3d const &r0,
        Eigen::Vector3d const &v0,
        OrbitalParams const &oe_final,
        ThrustVectorData *data,
        Eigen::MatrixXd control,
        Eigen::MatrixXd time)
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

    problem.name        		= "Thrust Vector Launch";
    problem.outfilename                 = "tvlaunch.txt";

////////////////////////////////////////////////////////////////////////////
///////////////////  Declare an instance of Constants structure /////////////
////////////////////////////////////////////////////////////////////////////

    problem.user_data = (void*)data;



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
    problem.phases(1).npath     = 1;
    problem.phases(1).nodes    << 50;

    psopt_level2_setup(problem, algorithm);

    Vector6d oe;
    oe(0) = oe_final.a;
    oe(1) = oe_final.e;
    oe(2) = oe_final.i;
    oe(3) = oe_final.Omega;
    oe(4) = oe_final.omega;
    oe(5) = oe_final.nu;

    Eigen::Vector3d rg;
    Eigen::Vector3d vg;
    oe2rv(oe, data->mu, rg, vg);

    aVector3d arg = rg.cast<adouble>();
    aVector3d avg = vg.cast<adouble>();
    adouble aoe[6];

    rv2oe(arg.data(), avg.data(), data->mu, aoe);

    std::cout << "a=" << aoe[0].value() << std::endl;
    std::cout << "e=" << aoe[1].value() << std::endl;
    std::cout << "i=" << aoe[2].value()*180/M_PI << std::endl;
    std::cout << "Om=" << aoe[3].value()*180/M_PI << std::endl;
    std::cout << "om=" << aoe[4].value()*180/M_PI << std::endl;
    std::cout << "nu=" << aoe[5].value()*180/M_PI << std::endl;

    double r_max = rg.norm()*2;
    double r_min = -r_max;

    double v_max = vg.norm()*2;
    double v_min = -v_max;

    double pdot_max = 0.5 * 2 * pi;
    double pdot_min = -pdot_max;


    ////////////////////////////////////////////////////////////////////////////
    ///////////////////  Enter problem bounds information //////////////////////
    ////////////////////////////////////////////////////////////////////////////

    // BOUNDS FOR PHASE 1
    problem.phases(1).bounds.lower.states(0) = r_min;
    problem.phases(1).bounds.lower.states(1) = r_min;
    problem.phases(1).bounds.lower.states(2) = r_min;

    problem.phases(1).bounds.lower.states(3) = v_min;
    problem.phases(1).bounds.lower.states(4) = v_min;
    problem.phases(1).bounds.lower.states(5) = v_min;

    problem.phases(1).bounds.upper.states(0) = r_max;
    problem.phases(1).bounds.upper.states(1) = r_max;
    problem.phases(1).bounds.upper.states(2) = r_max;

    problem.phases(1).bounds.upper.states(3) = v_max;
    problem.phases(1).bounds.upper.states(4) = v_max;
    problem.phases(1).bounds.upper.states(5) = v_max;

    Eigen::Matrix<double, 11, 1> events_init;
    events_init.segment<3>(0) = r0;
    events_init.segment<3>(3) = v0;
    events_init.segment<5>(6) = oe.head<5>();

    problem.phases(1).bounds.lower.events = events_init;
    problem.phases(1).bounds.upper.events = events_init;

    problem.phases(1).bounds.lower.controls(0) = -1;
    problem.phases(1).bounds.lower.controls(1) = -1;
    problem.phases(1).bounds.lower.controls(2) = -1;
    problem.phases(1).bounds.upper.controls(0) = 1;
    problem.phases(1).bounds.upper.controls(1) = 1;
    problem.phases(1).bounds.upper.controls(2) = 1;

    problem.phases(1).bounds.lower.path(0) = 1.0;
    problem.phases(1).bounds.upper.path(0) = 1.0;

    problem.phases(1).bounds.lower.StartTime = 0.0;
    problem.phases(1).bounds.upper.StartTime = 0.0;

    problem.phases(1).bounds.lower.EndTime = 300.0;
    problem.phases(1).bounds.upper.EndTime = 1000.0;

    ////////////////////////////////////////////////////////////////////////////
    ///////////////////  Define & register initial guess ///////////////////////
    ////////////////////////////////////////////////////////////////////////////

    int nnodes = problem.phases(1).nodes(0);

    MatrixXd state_guess(6,nnodes);
    state_guess.row(0) = linspace(r0(0), rg(0), nnodes);
    state_guess.row(1) = linspace(r0(1), rg(1), nnodes);
    state_guess.row(2) = linspace(r0(2), rg(2), nnodes);
    state_guess.row(3) = linspace(v0(0), vg(0), nnodes);
    state_guess.row(4) = linspace(v0(1), vg(1), nnodes);
    state_guess.row(5) = linspace(v0(2), vg(2), nnodes);

    Eigen::MatrixXd control_guess(3,nnodes);
    control_guess.row(0) = ones(1,nnodes);
    control_guess.row(1) = zeros(1,nnodes);
    control_guess.row(2) = zeros(1,nnodes);

    Eigen::MatrixXd time_guess(1,nnodes);
    time_guess = linspace(0.0, 600.0, nnodes);

    problem.phases(1).guess.states = state_guess;
    problem.phases(1).guess.controls = control_guess;
    problem.phases(1).guess.time = time_guess;

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

    control = solution.get_controls_in_phase(1);
    time = solution.get_time_in_phase(1);
    // Eigen::MatrixXd x = solution.get_states_in_phase(1);
    // Save(x,"tvstates.dat");
    // Eigen::MatrixXd de = solution.get_dual_events_in_phase(1);
    // Save(de,"tvdual_events.dat");
    plot(time, control, problem.name, "time(s)", "u");

}

void rv2oe(adouble* rv, adouble* vv, double mu, adouble* oe)
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

void oe2rv(Vector6d& oe, double mu, Eigen::Vector3d &ri, Eigen::Vector3d &vi)
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





