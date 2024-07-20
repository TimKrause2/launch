#include "two_phase.h"
#include <iostream>
#include <math.h>

using namespace PSOPT;

typedef Eigen::Matrix<adouble, 3, 1> aVector3d;
typedef Eigen::Matrix<adouble, 3, 3> aMatrix3d;

#define RAD_PER_DEG (M_PI/180.0)

struct Data
{
    double mu;
    double thrust;
};

static Data g_data;

//////////////////////////////////////////////////////////////////////////
///////////////////  Define the end point (Mayer) cost function //////////
//////////////////////////////////////////////////////////////////////////

static adouble endpoint_cost(adouble* initial_states, adouble* final_states,
                      adouble* parameters,adouble& t0, adouble& tf,
                      adouble* xad, int iphase, Workspace* workspace)
{

    if (iphase == 2) {
        return tf;
    }
    else {
        return (0);
    }
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
    // Local integers
    int i, j;

    // Extract individual variables
    double mu = g_data.mu;
    double thrust = g_data.thrust;

    adouble p = states[ 0 ];
    adouble f = states[ 1 ];
    adouble g = states[ 2 ];
    adouble h = states[ 3 ];
    adouble k = states[ 4 ];
    adouble L = states[ 5 ];

    adouble* u  = controls;

    // Define some dependent variables

    adouble q      =  1.0 + f*cos(L) + g*sin(L);
    adouble r      =  p/q;
    adouble alpha2 = h*h - k*k;
    adouble X      = sqrt( h*h + k*k );
    adouble s2     = 1 + X*X;

    adouble delta1;
    adouble delta2;
    adouble delta3;


    if(iphase==1){
        delta1 = 0.0;
        delta2 = 0.0;
        delta3 = 0.0;
    }else{
        adouble latitude  = u[ 0 ];
        adouble longitude = u[ 1 ];

        adouble Tvec_tangent = thrust*cos(latitude)*cos(longitude);
        adouble Tvec_normal = thrust*cos(latitude)*sin(longitude);
        adouble Tvec_radial = thrust*sin(latitude);

        delta1= Tvec_radial;
        delta2= Tvec_tangent;
        delta3= Tvec_normal;
    }

    // derivatives

    adouble pdot =                                      2*p/q*sqrt(p/mu)               * delta2;
    adouble fdot =  sqrt(p/mu)*sin(L) * delta1 + sqrt(p/mu)*(1.0/q)*((q+1.0)*cos(L)+f) * delta2
                   -  sqrt(p/mu)*(g/q)*(h*sin(L)-k*cos(L))  * delta3;
    adouble gdot = -sqrt(p/mu)*cos(L) * delta1 + sqrt(p/mu)*(1.0/q)*((q+1.0)*sin(L)+g) * delta2
                   +  sqrt(p/mu)*(f/q)*(h*sin(L)-k*cos(L))  * delta3;
    adouble hdot =  sqrt(p/mu)*s2*cos(L)/(2.0*q)             * delta3;
    adouble kdot =  sqrt(p/mu)*s2*sin(L)/(2.0*q)             * delta3;
    adouble Ldot =  sqrt(p/mu)*(1.0/q)*(h*sin(L)-k*cos(L))* delta3     + sqrt(mu*p)*pow( (q/p),2.);

    derivatives[ 0 ] = pdot;
    derivatives[ 1 ] = fdot;
    derivatives[ 2 ] = gdot;
    derivatives[ 3 ] = hdot;
    derivatives[ 4 ] = kdot;
    derivatives[ 5 ] = Ldot;

}

////////////////////////////////////////////////////////////////////////////
///////////////////  Define the events function ////////////////////////////
////////////////////////////////////////////////////////////////////////////

static void events(adouble* e, adouble* initial_states, adouble* final_states,
            adouble* parameters,adouble& t0, adouble& tf, adouble* xad,
            int iphase, Workspace* workspace)

{
    if(iphase==1){
        Eigen::Map<Eigen::Matrix<adouble, 6, 1>> x_i(initial_states, 6);
        Eigen::Map<Eigen::Matrix<adouble, 6, 1>> e_vec(e, 6);
        e_vec = x_i;
    }else if(iphase==2){
        Eigen::Map<Eigen::Matrix<adouble, 6, 1>> x_f(final_states, 6);
        Eigen::Map<Eigen::Matrix<adouble, 5, 1>> e_vec(e, 5);
        e_vec = x_f.head(5);
    }
}

///////////////////////////////////////////////////////////////////////////
///////////////////  Define the phase linkages function ///////////////////
///////////////////////////////////////////////////////////////////////////

static void linkages( adouble* linkages, adouble* xad, Workspace* workspace)
{
//    int index = 0;
//    auto_link(linkages, &index, xad, 1, 2, workspace);
    adouble xf[6], xi[6];
    get_final_states( xf, xad, 1, workspace );
    get_initial_states( xi, xad, 2, workspace );

    for(int i=0;i<6;i++)
        linkages[i] = xf[i] - xi[i];
}

////////////////////////////////////////////////////////////////////////////
///////////////////  Define the main routine ///////////////////////////////
////////////////////////////////////////////////////////////////////////////

void two_phase_optimize_constant_longitude(
    ModEquOrbitalParams &meparams0,
    ModEquOrbitalParams &meparams1,
    Eigen::MatrixXd &phase1_time, // initial time vector for phase 1 and result
    Eigen::MatrixXd &phase2_time, // initial time vector for phase 2 and result
    Eigen::MatrixXd &phase2_control, // control guess and result on return
    double mu,
    double thrust)
{
    // check the dimensions of the matricies
    if(phase1_time.rows()!=1){
        std::cout << "two_phase_optimize: phase1_time not a row vector." << std::endl;
        return;
    }
    if(phase2_time.rows()!=1){
        std::cout << "two_phase_optimize: phase2_time not a row vector." << std::endl;
        return;
    }
    if(phase2_control.rows()!=2){
        std::cout << "two_phase_optimize: phase2_control doesn't have two rows." << std::endl;
        return;
    }
    if(phase2_control.cols()!=phase2_time.cols()){
        std::cout << "two_phase_optimize: phase2 columns don't match." << std::endl;
        return;
    }

    // initialize the data structure
    g_data.mu = mu;
    g_data.thrust = thrust;

    ////////////////////////////////////////////////////////////////////////////
    ///////////////////  Declare key structures ////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////

    Alg  algorithm;
    Sol  solution;
    Prob problem;

    ////////////////////////////////////////////////////////////////////////////
    ///////////////////  Register problem name  ////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////

    problem.name         = "Coast-Burn Maneuver";
    problem.outfilename  = "two_phase.txt";

    ////////////////////////////////////////////////////////////////////////////
    ////////////  Define problem level constants & do level 1 setup ////////////
    ////////////////////////////////////////////////////////////////////////////

    problem.nphases     = 2;
    problem.nlinkages   = 6;

    psopt_level1_setup(problem);

    /////////////////////////////////////////////////////////////////////////////
    /////////   Define phase related information & do level 2 setup  ////////////
    /////////////////////////////////////////////////////////////////////////////

    problem.phases(1).nstates   	= 6;
    problem.phases(1).ncontrols 	= 0;
    problem.phases(1).nparameters   = 0;
    problem.phases(1).nevents   	= 6;
    problem.phases(1).npath     	= 0;
    problem.phases(1).nodes         << phase1_time.cols();

    problem.phases(2).nstates   	= 6;
    problem.phases(2).ncontrols 	= 2;
    problem.phases(2).nparameters   = 0;
    problem.phases(2).nevents   	= 5;
    problem.phases(2).npath     	= 0;
    problem.phases(2).nodes         << phase2_time.cols();

    problem.bounds.lower.linkage = Eigen::VectorXd::Zero(6);
    problem.bounds.upper.linkage = Eigen::VectorXd::Zero(6);

    psopt_level2_setup(problem, algorithm);

    ////////////////////////////////////////////////////////////////////////////
    ///////////////////  Enter problem bounds information //////////////////////
    ////////////////////////////////////////////////////////////////////////////

    OrbitalParams oparams0;
    OrbitalParams oparams1;
    meparams0.ToOrbitalParams(oparams0);
    meparams1.ToOrbitalParams(oparams1);
    double ecc_max = fmax(oparams0.e, oparams1.e);
    double p_max = fmax(meparams0.p, meparams1.p)*2.0;
    double p_min = fmin(meparams0.p, meparams1.p)*0.5;
    double hk_max = fmax(tan(oparams0.i/2.0), tan(oparams1.i/2.0));
    double T_coast = phase1_time(phase1_time.cols()-1);
    double T_burn = phase2_time(phase2_time.cols()-1);
    double longitude = phase2_control.col(0)(1);

    // BOUNDS FOR PHASE 1


    problem.phases(1).bounds.lower.states(0) = p_min;
    problem.phases(1).bounds.lower.states(1) = -ecc_max;
    problem.phases(1).bounds.lower.states(2) = -ecc_max;
    problem.phases(1).bounds.lower.states(3) = -hk_max;
    problem.phases(1).bounds.lower.states(4) = -hk_max;
    problem.phases(1).bounds.lower.states(5) = -4*pi;


    problem.phases(1).bounds.upper.states(0) = p_max;
    problem.phases(1).bounds.upper.states(1) = ecc_max;
    problem.phases(1).bounds.upper.states(2) = ecc_max;
    problem.phases(1).bounds.upper.states(3) = hk_max;
    problem.phases(1).bounds.upper.states(4) = hk_max;
    problem.phases(1).bounds.upper.states(5) = 6*pi;


    Eigen::VectorXd phase1_events = meparams0.ToVector();
    problem.phases(1).bounds.lower.events  = phase1_events;
    problem.phases(1).bounds.upper.events  = phase1_events;

    problem.phases(1).bounds.lower.StartTime    = 0.0;
    problem.phases(1).bounds.upper.StartTime    = 0.0;

    problem.phases(1).bounds.lower.EndTime      = T_coast*0.9;
    problem.phases(1).bounds.upper.EndTime      = T_coast;

    // BOUNDS FOR PHASE 2

    problem.phases(2).bounds.lower.states =
        problem.phases(1).bounds.lower.states;

    problem.phases(2).bounds.upper.states =
        problem.phases(1).bounds.upper.states;

    problem.phases(2).bounds.lower.controls(0) = -pi;
    problem.phases(2).bounds.upper.controls(0) = pi;
    problem.phases(2).bounds.lower.controls(1) = longitude;
    problem.phases(2).bounds.upper.controls(1) = longitude;

    problem.phases(2).bounds.lower.StartTime    = 0;
    problem.phases(2).bounds.upper.StartTime    = 0;

    problem.phases(2).bounds.lower.EndTime      = T_burn*0.5;
    problem.phases(2).bounds.upper.EndTime      = T_burn*1.2;

    Eigen::VectorXd phase2_events = meparams1.ToVector();
    problem.phases(2).bounds.lower.events  = phase2_events.head(5);
    problem.phases(2).bounds.upper.events  = phase2_events.head(5);

    ////////////////////////////////////////////////////////////////////////////
    ///////////////////  Register problem functions  ///////////////////////////
    ////////////////////////////////////////////////////////////////////////////


    problem.integrand_cost = &integrand_cost;
    problem.endpoint_cost  = &endpoint_cost;
    problem.dae            = &dae;
    problem.events         = &events;
    problem.linkages       = &linkages;

    ////////////////////////////////////////////////////////////////////////////
    ///////////////////  Define & register initial guess ///////////////////////
    ////////////////////////////////////////////////////////////////////////////

    Eigen::MatrixXd state1_ini;
    state1_ini = meparams0.ToVector();
    Eigen::MatrixXd state1_guess;
    Eigen::MatrixXd state1_control;
    Eigen::MatrixXd state1_parameters;
    rk4_propagate( dae, state1_control, phase1_time, state1_ini,
                  state1_parameters, problem, 1, state1_guess, NULL);

    Eigen::MatrixXd state2_ini;
    state2_ini = state1_guess.col(state1_guess.cols()-1);
    Eigen::MatrixXd state2_guess;
    Eigen::MatrixXd state2_parameters;
    rk4_propagate( dae, phase2_control, phase2_time, state2_ini,
                  state2_parameters, problem, 2, state2_guess, NULL);

    problem.phases(1).guess.states = state1_guess;
    problem.phases(1).guess.time = phase1_time;

    problem.phases(2).guess.states = state2_guess;
    problem.phases(2).guess.time = phase2_time;
    problem.phases(2).guess.controls = phase2_control;

    ModEquOrbitalParams meparams_guess1;
    meparams_guess1.FromVector(state2_ini);
    OrbitalParams oparams_guess1;
    meparams_guess1.ToOrbitalParams(oparams_guess1);
    double nu = oparams_guess1.nu/RAD_PER_DEG;


    ModEquOrbitalParams meparams_guess;
    meparams_guess.FromVector(state2_guess.col(state2_guess.cols()-1));
    OrbitalParams oparams_guess;
    meparams_guess.ToOrbitalParams(oparams_guess);

    ////////////////////////////////////////////////////////////////////////////
    ///////////////////  Enter algorithm options  //////////////////////////////
    ////////////////////////////////////////////////////////////////////////////


    algorithm.nlp_iter_max                = 1000;
    algorithm.nlp_tolerance               = 1.e-6;
    algorithm.nlp_method                  = "IPOPT";
    algorithm.scaling                     = "automatic";
    algorithm.derivatives                 = "automatic";
    algorithm.defect_scaling              = "jacobian-based";
    algorithm.jac_sparsity_ratio          =  0.11;
    algorithm.collocation_method          = "trapezoidal";
    //    algorithm.diff_matrix                 = "central-differences";
    algorithm.mesh_refinement             = "automatic";
    algorithm.mr_max_iterations           = 4;
    algorithm.ode_tolerance               = 1.0e-6;

    ////////////////////////////////////////////////////////////////////////////
    ///////////////////  Now call PSOPT to solve the problem   //////////////////
    ////////////////////////////////////////////////////////////////////////////

    psopt(solution, problem, algorithm);

    ////////////////////////////////////////////////////////////////////////////
    ///////////  Extract relevant variables from solution structure   //////////
    ////////////////////////////////////////////////////////////////////////////

    if(!solution.error_flag){
        phase1_time = solution.get_time_in_phase(1);
        phase2_time = solution.get_time_in_phase(2);
        phase2_control = solution.get_controls_in_phase(2);

        Eigen::MatrixXd opt_latitude = phase2_control.row(0)/RAD_PER_DEG;
        Eigen::MatrixXd opt_longitude = phase2_control.row(1)/RAD_PER_DEG;

        plot(phase2_time, opt_latitude, problem.name+": latitude","time (s)","latitude (deg)", "latitude");
        plot(phase2_time, opt_longitude, problem.name+": longitude", "time (s)", "longitude (deg)", "longitude");

    }else{
        std::cout << "Optimizer error!!!\n";
        std::cout << "error_msg:" << solution.error_msg << std::endl;
    }
}

void two_phase_optimize(
    ModEquOrbitalParams &meparams0,
    ModEquOrbitalParams &meparams1,
    Eigen::MatrixXd &phase1_time, // initial time vector for phase 1 and result
    Eigen::MatrixXd &phase2_time, // initial time vector for phase 2 and result
    Eigen::MatrixXd &phase2_control, // control guess and result on return
    double mu,
    double thrust)
{
    // check the dimensions of the matricies
    if(phase1_time.rows()!=1){
        std::cout << "two_phase_optimize: phase1_time not a row vector." << std::endl;
        return;
    }
    if(phase2_time.rows()!=1){
        std::cout << "two_phase_optimize: phase2_time not a row vector." << std::endl;
        return;
    }
    if(phase2_control.rows()!=2){
        std::cout << "two_phase_optimize: phase2_control doesn't have two rows." << std::endl;
        return;
    }
    if(phase2_control.cols()!=phase2_time.cols()){
        std::cout << "two_phase_optimize: phase2 columns don't match." << std::endl;
        return;
    }

    // initialize the data structure
    g_data.mu = mu;
    g_data.thrust = thrust;

    ////////////////////////////////////////////////////////////////////////////
    ///////////////////  Declare key structures ////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////

    Alg  algorithm;
    Sol  solution;
    Prob problem;

    ////////////////////////////////////////////////////////////////////////////
    ///////////////////  Register problem name  ////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////

    problem.name         = "Coast-Burn Maneuver";
    problem.outfilename  = "two_phase.txt";

    ////////////////////////////////////////////////////////////////////////////
    ////////////  Define problem level constants & do level 1 setup ////////////
    ////////////////////////////////////////////////////////////////////////////

    problem.nphases     = 2;
    problem.nlinkages   = 6;

    psopt_level1_setup(problem);

    /////////////////////////////////////////////////////////////////////////////
    /////////   Define phase related information & do level 2 setup  ////////////
    /////////////////////////////////////////////////////////////////////////////

    problem.phases(1).nstates   	= 6;
    problem.phases(1).ncontrols 	= 0;
    problem.phases(1).nparameters   = 0;
    problem.phases(1).nevents   	= 6;
    problem.phases(1).npath     	= 0;
    problem.phases(1).nodes         << phase1_time.cols();

    problem.phases(2).nstates   	= 6;
    problem.phases(2).ncontrols 	= 2;
    problem.phases(2).nparameters   = 0;
    problem.phases(2).nevents   	= 5;
    problem.phases(2).npath     	= 0;
    problem.phases(2).nodes         << phase2_time.cols();

    problem.bounds.lower.linkage = Eigen::VectorXd::Zero(6);
    problem.bounds.upper.linkage = Eigen::VectorXd::Zero(6);

    psopt_level2_setup(problem, algorithm);

    ////////////////////////////////////////////////////////////////////////////
    ///////////////////  Enter problem bounds information //////////////////////
    ////////////////////////////////////////////////////////////////////////////

    OrbitalParams oparams0;
    OrbitalParams oparams1;
    meparams0.ToOrbitalParams(oparams0);
    meparams1.ToOrbitalParams(oparams1);
    double ecc_max = fmax(oparams0.e, oparams1.e);
    double p_max = fmax(meparams0.p, meparams1.p)*2.0;
    double p_min = fmin(meparams0.p, meparams1.p)*0.5;
    double hk_max = fmax(tan(oparams0.i/2.0), tan(oparams1.i/2.0));
    double T_coast = phase1_time(phase1_time.cols()-1);
    double T_burn = phase2_time(phase2_time.cols()-1);
    double longitude = phase2_control.col(0)(1);

    // BOUNDS FOR PHASE 1


    problem.phases(1).bounds.lower.states(0) = p_min;
    problem.phases(1).bounds.lower.states(1) = -ecc_max;
    problem.phases(1).bounds.lower.states(2) = -ecc_max;
    problem.phases(1).bounds.lower.states(3) = -hk_max;
    problem.phases(1).bounds.lower.states(4) = -hk_max;
    problem.phases(1).bounds.lower.states(5) = -4*pi;


    problem.phases(1).bounds.upper.states(0) = p_max;
    problem.phases(1).bounds.upper.states(1) = ecc_max;
    problem.phases(1).bounds.upper.states(2) = ecc_max;
    problem.phases(1).bounds.upper.states(3) = hk_max;
    problem.phases(1).bounds.upper.states(4) = hk_max;
    problem.phases(1).bounds.upper.states(5) = 6*pi;


    Eigen::VectorXd phase1_events = meparams0.ToVector();
    problem.phases(1).bounds.lower.events  = phase1_events;
    problem.phases(1).bounds.upper.events  = phase1_events;

    problem.phases(1).bounds.lower.StartTime    = 0.0;
    problem.phases(1).bounds.upper.StartTime    = 0.0;

    problem.phases(1).bounds.lower.EndTime      = T_coast*0.9;
    problem.phases(1).bounds.upper.EndTime      = T_coast;

    // BOUNDS FOR PHASE 2

    problem.phases(2).bounds.lower.states =
        problem.phases(1).bounds.lower.states;

    problem.phases(2).bounds.upper.states =
        problem.phases(1).bounds.upper.states;

    problem.phases(2).bounds.lower.controls(0) = -pi;
    problem.phases(2).bounds.upper.controls(0) = pi;
    problem.phases(2).bounds.lower.controls(1) = -pi;
    problem.phases(2).bounds.upper.controls(1) = pi;

    problem.phases(2).bounds.lower.StartTime    = 0;
    problem.phases(2).bounds.upper.StartTime    = 0;

    problem.phases(2).bounds.lower.EndTime      = T_burn*0.5;
    problem.phases(2).bounds.upper.EndTime      = T_burn*1.2;

    Eigen::VectorXd phase2_events = meparams1.ToVector();
    problem.phases(2).bounds.lower.events  = phase2_events.head(5);
    problem.phases(2).bounds.upper.events  = phase2_events.head(5);

    ////////////////////////////////////////////////////////////////////////////
    ///////////////////  Register problem functions  ///////////////////////////
    ////////////////////////////////////////////////////////////////////////////


    problem.integrand_cost = &integrand_cost;
    problem.endpoint_cost  = &endpoint_cost;
    problem.dae            = &dae;
    problem.events         = &events;
    problem.linkages       = &linkages;

    ////////////////////////////////////////////////////////////////////////////
    ///////////////////  Define & register initial guess ///////////////////////
    ////////////////////////////////////////////////////////////////////////////

    Eigen::MatrixXd state1_ini;
    state1_ini = meparams0.ToVector();
    Eigen::MatrixXd state1_guess;
    Eigen::MatrixXd state1_control;
    Eigen::MatrixXd state1_parameters;
    rk4_propagate( dae, state1_control, phase1_time, state1_ini,
                  state1_parameters, problem, 1, state1_guess, NULL);

    Eigen::MatrixXd state2_ini;
    state2_ini = state1_guess.col(state1_guess.cols()-1);
    Eigen::MatrixXd state2_guess;
    Eigen::MatrixXd state2_parameters;
    rk4_propagate( dae, phase2_control, phase2_time, state2_ini,
                  state2_parameters, problem, 2, state2_guess, NULL);

    problem.phases(1).guess.states = state1_guess;
    problem.phases(1).guess.time = phase1_time;

    problem.phases(2).guess.states = state2_guess;
    problem.phases(2).guess.time = phase2_time;
    problem.phases(2).guess.controls = phase2_control;

    ModEquOrbitalParams meparams_guess1;
    meparams_guess1.FromVector(state2_ini);
    OrbitalParams oparams_guess1;
    meparams_guess1.ToOrbitalParams(oparams_guess1);
    double nu = oparams_guess1.nu/RAD_PER_DEG;


    ModEquOrbitalParams meparams_guess;
    meparams_guess.FromVector(state2_guess.col(state2_guess.cols()-1));
    OrbitalParams oparams_guess;
    meparams_guess.ToOrbitalParams(oparams_guess);

    ////////////////////////////////////////////////////////////////////////////
    ///////////////////  Enter algorithm options  //////////////////////////////
    ////////////////////////////////////////////////////////////////////////////


    algorithm.nlp_iter_max                = 1000;
    algorithm.nlp_tolerance               = 1.e-6;
    algorithm.nlp_method                  = "IPOPT";
    algorithm.scaling                     = "automatic";
    algorithm.derivatives                 = "automatic";
    algorithm.defect_scaling              = "jacobian-based";
    algorithm.jac_sparsity_ratio          =  0.11;
    //algorithm.collocation_method          = "trapezoidal";
    //    algorithm.diff_matrix                 = "central-differences";
    algorithm.mesh_refinement             = "automatic";
    algorithm.mr_max_iterations           = 4;
    algorithm.ode_tolerance               = 1.0e-6;

    ////////////////////////////////////////////////////////////////////////////
    ///////////////////  Now call PSOPT to solve the problem   //////////////////
    ////////////////////////////////////////////////////////////////////////////

    psopt(solution, problem, algorithm);

    ////////////////////////////////////////////////////////////////////////////
    ///////////  Extract relevant variables from solution structure   //////////
    ////////////////////////////////////////////////////////////////////////////

    if(!solution.error_flag){
        phase1_time = solution.get_time_in_phase(1);
        phase2_time = solution.get_time_in_phase(2);
        phase2_control = solution.get_controls_in_phase(2);

        Eigen::MatrixXd opt_latitude = phase2_control.row(0)/RAD_PER_DEG;
        Eigen::MatrixXd opt_longitude = phase2_control.row(1)/RAD_PER_DEG;

        plot(phase2_time, opt_latitude, problem.name+": latitude","time (s)","latitude (deg)", "latitude");
        plot(phase2_time, opt_longitude, problem.name+": longitude", "time (s)", "longitude (deg)", "longitude");

    }else{
        std::cout << "Optimizer error!!!\n";
        std::cout << "error_msg:" << solution.error_msg << std::endl;
    }
}

