#include "two_phase.h"
#include <iostream>
#include <math.h>

using namespace PSOPT;

typedef Eigen::Matrix<adouble, 3, 1> aVector3d;
typedef Eigen::Matrix<adouble, 3, 3> aMatrix3d;

#define RAD_PER_DEG (M_PI/180.0)

struct MyData
{
    double mu;
    double thrust;
    MatrixXd phase2_u;
    MatrixXd phase2_t;
};

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

static adouble endpoint_cost_guess(adouble* initial_states, adouble* final_states,
                             adouble* parameters,adouble& t0, adouble& tf,
                             adouble* xad, int iphase, Workspace* workspace)
{
    return (0);
}

//////////////////////////////////////////////////////////////////////////
///////////////////  Define the integrand (Lagrange) cost function  //////
//////////////////////////////////////////////////////////////////////////

static adouble integrand_cost(adouble* states, adouble* controls, adouble* parameters,
                       adouble& time, adouble* xad, int iphase, Workspace* workspace)
{
    return  0.0;
}

static adouble integrand_cost_guess(adouble* states, adouble* controls, adouble* parameters,
                              adouble& time, adouble* xad, int iphase, Workspace* workspace)
{
    if(iphase==2){
        MatrixXd usm1, usm2;
        adouble us_1, us_2;
        
        MyData *md = (MyData*) workspace->problem->user_data;
        
        adouble u1 = controls[0];
        adouble u2 = controls[1];
        
        usm1 = md->phase2_u.row(0);
        usm2 = md->phase2_u.row(1);
        
        spline_interpolation(&us_1, time, md->phase2_t, usm1, usm1.cols());
        spline_interpolation(&us_2, time, md->phase2_t, usm2, usm2.cols());
        
        return pow(u1-us_1, 2.0) + pow(u2-us_2, 2.0);
        
    }else{
        return 0.0;
    }
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
    MyData *md = (MyData*) workspace->problem->user_data;
    double mu = md->mu;
    double thrust = md->thrust;

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

static void events_guess(adouble* e, adouble* initial_states, adouble* final_states,
                   adouble* parameters,adouble& t0, adouble& tf, adouble* xad,
                   int iphase, Workspace* workspace)

{
    if(iphase==1){
        Eigen::Map<Eigen::Matrix<adouble, 6, 1>> x_i(initial_states, 6);
        Eigen::Map<Eigen::Matrix<adouble, 6, 1>> e_vec(e, 6);
        e_vec = x_i;
    }
}

///////////////////////////////////////////////////////////////////////////
///////////////////  Define the phase linkages function ///////////////////
///////////////////////////////////////////////////////////////////////////

static void linkages( adouble* linkages, adouble* xad, Workspace* workspace)
{
    adouble xf[6], xi[6];
    get_final_states( xf, xad, 1, workspace );
    get_initial_states( xi, xad, 2, workspace );

    for(int i=0;i<6;i++)
        linkages[i] = xf[i] - xi[i];
}

////////////////////////////////////////////////////////////////////////////
///////////////////  Define the main routine ///////////////////////////////
////////////////////////////////////////////////////////////////////////////

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
    MyData md;
    md.mu = mu;
    md.thrust = thrust;
    md.phase2_t = phase2_time;
    md.phase2_u = phase2_control;
    
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
    Eigen::VectorXd phase1_events = meparams0.ToVector();
    Eigen::VectorXd phase2_events = meparams1.ToVector();

    ////////////////////////////////////////////////////////////////////////////
    ///////////////////  Declare key structures ////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////

    Alg  guess_algorithm;
    Sol  guess_solution;
    Prob guess_problem;

    ////////////////////////////////////////////////////////////////////////////
    ///////////////////  Register problem name  ////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////

    guess_problem.name         = "Coast-Burn Maneuver Guess";
    guess_problem.outfilename  = "two_phase_guess.txt";

    ////////////////////////////////////////////////////////////////////////////
    ////////////  Define problem level constants & do level 1 setup ////////////
    ////////////////////////////////////////////////////////////////////////////

    guess_problem.nphases     = 2;
    guess_problem.nlinkages   = 6;

    guess_problem.user_data = (void*)&md;
    
    psopt_level1_setup(guess_problem);

    /////////////////////////////////////////////////////////////////////////////
    /////////   Define phase related information & do level 2 setup  ////////////
    /////////////////////////////////////////////////////////////////////////////

    guess_problem.phases(1).nstates   	= 6;
    guess_problem.phases(1).ncontrols 	= 0;
    guess_problem.phases(1).nparameters   = 0;
    guess_problem.phases(1).nevents   	= 6;
    guess_problem.phases(1).npath     	= 0;
    guess_problem.phases(1).nodes         << phase1_time.cols();

    guess_problem.phases(2).nstates   	= 6;
    guess_problem.phases(2).ncontrols 	= 2;
    guess_problem.phases(2).nparameters   = 0;
    guess_problem.phases(2).nevents   	= 0;
    guess_problem.phases(2).npath     	= 0;
    guess_problem.phases(2).nodes         << phase2_time.cols();

    guess_problem.bounds.lower.linkage = Eigen::VectorXd::Zero(6);
    guess_problem.bounds.upper.linkage = Eigen::VectorXd::Zero(6);

    psopt_level2_setup(guess_problem, guess_algorithm);

    ////////////////////////////////////////////////////////////////////////////
    ///////////////////  Enter problem bounds information //////////////////////
    ////////////////////////////////////////////////////////////////////////////

    // BOUNDS FOR PHASE 1


    guess_problem.phases(1).bounds.lower.states(0) = p_min;
    guess_problem.phases(1).bounds.lower.states(1) = -ecc_max;
    guess_problem.phases(1).bounds.lower.states(2) = -ecc_max;
    guess_problem.phases(1).bounds.lower.states(3) = -hk_max;
    guess_problem.phases(1).bounds.lower.states(4) = -hk_max;
    guess_problem.phases(1).bounds.lower.states(5) = -4*pi;


    guess_problem.phases(1).bounds.upper.states(0) = p_max;
    guess_problem.phases(1).bounds.upper.states(1) = ecc_max;
    guess_problem.phases(1).bounds.upper.states(2) = ecc_max;
    guess_problem.phases(1).bounds.upper.states(3) = hk_max;
    guess_problem.phases(1).bounds.upper.states(4) = hk_max;
    guess_problem.phases(1).bounds.upper.states(5) = 6*pi;


    guess_problem.phases(1).bounds.lower.events  = phase1_events;
    guess_problem.phases(1).bounds.upper.events  = phase1_events;

    guess_problem.phases(1).bounds.lower.StartTime    = 0.0;
    guess_problem.phases(1).bounds.upper.StartTime    = 0.0;

    guess_problem.phases(1).bounds.lower.EndTime      = T_coast;
    guess_problem.phases(1).bounds.upper.EndTime      = T_coast;

    // BOUNDS FOR PHASE 2

    guess_problem.phases(2).bounds.lower.states =
        guess_problem.phases(1).bounds.lower.states;

    guess_problem.phases(2).bounds.upper.states =
        guess_problem.phases(1).bounds.upper.states;

    guess_problem.phases(2).bounds.lower.controls(0) = -pi;
    guess_problem.phases(2).bounds.upper.controls(0) = pi;
    guess_problem.phases(2).bounds.lower.controls(1) = -pi;
    guess_problem.phases(2).bounds.upper.controls(1) = pi;

    guess_problem.phases(2).bounds.lower.StartTime    = 0;
    guess_problem.phases(2).bounds.upper.StartTime    = 0;

    guess_problem.phases(2).bounds.lower.EndTime      = T_burn;
    guess_problem.phases(2).bounds.upper.EndTime      = T_burn;

    ////////////////////////////////////////////////////////////////////////////
    ///////////////////  Register problem functions  ///////////////////////////
    ////////////////////////////////////////////////////////////////////////////


    guess_problem.integrand_cost = &integrand_cost_guess;
    guess_problem.endpoint_cost  = &endpoint_cost_guess;
    guess_problem.dae            = &dae;
    guess_problem.events         = &events_guess;
    guess_problem.linkages       = &linkages;

    ////////////////////////////////////////////////////////////////////////////
    ///////////////////  Define & register initial guess ///////////////////////
    ////////////////////////////////////////////////////////////////////////////

    Eigen::MatrixXd state1_ini;
    state1_ini = meparams0.ToVector();
    Eigen::MatrixXd state1_guess(6,phase1_time.cols());
    for(int i=0;i<6;i++){
        state1_guess.row(i) = ones(1,phase1_time.cols())*state1_ini(i);
    }
    Eigen::MatrixXd state2_guess(6,phase2_time.cols());
    for(int i=0;i<6;i++){
        state2_guess.row(i) = ones(1,phase2_time.cols())*state1_ini(i);
    }
    
    guess_problem.phases(1).guess.states = state1_guess;
    guess_problem.phases(1).guess.time = phase1_time;

    guess_problem.phases(2).guess.states = state2_guess;
    guess_problem.phases(2).guess.time = phase2_time;
    guess_problem.phases(2).guess.controls = phase2_control;
    
    ////////////////////////////////////////////////////////////////////////////
    ///////////////////  Enter algorithm options  //////////////////////////////
    ////////////////////////////////////////////////////////////////////////////


    guess_algorithm.nlp_iter_max                = 1000;
    guess_algorithm.nlp_tolerance               = 1.e-6;
    guess_algorithm.nlp_method                  = "IPOPT";
    guess_algorithm.scaling                     = "automatic";
    guess_algorithm.derivatives                 = "numerical";
    //guess_algorithm.defect_scaling              = "jacobian-based";
    guess_algorithm.jac_sparsity_ratio          =  0.15;
    guess_algorithm.collocation_method          = "trapezoidal";
    //    guess_algorithm.diff_matrix                 = "central-differences";
    guess_algorithm.mesh_refinement             = "manual";
    guess_algorithm.ode_tolerance               = 1.0e-6;

    ////////////////////////////////////////////////////////////////////////////
    ///////////////////  Now call PSOPT to solve the problem   //////////////////
    ////////////////////////////////////////////////////////////////////////////

    psopt(guess_solution, guess_problem, guess_algorithm);

    ////////////////////////////////////////////////////////////////////////////
    ///////////  Extract relevant variables from solution structure   //////////
    ////////////////////////////////////////////////////////////////////////////

    Eigen::MatrixXd main_phase1_state_guess;
    Eigen::MatrixXd main_phase1_time_guess;
    Eigen::MatrixXd main_phase2_state_guess;
    Eigen::MatrixXd main_phase2_time_guess;
    Eigen::MatrixXd main_phase2_control_guess;
    
    if(!guess_solution.error_flag){
        main_phase1_state_guess = guess_solution.get_states_in_phase(1);
        main_phase1_time_guess = guess_solution.get_time_in_phase(1);
        main_phase2_state_guess = guess_solution.get_states_in_phase(2);
        main_phase2_time_guess = guess_solution.get_time_in_phase(2);
        main_phase2_control_guess = guess_solution.get_controls_in_phase(2);
    }else{
        std::cout << "Optimizer error in guess determination!!!\n";
        std::cout << "error_msg:" << guess_solution.error_msg << std::endl;
        return;
    }
    
    
    ////////////////////////////////////////////////////////////////////////////
    ///////////////////  Declare key structures ////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////
    
    Alg  main_algorithm;
    Sol  main_solution;
    Prob main_problem;
    
    ////////////////////////////////////////////////////////////////////////////
    ///////////////////  Register problem name  ////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////
    
    main_problem.name         = "Coast-Burn Maneuver";
    main_problem.outfilename  = "two_phase.txt";
    
    ////////////////////////////////////////////////////////////////////////////
    ////////////  Define problem level constants & do level 1 setup ////////////
    ////////////////////////////////////////////////////////////////////////////
    
    main_problem.nphases     = 2;
    main_problem.nlinkages   = 6;
    
    main_problem.user_data = (void*)&md;
    
    psopt_level1_setup(main_problem);
    
    /////////////////////////////////////////////////////////////////////////////
    /////////   Define phase related information & do level 2 setup  ////////////
    /////////////////////////////////////////////////////////////////////////////
    
    main_problem.phases(1).nstates   	= 6;
    main_problem.phases(1).ncontrols 	= 0;
    main_problem.phases(1).nparameters   = 0;
    main_problem.phases(1).nevents   	= 6;
    main_problem.phases(1).npath     	= 0;
    main_problem.phases(1).nodes         << phase1_time.cols();
    
    main_problem.phases(2).nstates   	= 6;
    main_problem.phases(2).ncontrols 	= 2;
    main_problem.phases(2).nparameters   = 0;
    main_problem.phases(2).nevents   	= 5;
    main_problem.phases(2).npath     	= 0;
    main_problem.phases(2).nodes         << phase2_time.cols();
    
    main_problem.bounds.lower.linkage = Eigen::VectorXd::Zero(6);
    main_problem.bounds.upper.linkage = Eigen::VectorXd::Zero(6);
    
    psopt_level2_setup(main_problem, main_algorithm);
    
    ////////////////////////////////////////////////////////////////////////////
    ///////////////////  Enter problem bounds information //////////////////////
    ////////////////////////////////////////////////////////////////////////////
    
    // BOUNDS FOR PHASE 1
    
    
    main_problem.phases(1).bounds.lower.states(0) = p_min;
    main_problem.phases(1).bounds.lower.states(1) = -ecc_max;
    main_problem.phases(1).bounds.lower.states(2) = -ecc_max;
    main_problem.phases(1).bounds.lower.states(3) = -hk_max;
    main_problem.phases(1).bounds.lower.states(4) = -hk_max;
    main_problem.phases(1).bounds.lower.states(5) = -4*pi;
    
    
    main_problem.phases(1).bounds.upper.states(0) = p_max;
    main_problem.phases(1).bounds.upper.states(1) = ecc_max;
    main_problem.phases(1).bounds.upper.states(2) = ecc_max;
    main_problem.phases(1).bounds.upper.states(3) = hk_max;
    main_problem.phases(1).bounds.upper.states(4) = hk_max;
    main_problem.phases(1).bounds.upper.states(5) = 6*pi;
    
    
    main_problem.phases(1).bounds.lower.events  = phase1_events;
    main_problem.phases(1).bounds.upper.events  = phase1_events;
    
    main_problem.phases(1).bounds.lower.StartTime    = 0.0;
    main_problem.phases(1).bounds.upper.StartTime    = 0.0;
    
    main_problem.phases(1).bounds.lower.EndTime      = T_coast*0.9;
    main_problem.phases(1).bounds.upper.EndTime      = T_coast;
    
    // BOUNDS FOR PHASE 2
    
    main_problem.phases(2).bounds.lower.states =
    main_problem.phases(1).bounds.lower.states;
    
    main_problem.phases(2).bounds.upper.states =
    main_problem.phases(1).bounds.upper.states;
    
    main_problem.phases(2).bounds.lower.controls(0) = -pi;
    main_problem.phases(2).bounds.upper.controls(0) = pi;
    main_problem.phases(2).bounds.lower.controls(1) = -pi;
    main_problem.phases(2).bounds.upper.controls(1) = pi;
    
    main_problem.phases(2).bounds.lower.StartTime    = 0;
    main_problem.phases(2).bounds.upper.StartTime    = 0;
    
    main_problem.phases(2).bounds.lower.EndTime      = T_burn*0.5;
    main_problem.phases(2).bounds.upper.EndTime      = T_burn*1.2;
    
    main_problem.phases(2).bounds.lower.events  = phase2_events.head(5);
    main_problem.phases(2).bounds.upper.events  = phase2_events.head(5);
    
    ////////////////////////////////////////////////////////////////////////////
    ///////////////////  Register problem functions  ///////////////////////////
    ////////////////////////////////////////////////////////////////////////////
    
    
    main_problem.integrand_cost = &integrand_cost;
    main_problem.endpoint_cost  = &endpoint_cost;
    main_problem.dae            = &dae;
    main_problem.events         = &events;
    main_problem.linkages       = &linkages;
    
    ////////////////////////////////////////////////////////////////////////////
    ///////////////////  Define & register initial guess ///////////////////////
    ////////////////////////////////////////////////////////////////////////////
    
    
    main_problem.phases(1).guess.states = main_phase1_state_guess;
    main_problem.phases(1).guess.time = main_phase1_time_guess;
    
    main_problem.phases(2).guess.states = main_phase2_state_guess;
    main_problem.phases(2).guess.time = main_phase2_time_guess;
    main_problem.phases(2).guess.controls = main_phase2_control_guess;
    
    ////////////////////////////////////////////////////////////////////////////
    ///////////////////  Enter algorithm options  //////////////////////////////
    ////////////////////////////////////////////////////////////////////////////
    
    
    main_algorithm.nlp_iter_max                = 1000;
    main_algorithm.nlp_tolerance               = 1.e-6;
    main_algorithm.nlp_method                  = "IPOPT";
    main_algorithm.scaling                     = "automatic";
    main_algorithm.derivatives                 = "automatic";
    main_algorithm.defect_scaling              = "jacobian-based";
    main_algorithm.jac_sparsity_ratio          =  0.15;
    //main_algorithm.collocation_method          = "trapezoidal";
    //    main_algorithm.diff_matrix                 = "central-differences";
    main_algorithm.mesh_refinement             = "automatic";
    main_algorithm.mr_max_iterations           = 8;
    main_algorithm.ode_tolerance               = 1.0e-6;
    
    ////////////////////////////////////////////////////////////////////////////
    ///////////////////  Now call PSOPT to solve the problem   //////////////////
    ////////////////////////////////////////////////////////////////////////////
    
    psopt(main_solution, main_problem, main_algorithm);
    
    ////////////////////////////////////////////////////////////////////////////
    ///////////  Extract relevant variables from solution structure   //////////
    ////////////////////////////////////////////////////////////////////////////
    
    if(!main_solution.error_flag){
        phase1_time = main_solution.get_time_in_phase(1);
        phase2_time = main_solution.get_time_in_phase(2);
        phase2_control = main_solution.get_controls_in_phase(2);
        
        Eigen::MatrixXd opt_latitude = phase2_control.row(0)/RAD_PER_DEG;
        Eigen::MatrixXd opt_longitude = phase2_control.row(1)/RAD_PER_DEG;
        
        plot(phase2_time, opt_latitude, main_problem.name+": latitude","time (s)","latitude (deg)", "latitude");
        plot(phase2_time, opt_longitude, main_problem.name+": longitude", "time (s)", "longitude (deg)", "longitude");
        
    }else{
        std::cout << "Optimizer error!!!\n";
        std::cout << "error_msg:" << main_solution.error_msg << std::endl;
    }
    

    
    
    
}

