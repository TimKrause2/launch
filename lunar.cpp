#include "psopt.h"
#define GL_GLEXT_PROTOTYPES 1
#include <GL/gl.h>
#include <GL/glext.h>
//#include <GL/glew.h>
#include "SDL.h"
#include <stdio.h>
#include <math.h>
#define GLM_FORCE_RADIANS
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <thread>
#include <iostream>
#include "font.h"
#include "orbit.h"
#include "gsim.h"
#include "trajectory.h"
#include "uv_sphere.h"
#include "event.h"
#include "satellite.h"

using namespace PSOPT;

typedef Eigen::Matrix<adouble, 3, 1> aVector3d;
typedef Eigen::Matrix<adouble, 3, 3> aMatrix3d;

bool mousing=false;
bool up_key=false;
bool down_key=false;
bool left_key=false;
bool right_key=false;
int last_x;
int last_y;
int modifiers;
Font font;

#define N_POINTS 3000
#define FONT_SIZE 16
#define RAD_PER_DEG (M_PI/180.0)

OrbitalParams param0;

double view_lat=45.0;
double view_long=90.0;

double mass_earth = 5.97219e24;
double radius_earth = 6.371e6;
double mass_moon = 7.3477e22;
double radius_moon = 1.7371e6;
double mass_vehicle = 1e3;
double radius_vehicle = 10.0;
double h_vehicle = 500e3;
double T_earth_moon = 27.321662*24*60*60;
double alpha_apogee = 0.975;

int window_width;
int window_height;

#define E_MAX 0.9
#define FOV_Y  25.0

Body* body_earth;
Body* body_moon;
Satellite* body_vehicle;

System sim_system;

Trajectory trajectory_vehicle(N_POINTS);
Trajectory trajectory_moon(N_POINTS);
Trajectory trajectory_earth(N_POINTS);

#define N_NODES_GUESS 10

//
// state vector format
// row(0-5) modified equinoctial orbital elements

// control vector format
// row(0) latitude
// row(1) longitude

Eigen::MatrixXd    guess_control(2,N_NODES_GUESS);
Eigen::MatrixXd    guess_time(1,N_NODES_GUESS);
double             guess_T_burn;

std::thread trajectory_thread;
int i_thread;

bool thread_finished;

enum display_state
{
    DS_INIT,
    DS_INIT_INSERT,
    DS_CALCULATING,
    DS_DISPLAY_ALL,
    DS_ANIMATE,
    DS_ANIMATE_PAUSE
};

int i_animate;

display_state dstate = DS_INIT;

double initialize_system(void){

    // initialize the earth moon system
    orbit_initialize(
                param0,
                mass_earth,
                mass_moon,
                body_earth->m_position,
                body_moon->m_position,
                body_earth->m_velocity,
                body_moon->m_velocity);

    // initialize the vehicle position and velocity
    double r_perigee = radius_earth+h_vehicle;
    Eigen::Vector3d r_vehicle;
    r_vehicle << -r_perigee, 0.0, 0.0;
    body_vehicle->m_position = body_earth->m_position +
            r_vehicle;
    double r_apogee = param0.a*alpha_apogee;
    double e_vehicle = (r_apogee-r_perigee)/(r_apogee+r_perigee);
    double mag_v_vehicle = sqrt(G_gravity*(mass_earth+mass_vehicle)*(1+e_vehicle)/r_perigee);
    Eigen::Vector3d v_vehicle;
    v_vehicle << 0.0, -mag_v_vehicle, 0.0;
    body_vehicle->m_velocity = body_earth->m_velocity +
            v_vehicle;
    double a_vehicle = (r_perigee+r_apogee)/2.0;
    double T_vehicle = period(a_vehicle,mass_earth,mass_vehicle);

    // delta-t to integrate the trajectory
    double dt = T_vehicle*3.0/(N_POINTS-1);
    return dt;
}

Eigen::Vector3d body_acceleration(
    Eigen::Vector3d &r,
    Body *body)
{
    Eigen::Vector3d r12 = r - body->m_position;
    double mag_r12 = r12.norm();
    Eigen::Vector3d n21 = -r12/mag_r12;
    mag_r12 *= mag_r12;
    double acc = body->m_mass*G_gravity/mag_r12;
    return acc*n21;
}

double gravity_factor(void){
    Eigen::Vector3d acc_earth = body_acceleration(
        body_vehicle->m_position, body_earth);
    Eigen::Vector3d acc_moon = body_acceleration(
        body_vehicle->m_position, body_moon);
    double m_acc_earth = acc_earth.norm();
    double m_acc_moon = acc_moon.norm();
    return m_acc_moon/(m_acc_earth+m_acc_moon);
}


void integrate_trajectory(double dt)
{
    for(int i_thread=0;i_thread<N_POINTS;i_thread++){
        trajectory_vehicle.SetPoint(i_thread, body_vehicle);
        trajectory_moon.SetPoint(i_thread, body_moon);
        trajectory_earth.SetPoint(i_thread, body_earth);
        sim_system.rkIntegrate(dt);
    }
    thread_finished=true;
}

//////////////////////////////////////////////////////////////////////////
///////////////////  Define the end point (Mayer) cost function //////////
//////////////////////////////////////////////////////////////////////////


adouble endpoint_cost(adouble* initial_states, adouble* final_states,
                      adouble* parameters,adouble& t0, adouble& tf,
                      adouble* xad, int iphase, Workspace* workspace)
{
    return tf;
}

//////////////////////////////////////////////////////////////////////////
///////////////////  Define the integrand (Lagrange) cost function  //////
//////////////////////////////////////////////////////////////////////////

adouble integrand_cost(adouble* states, adouble* controls, adouble* parameters,
                       adouble& time, adouble* xad, int iphase, Workspace* workspace)
{
    return 0;
}


//////////////////////////////////////////////////////////////////////////
///////////////////  Define the DAE's ////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

aVector3d perturbation_acceleration(
    aVector3d r,
    Eigen::Vector3d pri_body,
    Eigen::Vector3d sec_body,
    double mu)
{
    Eigen::Vector3d s_e = sec_body - pri_body;
    aVector3d s;
    s(0) = s_e(0);
    s(1) = s_e(1);
    s(2) = s_e(2);
    aVector3d d = r - s;
    adouble mag_d = d.norm();
    adouble mag_s = d.norm();
    aVector3d t = -mu*(d/mag_d/mag_d/mag_d + s/mag_s/mag_s/mag_s);
    return t;
}

void dae(adouble* derivatives, adouble* path, adouble* states,
         adouble* controls, adouble* parameters, adouble& time,
         adouble* xad, int iphase, Workspace* workspace)
{
    // Local integers
    int i, j;

    double mu = G_gravity*(mass_moon+mass_vehicle);
    double mu_earth = G_gravity*mass_earth;


    // Extract individual variables

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

    // r and v
    aVector3d rvec;
    rvec(0) = r/s2*( cos(L) + alpha2*cos(L) + 2*h*k*sin(L));
    rvec(1) = r/s2*( sin(L) - alpha2*sin(L) + 2*h*k*cos(L));
    rvec(2) = 2*r/s2*( h*sin(L) - k*cos(L) );

    aVector3d vvec;
    vvec(0) = -(1.0/s2)*sqrt(mu/p)*(  sin(L) + alpha2*sin(L) - 2*h*k*cos(L) + g - 2*f*h*k + alpha2*g);
    vvec(1) = -(1.0/s2)*sqrt(mu/p)*( -cos(L) + alpha2*cos(L) + 2*h*k*sin(L) - f + 2*g*h*k + alpha2*f);
    vvec(2) =  (2.0/s2)*sqrt(mu/p)*(h*cos(L) +      k*sin(L) + f*h + g*k);


    // compute Qr

    aVector3d ir, ith, ih;
    aVector3d rv;
    aVector3d rvr;

    cross( rvec.data(), vvec.data(), rv.data() );
    cross( rv.data(), rvec.data(), rvr.data() );

    ir = rvec.normalized();
    ith = rvr.normalized();
    ih = rv.normalized();

    aMatrix3d Qr;
    Qr.col(0) = ir;
    Qr.col(1) = ith;
    Qr.col(2) = ih;

    adouble latitude  = u[ 0 ];
    adouble longitude = u[ 1 ];

    adouble Tvec_tangent = THRUST*cos(latitude)*cos(longitude);
    adouble Tvec_normal = THRUST*cos(latitude)*sin(longitude);
    adouble Tvec_radial = THRUST*sin(latitude);

    aVector3d a_body = perturbation_acceleration(
        rvec, body_moon->m_position, body_earth->m_position, mu_earth);

    aVector3d delta_body = Qr.transpose()*a_body;

    adouble delta1= Tvec_radial  + delta_body(0);
    adouble delta2= Tvec_tangent + delta_body(1);
    adouble delta3= Tvec_normal  + delta_body(2);


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

void events(adouble* e, adouble* initial_states, adouble* final_states,
            adouble* parameters,adouble& t0, adouble& tf, adouble* xad,
            int iphase, Workspace* workspace)

{
    Eigen::Map<Eigen::Matrix<adouble, 6, 1>> x_i(initial_states, 6);
    Eigen::Map<Eigen::Matrix<adouble, 6, 1>> x_f(final_states, 6);
    Eigen::Map<Eigen::Matrix<adouble, 11, 1>> e_vec(e, 11);

    e_vec << x_i, x_f.head(5);
}



///////////////////////////////////////////////////////////////////////////
///////////////////  Define the phase linkages function ///////////////////
///////////////////////////////////////////////////////////////////////////

void linkages( adouble* linkages, adouble* xad, Workspace* workspace)
{



}

void insert_trajectory(double dt)
{
    // push the state of the system
    sim_system.PushState();

    // mark the beginning of the coast
    TimeSpec ts_coast_start = body_vehicle->GetCurrentTime();
    TimeSpec ts_coast_end;

    // integrate util the acceleration factor starts to deminish
    double acc_factor_last = 0.0;
    double acc_factor;
    bool found = false;
    for(int i=0;i<N_POINTS;i++){
        acc_factor = gravity_factor();
        if(acc_factor<acc_factor_last){
            // found the maximum
            found = true;
            break;
        }
        double dt_s;
        if(acc_factor>0.9){
            dt_s = 60.0;
        }else{
            dt_s = dt;
        }
        sim_system.rkIntegrate(dt_s);
        acc_factor_last = acc_factor;
    }

    if(!found || !(acc_factor > 0.95)){
        sim_system.RestoreState();
        sim_system.PopState();
        integrate_trajectory(dt);
        return;
    }

    // mark the end of the coast
    ts_coast_end = body_vehicle->GetCurrentTime();

    // compute the delta-v to enter a circular orbit
    Eigen::Vector3d r = body_vehicle->m_position - body_moon->m_position;
    Eigen::Vector3d v = body_vehicle->m_velocity - body_moon->m_velocity;
    Eigen::Matrix3d Av = vehicle_basis(r,v);

    // resolve the velocity in the vehicle frame
    Eigen::Vector3d v_vehicle = Av.transpose()*v;
    double mag_v1 = v_vehicle(0);
    double mag_v2 = v_circular(r.norm(), mass_vehicle, mass_moon);
    double delta_v = mag_v1 - mag_v2;

    // compute the burn time from the thrust
    guess_T_burn = delta_v/THRUST;

    // initialize the direction coordinate state and time vector
    double latitude = 0.0;
    double longitude = M_PI;
    guess_control.row(0) = latitude*ones(1,N_NODES_GUESS);
    guess_control.row(1) = longitude*ones(1,N_NODES_GUESS);

    guess_time = linspace(0.0, guess_T_burn, N_NODES_GUESS);

    // find the initial and final modified equinoctial orbital elements
    OrbitalParams oparams0;
    OrbitalParams oparams1;
    body_vehicle->primary_body = body_moon;
    body_vehicle->OrbitalElements(oparams0);
    oparams1 = oparams0;
    oparams1.a = r.norm();
    oparams1.e = 0.0;
    ModEquOrbitalParams meparams0;
    ModEquOrbitalParams meparams1;
    meparams0.FromOrbitalParams(oparams0);
    meparams1.FromOrbitalParams(oparams1);
    Eigen::VectorXd meparams0_vec = meparams0.ToVector();
    Eigen::VectorXd meparams1_vec = meparams1.ToVector();


    ////////////////////////////////////////////////////////////////////////////
    ///////////////////  Declare key structures ////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////

    Alg  algorithm;
    Sol  solution;
    Prob problem;

    ////////////////////////////////////////////////////////////////////////////
    ///////////////////  Register problem name  ////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////

    problem.name        = "Insertion Burn";
    problem.outfilename = "insertion_burn.txt";

    ////////////////////////////////////////////////////////////////////////////
    ////////////  Define problem level constants & do level 1 setup ////////////
    ////////////////////////////////////////////////////////////////////////////

    problem.nphases   			= 1;
    problem.nlinkages                   = 0;

    psopt_level1_setup(problem);

    /////////////////////////////////////////////////////////////////////////////
    /////////   Define phase related information & do level 2 setup  ////////////
    /////////////////////////////////////////////////////////////////////////////

    problem.phases(1).nstates     = 6;
    problem.phases(1).ncontrols   = 2;
    problem.phases(1).nparameters = 0;
    problem.phases(1).nevents     = 11;
    problem.phases(1).npath       = 0;
    problem.phases(1).nodes      << N_NODES_GUESS;

    psopt_level2_setup(problem, algorithm);

    ////////////////////////////////////////////////////////////////////////////
    ///////////////////  Enter problem bounds information //////////////////////
    ////////////////////////////////////////////////////////////////////////////

    problem.phases(1).bounds.lower.controls(0) = -pi/4;
    problem.phases(1).bounds.upper.controls(0) = pi/4;
    problem.phases(1).bounds.lower.controls(1) = -2.0*pi;
    problem.phases(1).bounds.upper.controls(1) = 2.0*pi;

    problem.phases(1).bounds.lower.states(0) = 0.0; // eccentricity might be close to one
    problem.phases(1).bounds.lower.states(1) = -2; // maximum eccentricity
    problem.phases(1).bounds.lower.states(2) = -2; // maximum eccentricity
    problem.phases(1).bounds.lower.states(3) = -1; // inclination less then pi/2
    problem.phases(1).bounds.lower.states(4) = -1;
    problem.phases(1).bounds.lower.states(5) = -pi;


    problem.phases(1).bounds.upper.states(0) = r.norm()*20.0;
    problem.phases(1).bounds.upper.states(1) = 2;
    problem.phases(1).bounds.upper.states(2) = 2;
    problem.phases(1).bounds.upper.states(3) = 1;
    problem.phases(1).bounds.upper.states(4) = 1;
    problem.phases(1).bounds.upper.states(5) = 30*pi;

    Eigen::VectorXd event(11,1);
    event << meparams0_vec, meparams1_vec.head(5);
    problem.phases(1).bounds.lower.events = event;
    problem.phases(1).bounds.upper.events = event;

    problem.phases(1).bounds.lower.StartTime    = 0.0;
    problem.phases(1).bounds.upper.StartTime    = 0.0;

    problem.phases(1).bounds.lower.EndTime      = guess_T_burn*0.1;
    problem.phases(1).bounds.upper.EndTime      = guess_T_burn*1.5;


    ////////////////////////////////////////////////////////////////////////////
    ///////////////////  Register problem functions  ///////////////////////////
    ////////////////////////////////////////////////////////////////////////////


    problem.integrand_cost  = &integrand_cost;
    problem.endpoint_cost   = &endpoint_cost;
    problem.dae             = &dae;
    problem.events          = &events;
    problem.linkages        = &linkages;

    ////////////////////////////////////////////////////////////////////////////
    ///////////////////  Define & register initial guess ///////////////////////
    ////////////////////////////////////////////////////////////////////////////

    int nnodes;
    int nstates;
    int iphase;
    MatrixXd x_guess, param_guess, xini, xfinal;

    nnodes  = problem.phases(1).nodes(0);
    nstates = problem.phases(1).nstates;
    iphase = 1;

    xini = meparams0_vec;

    rk4_propagate( dae, guess_control, guess_time, xini, param_guess, problem, iphase, x_guess, NULL);

    xfinal = x_guess.col(x_guess.cols()-1);

    ModEquOrbitalParams meparams_pro;
    meparams_pro.FromVector(xfinal);
    OrbitalParams oparams_pro;
    meparams_pro.ToOrbitalParams(oparams_pro);





    problem.phases(1).guess.states   = x_guess;
    problem.phases(1).guess.time     = guess_time;
    problem.phases(1).guess.controls = guess_control;

    ////////////////////////////////////////////////////////////////////////////
    ///////////////////  Enter algorithm options  //////////////////////////////
    ////////////////////////////////////////////////////////////////////////////


    algorithm.nlp_iter_max                = 1000;
    algorithm.nlp_tolerance               = 1.e-6;
    algorithm.nlp_method                  = "IPOPT";
    algorithm.scaling                     = "automatic";
    algorithm.derivatives                 = "automatic";
    algorithm.defect_scaling              = "jacobian-based";
    algorithm.collocation_method          = "trapezoidal";
    algorithm.mesh_refinement             = "automatic";
    algorithm.ode_tolerance               = 1.e-5;
    algorithm.mr_max_iterations           = 2;



    ////////////////////////////////////////////////////////////////////////////
    ///////////////////  Now call PSOPT to solve the problem   //////////////////
    ////////////////////////////////////////////////////////////////////////////

    psopt(solution, problem, algorithm);

    ////////////////////////////////////////////////////////////////////////////
    ///////////  Extract relevant variables from solution structure   //////////
    ////////////////////////////////////////////////////////////////////////////

    sim_system.RestoreState();
    sim_system.PopState();

    if(!solution.error_flag){
        Eigen::MatrixXd opt_state = solution.get_states_in_phase(1);
        Eigen::MatrixXd opt_control = solution.get_controls_in_phase(1);
        Eigen::RowVectorXd opt_time = solution.get_time_in_phase(1);

        Eigen::VectorXd opt_final = opt_state.col(opt_state.cols()-1);
        ModEquOrbitalParams meparams_opt;
        meparams_opt.FromVector(opt_final);
        OrbitalParams oparams_opt;
        meparams_opt.ToOrbitalParams(oparams_opt);

        Save(opt_time, "opt_time.dat");

        Eigen::MatrixXd opt_latitude = opt_control.row(0)/RAD_PER_DEG;
        Eigen::MatrixXd opt_longitude = opt_control.row(1)/RAD_PER_DEG;

        plot(opt_time, opt_latitude, problem.name+": latitude","time (s)","latitude (deg)", "latitude");
        plot(opt_time, opt_longitude, problem.name+": longitude", "time (s)", "longitude (deg)", "longitude");

        double T_coast = (ts_coast_end - ts_coast_start).to_double();

        body_vehicle->ScheduleDirectedAccelerationManeuver(
            opt_control, opt_time, THRUST, T_coast);
    }

    integrate_trajectory(dt);

}

void display(void)
{
    // calculate the radius of the scene
    double r_scene = param0.a*mass_earth/(mass_earth+mass_moon)*(1+param0.e);

    // calculate the distance of the camera
    double r_camera = r_scene/sin(FOV_Y*RAD_PER_DEG/2.0);

    // calculate the near and far planes
    float near = r_camera - r_scene;
    float far = r_camera + r_scene;

    // calculate the orientation of the camera
    glm::dmat4 A_lat = glm::rotate(
                glm::dmat4(1.0),
                (90.0-view_lat)*RAD_PER_DEG,
                glm::dvec3(1.0,0.0,0.0));
    glm::dmat4 A_long = glm::rotate(
                glm::dmat4(1.0),
                view_long*RAD_PER_DEG,
                glm::dvec3(0.0,0.0,1.0));

    glm::dmat4 A_cam_rot = A_long*A_lat;

    glm::dvec4 T_cam_v4 = A_cam_rot*glm::dvec4(0.0,0.0,r_camera,0.0);
    glm::dvec3 T_cam(T_cam_v4);

    glm::dmat4 A_cam_trans = glm::translate(glm::dmat4(1.0),T_cam);

    // calculate the camera space to global space transform
    glm::dmat4 A_cam = A_cam_trans*A_cam_rot;

    // initialize the perspective transfrom
    glm::mat4 A_pers = glm::perspective(
        (float)(FOV_Y*RAD_PER_DEG),
        (float)window_width/(float)window_height,
        near,
        far);

    glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);

    glMatrixMode(GL_PROJECTION);
    glLoadMatrixf( glm::value_ptr(A_pers) );

    glMatrixMode(GL_MODELVIEW);
    glm::dmat4 A_cam_inv = glm::inverse(A_cam);
    glLoadMatrixd( glm::value_ptr(A_cam_inv) );
    glm::mat4 view(A_cam_inv);
    orbit_draw(param0,mass_earth,mass_moon);

    double dt;
    switch(dstate){
    case DS_INIT:
        dt = initialize_system();
        thread_finished = false;
        trajectory_thread = std::thread(integrate_trajectory,dt);
        dstate = DS_CALCULATING;
        break;

    case DS_INIT_INSERT:
        dt = initialize_system();
        thread_finished = false;
        trajectory_thread = std::thread(insert_trajectory,dt);
        dstate = DS_CALCULATING;
        break;

    case DS_CALCULATING:
        // draw the progress of the trajectory
        glColor3f(0.5f,0.5f,1.0f);
        trajectory_vehicle.Draw(i_thread);
        glColor3f(1.0f,1.0f,0.5f);
        trajectory_moon.Draw(i_thread);
        glColor3f(0.5f,1.0f,1.0f);
        trajectory_earth.Draw(i_thread);
        if(thread_finished){
            trajectory_thread.join();
            dstate = DS_DISPLAY_ALL;
        }
        break;
    case DS_DISPLAY_ALL:
        glColor3f(0.5f,0.5f,1.0f);
        trajectory_vehicle.Draw(N_POINTS);
        glColor3f(1.0f,1.0f,0.5f);
        trajectory_moon.Draw(N_POINTS);
        glColor3f(0.5f,1.0f,1.0f);
        trajectory_earth.Draw(N_POINTS);
        sim_system.draw(A_pers, view);
        break;

    case DS_ANIMATE:
        glColor3f(0.5f,0.5f,1.0f);
        trajectory_vehicle.Draw(i_animate);
        trajectory_vehicle.GetPoint(i_animate, body_vehicle);
        glColor3f(1.0f,1.0f,0.5f);
        trajectory_moon.Draw(i_animate);
        trajectory_moon.GetPoint(i_animate, body_moon);
        glColor3f(0.5f,1.0f,1.0f);
        trajectory_earth.Draw(i_animate);
        trajectory_earth.GetPoint(i_animate, body_earth);
        sim_system.draw(A_pers, view);
        i_animate++;
        i_animate%=N_POINTS;
        break;
    case DS_ANIMATE_PAUSE:
        glColor3f(0.5f,0.5f,1.0f);
        trajectory_vehicle.Draw(i_animate);
        trajectory_vehicle.GetPoint(i_animate, body_vehicle);
        glColor3f(1.0f,1.0f,0.5f);
        trajectory_moon.Draw(i_animate);
        trajectory_moon.GetPoint(i_animate, body_moon);
        glColor3f(0.5f,1.0f,1.0f);
        trajectory_earth.Draw(i_animate);
        trajectory_earth.GetPoint(i_animate, body_earth);
        sim_system.draw(A_pers, view);
        break;
    }

    glColor3f(1.0f,1.0f,1.0f);
    // print some stuff
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho( 0.0, (GLdouble)window_width,
            (GLdouble)window_height, 0.0, -1.0, 1.0 );
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    double x=FONT_SIZE;
    double y=FONT_SIZE;

    font.Printf(x,y,"alpha_apogee:%.4lf",alpha_apogee);
    y+=FONT_SIZE;
    font.Printf(x,y,"nu moon:%.4lf",param0.nu/RAD_PER_DEG);
    if(dstate == DS_ANIMATE || dstate == DS_ANIMATE_PAUSE){
        double factor = gravity_factor();
        y+=FONT_SIZE;
        font.Printf(x,y,"acceleration factor:%5.2lf", factor*100.0);
        Eigen::Vector3d v12 = body_vehicle->m_velocity -
                            body_moon->m_velocity;
        double m_v12 = v12.norm();
        Eigen::Vector3d r12 = body_vehicle->m_position -
                              body_moon->m_position;
        y+=FONT_SIZE;
        font.Printf(x,y,"relative velocity:%.2lf", m_v12);
        y+=FONT_SIZE;
        font.Printf(x,y,"moon altitude:%.2lf", r12.norm() - radius_moon);
    }

}

void init_system(void)
{

    param0.a = semi_major(T_earth_moon,mass_earth,mass_moon);
    param0.e = 0.0;
    param0.i = 0.0;
    param0.Omega = 0.0;
    param0.omega = 0.0;
    param0.nu = -65.75*RAD_PER_DEG;

    Eigen::Vector3d pos_earth;
    Eigen::Vector3d pos_moon;
    Eigen::Vector3d pos_vehicle;
    Eigen::Vector3d vel_earth;
    Eigen::Vector3d vel_moon;
    Eigen::Vector3d vel_vehicle;

    Eigen::Vector4d p_identity(1.0, 0.0, 0.0, 0.0);

    Eigen::Vector4d p_earth = p_identity;
    Eigen::Vector4d p_moon = p_identity;
    Eigen::Vector4d p_vehicle = p_identity;
    Eigen::Vector4d pdot_earth = Eigen::Vector4d::Zero();
    Eigen::Vector4d pdot_moon = Eigen::Vector4d::Zero();
    Eigen::Vector4d pdot_vehicle = Eigen::Vector4d::Zero();
    Eigen::Matrix3d Jp_earth = caams::J_p_sphere(mass_earth,radius_earth);
    Eigen::Matrix3d Jp_moon = caams::J_p_sphere(mass_moon,radius_moon);
    Eigen::Matrix3d Jp_vehicle = caams::J_p_cylinder_z_axis(mass_vehicle,radius_vehicle*0.1,radius_vehicle);
    Eigen::Vector3d rp_vehicle(0.0, 0.0, -radius_vehicle);

    body_earth = new Body(
                mass_earth,
                radius_earth,
                pos_earth,
                vel_earth,
                p_earth,
                pdot_earth,
                Jp_earth,
                (char*)"earth_8k.jpg");

    body_moon = new Body(
                mass_moon,
                radius_moon,
                pos_moon,
                vel_moon,
                p_moon,
                pdot_moon,
                Jp_moon,
                (char*)"moon_8k.tif");

    body_vehicle = new Satellite(
                mass_vehicle,
                radius_vehicle,
                pos_vehicle,
                vel_vehicle,
                p_vehicle,
                pdot_vehicle,
                Jp_vehicle,
                (char*)"sat_tex.jpg",
                body_earth,
                rp_vehicle);

    sim_system.AddBody(body_earth);
    sim_system.AddBody(body_moon);
    sim_system.AddBody(body_vehicle);

    glClearColor(0.0,0.0,0.0,0.0);
    glDrawBuffer(GL_BACK);

    font.LoadOutline(
                "/usr/share/fonts/truetype/ubuntu/UbuntuMono-R.ttf",
                FONT_SIZE,Pixel32(0, 0, 0),Pixel32(255, 255, 255), 1.0 );

    uv_sphere_init(64,128);


    OrbitalParams c0;

}


void keydown(SDL_Event *event)
{
    switch(event->key.keysym.sym){
    case SDLK_RETURN:
        if(dstate == DS_ANIMATE){
            dstate = DS_ANIMATE_PAUSE;
        }else if(dstate == DS_ANIMATE_PAUSE){
            dstate = DS_ANIMATE;
        }
        break;
    case SDLK_UP:
        if(dstate == DS_DISPLAY_ALL){
            alpha_apogee+=0.001;
            dstate = DS_INIT;
        }
        break;
    case SDLK_DOWN:
        if(dstate == DS_DISPLAY_ALL){
            alpha_apogee-=0.001;
            dstate = DS_INIT;
        }
        break;
    case SDLK_LEFT:
        if(dstate == DS_DISPLAY_ALL){
            param0.nu -= 0.1*RAD_PER_DEG;
            dstate = DS_INIT;
        }else if(dstate == DS_ANIMATE_PAUSE){
            if(i_animate){
                i_animate--;
            }else{
                i_animate = N_POINTS-1;
            }
        }
        break;
    case SDLK_RIGHT:
        if(dstate == DS_DISPLAY_ALL){
            param0.nu += 0.1*RAD_PER_DEG;
            dstate = DS_INIT;
        }else if(dstate == DS_ANIMATE_PAUSE){
            i_animate++;
            i_animate%=N_POINTS;
        }
        break;
    case SDLK_SPACE:
        if(dstate==DS_DISPLAY_ALL){
            i_animate = 0;
            dstate = DS_ANIMATE;
        }else if(dstate==DS_ANIMATE){
            dstate = DS_DISPLAY_ALL;
        }
        break;
    case SDLK_F1:
        if(dstate==DS_DISPLAY_ALL){
            dstate=DS_INIT_INSERT;
        }
        break;

    default:
        break;
    }
}

void keyup(SDL_Event *event)
{
    switch(event->key.keysym.sym){
    default:
        break;
    }
}

void mousemotion(SDL_Event *event)
{
    int x = event->motion.x;
    int y = event->motion.y;
    if( mousing ){
        int delta_x_int = x - last_x;
        int delta_y_int = y - last_y;
    }
    last_x = x;
    last_y = y;
}

void mousebuttondown(SDL_Event *event)
{
    if(event->button.button == SDL_BUTTON_LEFT){
        mousing = true;
        last_x = event->button.x;
        last_y = event->button.y;
    }
}

void mousebuttonup(SDL_Event *event)
{
    if(event->button.button == SDL_BUTTON_LEFT){
        mousing = false;
    }
}

int main(int argc, char **argv)
{
    if(SDL_Init(SDL_INIT_VIDEO)!=0){
        SDL_Log("Unable to initialize SDL: %s", SDL_GetError());
        return 1;
    }

    SDL_Window *window;
    window = SDL_CreateWindow(
        "Lunar Approach",
        SDL_WINDOWPOS_UNDEFINED,
        SDL_WINDOWPOS_UNDEFINED,
        1024,
        720,
        SDL_WINDOW_OPENGL|SDL_WINDOW_RESIZABLE);

    if(window==NULL){
        printf("Couldn't create window: %s\n",SDL_GetError());
        SDL_Quit();
        return 1;
    }

    SDL_GL_SetAttribute( SDL_GL_RED_SIZE, 8 );
    SDL_GL_SetAttribute( SDL_GL_GREEN_SIZE, 8 );
    SDL_GL_SetAttribute( SDL_GL_BLUE_SIZE, 8 );
    SDL_GL_SetAttribute( SDL_GL_ALPHA_SIZE, 8 );
    SDL_GL_SetAttribute( SDL_GL_DOUBLEBUFFER, 1 );
    SDL_GL_SetAttribute( SDL_GL_DEPTH_SIZE, 24 );
    SDL_GL_SetAttribute( SDL_GL_CONTEXT_PROFILE_MASK, SDL_GL_CONTEXT_PROFILE_COMPATIBILITY );
    SDL_GL_SetAttribute( SDL_GL_CONTEXT_MAJOR_VERSION, 3 );
    SDL_GL_SetAttribute( SDL_GL_CONTEXT_MINOR_VERSION, 0 );

    SDL_GLContext glcontext = SDL_GL_CreateContext(window);

    if(glcontext==NULL){
        printf("Couldn't create GL context: %s\n", SDL_GetError());
        SDL_DestroyWindow(window);
        SDL_Quit();
        return 1;
    }

    // GLenum err = glewInit();
    // if(err!=GLEW_OK){
    //     printf("glewInit error:%s\n",glewGetErrorString(err));
    //     return -1;
    // }

    SDL_ShowWindow(window);

    init_system();

    SDL_GetWindowSize(window, &window_width, &window_height);

    int running = 1;
    SDL_Event event;
    while(running){
        while(SDL_PollEvent(&event)){
            switch(event.type){
            case SDL_KEYDOWN:
                if(event.key.keysym.scancode == SDL_SCANCODE_ESCAPE)
                    running=0;
                else
                    keydown(&event);
                break;
            case SDL_KEYUP:
                keyup(&event);
                break;
            case SDL_MOUSEMOTION:
                mousemotion(&event);
                break;
            case SDL_MOUSEBUTTONDOWN:
                mousebuttondown(&event);
                break;
            case SDL_MOUSEBUTTONUP:
                mousebuttonup(&event);
                break;
            case SDL_QUIT:
                running=0;
                break;
            case SDL_WINDOWEVENT:
                switch(event.window.event){
                case SDL_WINDOWEVENT_RESIZED:
                case SDL_WINDOWEVENT_SIZE_CHANGED:
                    window_width = event.window.data1;
                    window_height = event.window.data2;
                    break;
                default:
                    break;
                }
            default:
                break;
            }
        }

        if(!running)
            break;

        glViewport(0, 0, window_width, window_height);

        display();

        SDL_GL_SwapWindow(window);
    }
    SDL_GL_DeleteContext(glcontext);
    SDL_DestroyWindow(window);
    SDL_Quit();
    return 0;
}


