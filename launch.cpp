//#include <epoxy/gl.h>
//#include <GL/gl.h>
#include <GLES3/gl32.h>
#include <GLES3/gl3ext.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#define GLM_FORCE_RADIANS
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtx/string_cast.hpp>
#include "launch.h"
#include "launchapp.h"

LaunchData ld;


bool mousing=false;
int last_x;
int last_y;
int modifiers;

#define RAD_PER_DEG (M_PI/180.0)

double view_lat=0.0;
double view_long=-90.0;
double mass_earth = 5.97219e24;
double radius_earth = 6.371e6;
double g_earth = 9.80665;
double a_WGS84 = 6378137;
double one_over_f = 298.257223563;
double mag_omega_earth = (2*M_PI/86164.0905);
double mass_sat = 1e3;
double radius_sat = 5.0;
double launch_latitude = 52.27024;
double launch_longitude = -113.813697;
double launch_heading = 90.0;

int window_width;
int window_height;
double time_factor=1.0;

#define FOV_Y         60.0
#define CAM_NEAR      20.0f
#define CAM_FAR       (6.371e6f*2.0f + 210e6f)
#define CAM_DISTANCE  46.0
#define CAM_DISTANCE_MIN 25.0
#define CAM_DISTANCE_MAX 1E6
#define CAM_DISTANCE_FACTOR 1.1
#define THRUST        (9.819609*2.0)
#define DELTA_T       (1.0/60.0)
#define FONT_SIZE     24
#define BL_LAT_VEHICLE 48.935333180089664*RAD_PER_DEG
#define BL_LONG_VEHICLE -54.58524740913247*RAD_PER_DEG
#define BL_LAT_TARGET 35.94248310361987*RAD_PER_DEG
#define BL_LONG_TARGET -5.617076200801353*RAD_PER_DEG
#define BL_T1 500
#define BL_T2 1800

float cam_distance = CAM_DISTANCE;


Spheroid* body_earth;
Satellite* body_sat;
System earth_sat_system;

Eigen::Vector3d n_h_launch;

OrbitalParams params_new;
OrbitalParams params_dl;

LaunchData::LaunchData()
{
    view_lat = 0.0;
    view_long = -90.0;
    time_factor = 1.0;
    cam_distance = CAM_DISTANCE;
}

void random_orbit(void)
{
    double r_p = radius_earth + 499980.0;
    double r_a = radius_earth + 500020.0 + 10000000.0*drand48();
    params_new.a = (r_a+r_p)/2.0;
    params_new.e = (r_a-r_p)/(r_a+r_p);
    params_new.i = (0.001 + 89.0*drand48())*RAD_PER_DEG;
    params_new.Omega = (-180.0 + 360.0*drand48())*RAD_PER_DEG;
    params_new.omega = (-180.0 + 360.0*drand48())*RAD_PER_DEG;
    params_new.nu =    1.0*RAD_PER_DEG;//(-180.0 + 360.0*drand48())*RAD_PER_DEG;
}

void LaunchData::display(void)
{
    earth_sat_system.origin( body_sat->m_position );

    Eigen::Vector3d r_sat_earth = body_sat->m_position - body_earth->m_position;
    Eigen::Vector3d n_r = r_sat_earth.normalized();
    Eigen::Vector3d v_sat_earth = body_sat->m_velocity - body_earth->m_velocity;
    Eigen::Vector3d h = r_sat_earth.cross(v_sat_earth);
    Eigen::Vector3d n_h;
    double mag_h = h.norm();
    if(mag_h==0.0){
        n_h = n_h_launch;
    }else{
        n_h = h.normalized();
    }
    Eigen::Vector3d n_theta = n_h.cross(n_r);
    Eigen::Matrix3d cam_rot0;
    cam_rot0.col(0) = n_h;
    cam_rot0.col(1) = n_r;
    cam_rot0.col(2) = n_theta;

    Eigen::Matrix3d cam_rot_lat = caams::AAA(view_lat*RAD_PER_DEG,-n_h);
    Eigen::Matrix3d cam_rot_long = caams::AAA(view_long*RAD_PER_DEG,n_r);

    Eigen::Matrix3d cam_rot = cam_rot_long*cam_rot_lat*cam_rot0;
    Eigen::Vector3d cam_tran0;
    cam_tran0 << 0.0, 0.0, cam_distance;
    Eigen::Vector3d cam_tran = cam_rot*cam_tran0 + body_sat->m_position;

    Eigen::Matrix3d sat_rot;
    sat_rot.col(0) = n_theta;
    sat_rot.col(1) = n_h;
    sat_rot.col(2) = n_r;

    body_sat->p = caams::pA(sat_rot);

    glm::dmat3 glm_cam_rot3 = E2GLM(cam_rot);
    glm::dmat4 glm_cam_rot4(glm_cam_rot3);
    glm::dvec3 glm_cam_tranv = E2GLM(cam_tran);
    glm::dmat4 glm_cam_tranm = glm::translate(glm::dmat4(1.0),glm_cam_tranv);

    glm::dmat4 A_cam = glm_cam_tranm*glm_cam_rot4;


    float dist = (float)r_sat_earth.norm() + radius_earth + cam_distance;

    glm::mat4 A_pers = glm::perspective(
        (float)(FOV_Y*RAD_PER_DEG),
        (float)window_width/(float)window_height,
        CAM_NEAR,
        dist);

    // std::cout << "A_pers:" << std::endl;
    // std::cout << glm::to_string(A_pers) << std::endl;


    glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);


    glm::dmat4 A_cam_inv = glm::inverse(A_cam);

    glm::mat4 view(A_cam_inv);

    earth_sat_system.draw(A_pers,view);

    // print some stuff

    Eigen::Vector3d r = body_sat->m_position - body_earth->m_position;
    Eigen::Vector3d v = body_sat->m_velocity - body_earth->m_velocity;
    OrbitalParams params;

    orbital_elements(
                r,v,body_earth->m_mass,body_sat->m_mass,
                params);

    double x = FONT_SIZE;
    double y = window_height - FONT_SIZE;

    font.Printf(x, y, "Eccentricity:%.15lg", params.e );
    y-=FONT_SIZE;
    font.Printf(x, y, "a:%.15lg h_peri:%12.1lf h_apo:%12.1lf",
                params.a,
                params.a*(1-params.e) - body_earth->m_radius,
                params.a*(1+params.e) - body_earth->m_radius);
    y-=FONT_SIZE;
    font.Printf(x, y, "inclination:%.15lg", params.i*180/M_PI);
    y-=FONT_SIZE;
    font.Printf(x, y, "Omega:%.15lg", params.Omega*180/M_PI);
    y-=FONT_SIZE;
    font.Printf(x, y, "omega:%.15lg", params.omega*180/M_PI);
    y-=FONT_SIZE;
    font.Printf(x, y, "nu:%.7lg", params.nu*180/M_PI);
    y-=FONT_SIZE;
    font.Printf(x, y, "altitude:%.7lg", r.norm()-body_earth->m_radius);
    y-=FONT_SIZE;
    font.Printf(x, y, "time_factor:%4lf",time_factor);
    y-=FONT_SIZE;
    if(body_sat->eventQueue.Empty()){
        font.Printf(x, y, "Next event : Empty");
    }else{
        font.Printf(x, y, "Next event : %lf", body_sat->eventQueue.NextEvent());
    }

    earth_sat_system.rkIntegrate( DELTA_T*time_factor );

}

void LaunchData::init(void)
{
    Eigen::Vector3d pos_earth = Eigen::Vector3d::Zero();
    Eigen::Vector3d vel_earth = Eigen::Vector3d::Zero();

    Eigen::Vector3d pos_sat = Eigen::Vector3d::Zero();
    Eigen::Vector3d vel_sat = Eigen::Vector3d::Zero();

    Eigen::Vector4d p_identity;
    p_identity << 1.0, 0.0, 0.0, 0.0;
    Eigen::Vector4d p_earth = p_identity;
    Eigen::Vector4d p_sat = p_identity;
    Eigen::Vector3d omega_p_earth(0.0, 0.0, mag_omega_earth);
    Eigen::Vector4d pdot_earth = caams::p_dot_omega_p(p_earth, omega_p_earth);
    Eigen::Vector4d pdot_sat = Eigen::Vector4d::Zero();
    Eigen::Matrix3d Jp_earth = caams::J_p_sphere(mass_earth,radius_earth);
    Eigen::Matrix3d Jp_sat = caams::J_p_cylinder_z_axis(mass_sat,radius_sat*0.1,radius_sat);
    Eigen::Vector3d rp_sat(0.0, 0.0, -radius_sat);

    body_earth = new Spheroid(
                mass_earth,
                radius_earth,
                pos_earth,
                vel_earth,
                p_earth,
                pdot_earth,
                Jp_earth,
                (char*)"earthmap_large.jpg",
                mass_earth*1e-5,
                radius_earth*0.9);

    body_earth->SetMass(
                g_earth, a_WGS84, one_over_f);

    body_sat = new Satellite(
                mass_sat,
                radius_sat,
                pos_sat,
                vel_sat,
                p_sat,
                pdot_sat,
                Jp_sat,
                (char*)"sat_tex.jpg",
                body_earth,
                rp_sat);

    // n_h_launch = body_sat->ScheduleLaunchManeuver(
    //             launch_latitude*RAD_PER_DEG,
    //             launch_longitude*RAD_PER_DEG,
    //             launch_heading*RAD_PER_DEG);
    body_sat->ScheduleICBMLaunch(
                30*RAD_PER_DEG, -90*RAD_PER_DEG,
                -10*RAD_PER_DEG, -90*RAD_PER_DEG,
                400, 1000, 20);

    earth_sat_system.AddBody(body_earth);
    earth_sat_system.AddBody(body_sat);

    glm::vec4 fontcolor(0.0,0.0,0.0,1.0);
    glm::vec4 outlinecolor(1.0,1.0,1.0,1.0);


    font.LoadOutline(
                "/usr/share/fonts/truetype/ubuntu/UbuntuMono-R.ttf",
                FONT_SIZE,fontcolor,outlinecolor, 1.0 );

    uv_sphere_init(640,1280);
}

void LaunchData::update_view_direction(int dx, int dy)
{
    view_long -= (double)dx*0.25;
    view_lat += (double)dy*0.25;
    if(view_long>180.0)view_long-=360.0;
    if(view_long<-180.0)view_long+=360.0;
    if(view_lat>90.0)view_lat=90.0;
    if(view_lat<0.0)view_lat=0.0;
}

void LaunchData::update_view_distance(int dr)
{
    cam_distance *= pow(CAM_DISTANCE_FACTOR, dr);
    if(cam_distance < CAM_DISTANCE_MIN){
        cam_distance = CAM_DISTANCE_MIN;
    }else if(cam_distance > CAM_DISTANCE_MAX){
        cam_distance = CAM_DISTANCE_MAX;
    }
}

void LaunchData::update_time_factor(int dt)
{
    time_factor *= pow(10.0f, dt);
    if(time_factor<0.001)time_factor = 0.001;
    if(time_factor>1000.0)time_factor = 1000.0;
}

void LaunchData::update_viewport(int width, int height)
{
    window_width = width;
    window_height = height;
}

// void keydown(SDL_Event *event)
// {
//     switch(event->key.keysym.sym){
//     case SDLK_UP:
//         if(time_factor<1000.0)
//             time_factor *= 10.0;
//         break;
//     case SDLK_DOWN:
//         if(time_factor>1.0)
//             time_factor /= 10.0;
//         break;
//     case SDLK_F1:
//         body_sat->ScheduleCircularizeAtPeriapsis();
//         break;
//     case SDLK_F2:
//         body_sat->ScheduleCircularizeAtApoapsis();
//         break;
//     case SDLK_F3:
//         random_orbit();
//         body_sat->NewOrbit(params_new);
//         break;
//     case SDLK_F4:
//         body_sat->ScheduleAdjustAtPariapsis(radius_earth + 20e6);
//         break;
//     case SDLK_F5:
//         body_sat->ScheduleCircularizeAtRadius(radius_earth + 500e3);
//         break;
//     case SDLK_F6:
//         body_sat->ScheduleAdjustAtApoapsis(radius_earth + 500e3);
//         break;
//     case SDLK_F7:
//         body_sat->ScheduleInclinationManeuver(90.0*RAD_PER_DEG, 2.5*RAD_PER_DEG);
//         break;
//     case SDLK_F8:
//         params_dl.a = body_earth->m_radius + 500e3;
//         params_dl.e = 0.0001;
//         params_dl.i = 51.6 * RAD_PER_DEG;
//         params_dl.Omega = 285.0 * RAD_PER_DEG;
//         params_dl.omega = 90.0 * RAD_PER_DEG;
//         params_dl.nu = 0.0 * RAD_PER_DEG;
//         body_sat->ScheduleDelayedLaunch(
//                     28.608295715431893 * RAD_PER_DEG,
//                     -80.60414807270492 * RAD_PER_DEG,
//                     params_dl, THRUST);
//         break;
//     case SDLK_F9:
//         body_sat->ScheduleBallisticLaunch(
//                     BL_LAT_VEHICLE, BL_LONG_VEHICLE,
//                     BL_LAT_TARGET, BL_LONG_TARGET,
//                     BL_T1, BL_T2,
//                     THRUST);
//     default:
//         break;
//     }
// }

// void keyup(SDL_Event *event)
// {
//     switch(event->key.keysym.sym){
//     default:
//         break;
//     }
// }

// void mousemotion(SDL_Event *event)
// {
//     int x = event->motion.x;
//     int y = event->motion.y;
//     if( mousing ){
//         int delta_x_int = x - last_x;
//         int delta_y_int = y - last_y;
//         view_long -= (double)delta_x_int*0.25;
//         view_lat += (double)delta_y_int*0.25;
//         if(view_long>180.0)view_long-=360.0;
//         if(view_long<-180.0)view_long+=360.0;
//         if(view_lat>90.0)view_lat=90.0;
//         if(view_lat<0.0)view_lat=0.0;
//     }
//     last_x = x;
//     last_y = y;
// }

// void mousebuttondown(SDL_Event *event)
// {
//     if(event->button.button == SDL_BUTTON_LEFT){
//         mousing = true;
//         last_x = event->button.x;
//         last_y = event->button.y;
//     }
// }

// void mousebuttonup(SDL_Event *event)
// {
//     if(event->button.button == SDL_BUTTON_LEFT){
//         mousing = false;
//     }
// }

// void mousewheel(SDL_Event *event)
// {
//     cam_distance *= pow(CAM_DISTANCE_FACTOR, event->wheel.y);
//     if(cam_distance < CAM_DISTANCE_MIN){
//         cam_distance = CAM_DISTANCE_MIN;
//     }else if(cam_distance > CAM_DISTANCE_MAX){
//         cam_distance = CAM_DISTANCE_MAX;
//     }
// }

int main(int argc, char **argv)
{
    // if(SDL_Init(SDL_INIT_VIDEO)!=0){
    //     SDL_Log("Unable to initialize SDL: %s", SDL_GetError());
    //     return 1;
    // }

    // printf("Current video driver:%s\n",SDL_GetCurrentVideoDriver());

    // SDL_Window *window;
    // window = SDL_CreateWindow(
    //     "Maneuver Tester",
    //     SDL_WINDOWPOS_UNDEFINED,
    //     SDL_WINDOWPOS_UNDEFINED,
    //     800,
    //     600,
    //     SDL_WINDOW_OPENGL|SDL_WINDOW_RESIZABLE);

    // if(window==NULL){
    //     printf("Couldn't create window: %s\n",SDL_GetError());
    //     SDL_Quit();
    //     return 1;
    // }

    // SDL_GL_SetAttribute( SDL_GL_RED_SIZE, 8 );
    // SDL_GL_SetAttribute( SDL_GL_GREEN_SIZE, 8 );
    // SDL_GL_SetAttribute( SDL_GL_BLUE_SIZE, 8 );
    // SDL_GL_SetAttribute( SDL_GL_ALPHA_SIZE, 8 );
    // SDL_GL_SetAttribute( SDL_GL_DOUBLEBUFFER, 1 );
    // SDL_GL_SetAttribute( SDL_GL_DEPTH_SIZE, 24 );
    // SDL_GL_SetAttribute( SDL_GL_CONTEXT_PROFILE_MASK, SDL_GL_CONTEXT_PROFILE_COMPATIBILITY );
    // SDL_GL_SetAttribute( SDL_GL_CONTEXT_MAJOR_VERSION, 3 );
    // SDL_GL_SetAttribute( SDL_GL_CONTEXT_MINOR_VERSION, 0 );

    // SDL_GLContext glcontext = SDL_GL_CreateContext(window);

    // if(glcontext==NULL){
    //     printf("Couldn't create GL context: %s\n", SDL_GetError());
    //     SDL_DestroyWindow(window);
    //     SDL_Quit();
    //     return 1;
    // }

    // // GLenum err = glewInit();
    // // if(err!=GLEW_OK){
    // //     printf("glewInit error:%s\n",glewGetErrorString(err));
    // //     return -1;
    // // }

    // SDL_ShowWindow(window);

    // init_program();

    // SDL_GetWindowSize(window, &window_width, &window_height);

    // int running = 1;
    // SDL_Event event;
    // while(running){
    //     while(SDL_PollEvent(&event)){
    //         switch(event.type){
    //         case SDL_KEYDOWN:
    //             if(event.key.keysym.scancode == SDL_SCANCODE_ESCAPE)
    //                 running=0;
    //             else
    //                 keydown(&event);
    //             break;
    //         case SDL_KEYUP:
    //             keyup(&event);
    //             break;
    //         case SDL_MOUSEMOTION:
    //             mousemotion(&event);
    //             break;
    //         case SDL_MOUSEBUTTONDOWN:
    //             mousebuttondown(&event);
    //             break;
    //         case SDL_MOUSEBUTTONUP:
    //             mousebuttonup(&event);
    //             break;
    //         case SDL_MOUSEWHEEL:
    //             mousewheel(&event);
    //             break;
    //         case SDL_QUIT:
    //             running=0;
    //             break;
    //         case SDL_WINDOWEVENT:
    //             switch(event.window.event){
    //             case SDL_WINDOWEVENT_RESIZED:
    //             case SDL_WINDOWEVENT_SIZE_CHANGED:
    //                 window_width = event.window.data1;
    //                 window_height = event.window.data2;
    //                 break;
    //             default:
    //                 break;
    //             }
    //         default:
    //             break;
    //         }
    //     }

    //     if(!running)
    //         break;

    //     glViewport(0, 0, window_width, window_height);

    //     display();

    //     SDL_GL_SwapWindow(window);
    // }
    // SDL_GL_DeleteContext(glcontext);
    // SDL_DestroyWindow(window);
    // SDL_Quit();

    auto application = LaunchApplication::create(ld);

    return application->run(argc, argv);
}


