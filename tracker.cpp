#define GL_GLEXT_PROTOTYPES 1
#include <GL/gl.h>
#include <GL/glext.h>
//#include <GL/glew.h>
#include "SDL.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#define GLM_FORCE_RADIANS
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtx/string_cast.hpp>
#include "font.h"
#include "orbit.h"
#include "gsim.h"
#include "uv_sphere.h"
#include "tracker.h"

bool mousing=false;
int last_x;
int last_y;
int modifiers;
Font font;

#define RAD_PER_DEG (M_PI/180.0)

#define TORQUE_PER_RAD  5000.0
#define TORQUE_PER_RAD_PER_S 5000.0
#define THRUST (9.819609*2.0)
#define THETA_DEFLECTION_MAX (M_PI/4)
#define THETA_TRANSITION (4.0*RAD_PER_DEG)
#define THETA_MAX  (M_PI)



double track_lat=0.0;
double track_long=0.0;
double mass_tracker = 1e3;
double radius_tracker = 5.0;

int window_width;
int window_height;

#define FOV_Y         60.0
#define CAM_NEAR      20.0f
#define CAM_DISTANCE  46.0
#define CAM_DISTANCE_MIN 25.0
#define CAM_DISTANCE_MAX 1E6
#define CAM_DISTANCE_FACTOR 1.1
#define THRUST        (9.819609*2.0)
#define DELTA_T       (1.0/60.0)
#define FONT_SIZE     24

float cam_distance = CAM_DISTANCE;

Tracker *tracker;
System tracker_system;

Tracker::Tracker(
        double mass,
        double radius,
        Eigen::Vector3d &position,
        Eigen::Vector3d &velocity,
        Eigen::Vector4d &p,
        Eigen::Vector4d &pdot,
        Eigen::Matrix3d &Jp,
        char* tex_image) :
    Body(mass, radius, position, velocity,
         p, pdot, Jp, tex_image)
{
    vTrack << 0, 0, 1;
}

void Tracker::Prepare(void)
{
    Eigen::Vector3d z_axis(0,0,1);
    Eigen::Matrix3d B = caams::Ap(p);
    Eigen::Vector3d vTrack_b = B.transpose()*vTrack;
    Eigen::Vector3d n_torque_b = z_axis.cross(vTrack_b);
    Eigen::Vector3d torque_track;
    Eigen::Vector3d torque_friction;
    double mag_n = n_torque_b.norm();
    double alpha_friction = 1.0;
    if(mag_n==0){
        torque_track = Eigen::Vector3d::Zero();
    }else{
        n_torque_b.normalize();
        double vdot = vTrack_b(2);
        if(vdot>1.0)vdot = 1.0;
        if(vdot<-1.0)vdot = -1.0;
        double theta = acos(vdot);
        std::cout << "theta:" << theta*180/M_PI << std::endl;
        if(theta >= THETA_TRANSITION ){
            theta = THETA_MAX;
        }else{
            double alpha = (THETA_TRANSITION - theta)/THETA_TRANSITION;
            theta = THETA_MAX*(1-alpha);
            alpha_friction = 4.0;
        }
        std::cout << "theta after:" << theta << std::endl;
        torque_track = n_torque_b * theta * TORQUE_PER_RAD;
    }

    //std::cout << "torque_track:\n" << torque_track << std::endl;

    // add the torque due to friction in the world system
    Eigen::Vector3d omega_p = caams::omega_p_p_dot(p, pdot);
    //std::cout << "omega:" << omega_p.norm() << "\n" << omega_p << std::endl;
    torque_friction = omega_p * -TORQUE_PER_RAD_PER_S;
    torque_friction.head<2>() *= alpha_friction;

    //std::cout << "torque_friction:\n" << torque_friction << std::endl;

    // transform world torque to body torque
    Eigen::Vector3d torque_b = torque_track + torque_friction;
    Eigen::Vector2d v_theta;
    v_theta = torque_b.head<2>() * -(THETA_DEFLECTION_MAX/2000);
    theta_x = v_theta(0);
    theta_y = v_theta(1);
    if(theta_x>THETA_DEFLECTION_MAX)theta_x = THETA_DEFLECTION_MAX;
    if(theta_x<-THETA_DEFLECTION_MAX)theta_x = -THETA_DEFLECTION_MAX;
    if(theta_y>THETA_DEFLECTION_MAX)theta_y = THETA_DEFLECTION_MAX;
    if(theta_y<-THETA_DEFLECTION_MAX)theta_y = -THETA_DEFLECTION_MAX;
    torque_b_z = torque_b(2);
    std::cout << "theta_x:" << theta_x << "\n theta_y:" << theta_y << std::endl;

}


void Tracker::ForceAndTorque(
        double dt,
        Eigen::Vector3d &force,  // force in global coordinate space
        Eigen::Vector3d &torque) // torque in body coordinate space
{
    force = Eigen::Vector3d::Zero();

    Eigen::Vector3d nFb(cos(theta_x)*sin(theta_y),
                        -sin(theta_x),
                        cos(theta_x)*cos(theta_y));

    Eigen::Vector3d Fb = nFb * THRUST * m_mass;

    Eigen::Vector3d sp(0,0,-m_radius);

    torque = sp.cross(Fb);
    torque(2) = torque_b_z;
    if(std::isnan(torque(0))) exit(1);

    //std::cout << "torque:\n" << torque << std::endl;
}

void Tracker::setVTrack(Eigen::Vector3d const &vTrack)
{
    Tracker::vTrack = vTrack.normalized();
}

void display(void)
{
    glm::dmat4 glm_cam_rot(1.0);
    glm::dvec3 glm_cam_tranv(0, 0, cam_distance);
    glm::dmat4 glm_cam_tranm = glm::translate(glm::dmat4(1.0),glm_cam_tranv);

    glm::dmat4 A_cam = glm_cam_tranm*glm_cam_rot;


    float dist = 20.0f + cam_distance;

    glm::mat4 A_pers = glm::perspective(
        (float)(FOV_Y*RAD_PER_DEG),
        (float)window_width/(float)window_height,
        CAM_NEAR,
        dist);

    glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);


    glMatrixMode(GL_PROJECTION);
    glLoadMatrixf( glm::value_ptr(A_pers) );

    glMatrixMode(GL_MODELVIEW);
    glm::dmat4 A_cam_inv = glm::inverse(A_cam);
    glLoadMatrixd( glm::value_ptr(A_cam_inv) );

    glm::mat4 view(A_cam_inv);

    tracker_system.draw(A_pers,view);

    Eigen::Vector3d vTrack(cos(track_lat*RAD_PER_DEG)*sin(track_long*RAD_PER_DEG),
                           -sin(track_lat*RAD_PER_DEG),
                           cos(track_lat*RAD_PER_DEG)*cos(track_long*RAD_PER_DEG));
    tracker->setVTrack(vTrack);

    glm::vec3 vTrack0(0, 0, 0);
    glm::vec3 vTrack1 = E2GLM(vTrack);
    vTrack1 *= radius_tracker * 2.0;

    glEnable(GL_DEPTH_TEST);
    glBegin(GL_LINES);
    glColor3f(0.0, 1.0, 0.0);
    glVertex3fv(glm::value_ptr(vTrack0));
    glVertex3fv(glm::value_ptr(vTrack1));
    glEnd();

    for(int i=0;i<4;i++)
        tracker_system.rkIntegrate( DELTA_T/4 );

}

void init_program(void)
{
    Eigen::Vector3d pos_tracker(0,0,0);
    Eigen::Vector3d vel_tracker(0,0,0);
    Eigen::Vector4d p_tracker(1.0, 0, 0, 0);
    Eigen::Vector4d pdot_tracker(0, 0, 0, 0);
    Eigen::Matrix3d Jp_tracker = caams::J_p_cylinder_z_axis(mass_tracker,radius_tracker*0.1,radius_tracker);

    tracker = new Tracker(
                mass_tracker,
                radius_tracker,
                pos_tracker,
                vel_tracker,
                p_tracker,
                pdot_tracker,
                Jp_tracker,
                (char*)"sat_tex.jpg");

    tracker_system.AddBody(tracker);

    glClearColor(0.0, 0.0, 0.0, 0.0);
    glDrawBuffer(GL_BACK);

    font.LoadOutline(
                "/usr/share/fonts/truetype/ubuntu/UbuntuMono-R.ttf",
                FONT_SIZE,Pixel32(0, 0, 0),Pixel32(255, 255, 255), 1.0 );

    uv_sphere_init(64,128);
}

void keydown(SDL_Event *event)
{
    switch(event->key.keysym.sym){
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
        track_long += (double)delta_x_int*0.25;
        track_lat += (double)delta_y_int*0.25;
        if(track_long>180.0)track_long-=360.0;
        if(track_long<-180.0)track_long+=360.0;
        if(track_lat>90.0)track_lat=90.0;
        if(track_lat<-90.0)track_lat=-90.0;
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

void mousewheel(SDL_Event *event)
{
}

int main(int argc, char **argv)
{
    if(SDL_Init(SDL_INIT_VIDEO)!=0){
        SDL_Log("Unable to initialize SDL: %s", SDL_GetError());
        return 1;
    }

    printf("Current video driver:%s\n",SDL_GetCurrentVideoDriver());

    SDL_Window *window;
    window = SDL_CreateWindow(
        "Direction Tracker Tester",
        SDL_WINDOWPOS_UNDEFINED,
        SDL_WINDOWPOS_UNDEFINED,
        800,
        600,
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

    init_program();

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
            case SDL_MOUSEWHEEL:
                mousewheel(&event);
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


