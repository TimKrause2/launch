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
#include "font.h"
#include "orbit.h"
#include "gsim.h"
#include "uv_sphere.h"

bool mousing=false;
bool up_key=false;
bool down_key=false;
bool left_key=false;
bool right_key=false;
int last_x;
int last_y;
int modifiers;
Font font;

#define RAD_PER_DEG (M_PI/180.0)
#define PLANET_SCALE 10.0
#define E_MAX 0.9
#define FOV_Y  25.0
#define FONT_SIZE 32


OrbitalParams params0;
double view_lat=0.0;
double view_long=0.0;
double mass_earth = 5.97219e24;
double radius_earth = 6.371e6*PLANET_SCALE;
double mass_moon = 7.3477e22;
double radius_moon = 1.7371e6*PLANET_SCALE;
double T_earth_moon = (28*24*60*60);


int window_width;
int window_height;

Body* body_earth;
Body* body_moon;

System earth_moon_system;

enum edit_mode_t
{
    EDIT_MODE_NONE,
    EDIT_MODE_ECCENTRICITY,
    EDIT_MODE_INCLINATION,
    EDIT_MODE_LONGITUDE,
    EDIT_MODE_PERIAPSIS,
    EDIT_MODE_ANOMOLY,
    EDIT_MODE_VIEW,
    EDIT_MODE_ANIMATE
};

edit_mode_t edit_mode = EDIT_MODE_NONE;

void display(void)
{
    // calculate the radius of the scene
    double r_scene = params0.a*mass_earth/(mass_earth+mass_moon)*(1+params0.e);

    // calculate the distance of the camera
    double r_camera = r_scene/sin(FOV_Y*RAD_PER_DEG/2.0);

    // calculate the near and far planes
    float near = r_camera - r_scene*1.5;
    float far = r_camera + r_scene*1.5;

    // calculate the orientation of the camera
    glm::dmat4 A_lat = glm::rotate(
                glm::dmat4(1.0),
                (90.0-view_lat)*RAD_PER_DEG,
                glm::dvec3(1.0,0.0,0.0));
    glm::dmat4 A_long = glm::rotate(
                glm::dmat4(1.0),
                -view_long*RAD_PER_DEG,
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

    orbit_draw(params0,mass_earth,mass_moon);

    if(edit_mode!=EDIT_MODE_ANIMATE){
        orbit_initialize(params0,
                         mass_earth,mass_moon,
                         body_earth->m_position,
                         body_moon->m_position,
                         body_earth->m_velocity,
                         body_moon->m_velocity);
    }

    earth_moon_system.draw(A_pers, view);

    if(edit_mode==EDIT_MODE_ANIMATE){
        for(int n=0;n<24;n++){
            earth_moon_system.rkIntegrate(60.0);
        }
    }

    ModEquOrbitalParams eq;
    eq.FromOrbitalParams(params0);
    OrbitalParams params1;
    eq.ToOrbitalParams(params1);

    Eigen::Vector3d r;
    Eigen::Vector3d v;
    orbit_state(params0, mass_earth, mass_moon, r, v);
    OrbitalParams params_rv;
    orbital_elements(r, v, mass_earth, mass_moon, params_rv);

    double x=FONT_SIZE;
    double y=FONT_SIZE;

    font.Printf(x,y,"a:%le a1:%le",params0.a,params1.a);
    y+=FONT_SIZE;
    font.Printf(x,y,"e:%le e1:%le",params0.e,params1.e);
    y+=FONT_SIZE;
    font.Printf(x,y,"i:% 7.2lf i1:% 7.2lf",params0.i/RAD_PER_DEG,params1.i/RAD_PER_DEG);
    y+=FONT_SIZE;
    font.Printf(x,y,"O:% 7.2lf O1:% 7.2lf",params0.Omega/RAD_PER_DEG,params1.Omega/RAD_PER_DEG);
    y+=FONT_SIZE;
    font.Printf(x,y,"o:% 7.2lf o1:% 7.2lf",params0.omega/RAD_PER_DEG,params1.omega/RAD_PER_DEG);
    y+=FONT_SIZE;
    font.Printf(x,y,"n:% 7.2lf n1:% 7.2lf",params0.nu/RAD_PER_DEG,params1.nu/RAD_PER_DEG);

    y+=2*FONT_SIZE;
    font.Printf(x,y,"a_rv:%le", params_rv.a);
    y+=FONT_SIZE;
    font.Printf(x,y,"e_rv:%le", params_rv.e);
    y+=FONT_SIZE;
    font.Printf(x,y,"i_rv:% 7.2lf", params_rv.i/RAD_PER_DEG);
    y+=FONT_SIZE;
    font.Printf(x,y,"O_rv:% 7.2lf", params_rv.Omega/RAD_PER_DEG);
    y+=FONT_SIZE;
    font.Printf(x,y,"o_rv:% 7.2lf", params_rv.omega/RAD_PER_DEG);
    y+=FONT_SIZE;
    font.Printf(x,y,"n_rv:% 7.2lf", params_rv.nu/RAD_PER_DEG);



    x = window_width/2;
    y = FONT_SIZE;

    font.Printf(x,y,"p:%le",eq.p);
    y+=FONT_SIZE;
    font.Printf(x,y,"f:%le",eq.f);
    y+=FONT_SIZE;
    font.Printf(x,y,"g:%le",eq.g);
    y+=FONT_SIZE;
    font.Printf(x,y,"h:%le",eq.h);
    y+=FONT_SIZE;
    font.Printf(x,y,"k:%le",eq.k);
    y+=FONT_SIZE;
    font.Printf(x,y,"L:% 7.2lf",eq.L/RAD_PER_DEG);



}

void init_program(void)
{
    params0.a = semi_major(T_earth_moon,mass_earth,mass_moon);
    params0.e = 0.0;
    params0.i = 0.0;
    params0.Omega = 0.0;
    params0.omega = 0.0;
    params0.nu = 0.0;

    Eigen::Vector3d pos_earth;
    Eigen::Vector3d pos_moon;
    Eigen::Vector3d vel_earth;
    Eigen::Vector3d vel_moon;

    orbit_initialize(params0,
                     mass_earth,mass_moon,
                     pos_earth,pos_moon,
                     vel_earth,vel_moon);

    Eigen::Vector4d p_identity;
    p_identity << 1.0, 0.0, 0.0, 0.0;

    Eigen::Vector4d p_earth = p_identity;
    Eigen::Vector4d p_moon = p_identity;
    Eigen::Vector3d omega_earth;
    omega_earth << 0.0, 0.0, 2.0*M_PI/(24*60*60);
    Eigen::Vector4d pdot_earth
        = caams::p_dot_omega(p_earth, omega_earth);
    Eigen::Vector3d omega_moon;
    omega_moon << 0.0, 0.0, 2.0*M_PI/T_earth_moon;
    Eigen::Vector4d pdot_moon
        = caams::p_dot_omega(p_moon, omega_moon);
    Eigen::Matrix3d Jp_earth = caams::J_p_sphere(mass_earth,radius_earth);
    Eigen::Matrix3d Jp_moon = caams::J_p_sphere(mass_moon,radius_moon);

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

    earth_moon_system.AddBody(body_earth);
    earth_moon_system.AddBody(body_moon);

    glClearColor(0.0,0.0,0.0,0.0);
    glDrawBuffer(GL_BACK);

    uv_sphere_init(64,128);

    font.LoadOutline(
        "/usr/share/fonts/truetype/ubuntu/UbuntuMono-R.ttf",
        FONT_SIZE,Pixel32(0, 0, 0),Pixel32(255, 255, 255), 1.0 );
}

void keydown(SDL_Event *event)
{
    switch(event->key.keysym.sym){
    case SDLK_e:
        edit_mode = EDIT_MODE_ECCENTRICITY;
        break;
    case SDLK_i:
        edit_mode = EDIT_MODE_INCLINATION;
        break;
    case SDLK_o:
        edit_mode = EDIT_MODE_LONGITUDE;
        break;
    case SDLK_p:
        edit_mode = EDIT_MODE_PERIAPSIS;
        break;
    case SDLK_n:
        edit_mode = EDIT_MODE_ANOMOLY;
        break;
    case SDLK_v:
        edit_mode = EDIT_MODE_VIEW;
        break;
    case SDLK_SPACE:
        if(edit_mode == EDIT_MODE_NONE){
            edit_mode = EDIT_MODE_ANIMATE;
        }else if(edit_mode == EDIT_MODE_ANIMATE){
            edit_mode = EDIT_MODE_NONE;
        }
        break;
    case SDLK_UP:
        up_key = true;
        break;
    case SDLK_DOWN:
        down_key = true;
        break;
    case SDLK_LEFT:
        left_key = true;
        break;
    case SDLK_RIGHT:
        right_key = true;
        break;
    }
}

void keyup(SDL_Event *event)
{
    switch(event->key.keysym.sym){
    case SDLK_e:
    case SDLK_i:
    case SDLK_o:
    case SDLK_p:
    case SDLK_n:
    case SDLK_v:
        edit_mode = EDIT_MODE_NONE;
        break;
    case SDLK_SPACE:
        break;
    case SDLK_UP:
        up_key = false;
        break;
    case SDLK_DOWN:
        down_key = false;
        break;
    case SDLK_LEFT:
        left_key = false;
        break;
    case SDLK_RIGHT:
        right_key = false;
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
        if( delta_y_int != 0 || delta_x_int != 0 ){
            switch(edit_mode){
            case EDIT_MODE_NONE:
                break;
            case EDIT_MODE_ECCENTRICITY:
                params0.e += (double)delta_x_int*0.001;
                if(params0.e<0.0)params0.e=0.0;
                if(params0.e>E_MAX)params0.e=E_MAX;
                break;
            case EDIT_MODE_INCLINATION:
                params0.i += (double)delta_x_int*0.25*RAD_PER_DEG;
                if(params0.i<0.0)params0.i=0.0;
                if(params0.i>M_PI)params0.i=M_PI;
                break;
            case EDIT_MODE_LONGITUDE:
                params0.Omega += (double)delta_x_int*0.25*RAD_PER_DEG;
                if(params0.Omega<-M_PI)params0.Omega+=2*M_PI;
                if(params0.Omega>M_PI)params0.Omega-=2*M_PI;
                break;
            case EDIT_MODE_PERIAPSIS:
                params0.omega += (double)delta_x_int*0.25*RAD_PER_DEG;
                if(params0.omega<-M_PI)params0.omega+=2*M_PI;
                if(params0.omega>M_PI)params0.omega-=2*M_PI;
                break;
            case EDIT_MODE_ANOMOLY:
                params0.nu += (double)delta_x_int*0.25*RAD_PER_DEG;
                if(params0.nu<-M_PI)params0.nu+=2*M_PI;
                if(params0.nu>M_PI)params0.nu-=2*M_PI;
                break;
            case EDIT_MODE_VIEW:
                view_long += (double)delta_x_int*0.25;
                view_lat += (double)delta_y_int*0.25;
                if(view_long>180.0)view_long-=360.0;
                if(view_long<-180.0)view_long+=360.0;
                if(view_lat>90.0)view_lat=90.0;
                if(view_lat<-90.0)view_lat=-90.0;
                break;
            }
        }
        if(params0.i==0.0 || params0.i==M_PI){
            params0.Omega = 0.0;
        }
        if(params0.e==0.0){
            params0.omega = 0.0;
        }
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
        "Orbit Visualizer",
        SDL_WINDOWPOS_UNDEFINED,
        SDL_WINDOWPOS_UNDEFINED,
        1500,
        1000,
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


