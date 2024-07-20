#ifndef LAUNCH_H
#define LAUNCH_H

#include "esfont.h"
#include "orbit.h"
#include "gsim.h"
#include "uv_sphere.h"
#include "satellite.h"
#include "spheroid.h"

struct LaunchData
{
    LaunchData();
    void init(void);
    void display(void);
    void update_view_direction(int dx, int dy);
    void update_view_distance(int dr);
    void update_time_factor(int dt);
    void update_viewport(int width, int height);
    FreeTypeFont font;
    double view_lat;
    double view_long;
    int window_width;
    int window_height;
    double time_factor;
    float cam_distance;
    Spheroid* body_earth;
    Satellite* body_sat;
    System earth_sat_system;
    Eigen::Vector3d n_h_launch;

    OrbitalParams params_new;
    OrbitalParams params_dl;
};



#endif
