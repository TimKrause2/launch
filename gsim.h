#ifndef GSIM_H
#define GSIM_H

#include "caams.hpp"
//#define GL_GLEXT_PROTOTYPES 1
//#include <GL/gl.h>
//#include <GL/glext.h>
#include <GLES3/gl32.h>
#include <GLES3/gl3ext.h>

//#include <GL/glew.h>
#include <glm/glm.hpp>

// inserting into a list
#include <list>

#define G_gravity 6.67384e-11

#define SPHERE_N_LATITUDE 64
#define SPHERE_N_LONGITUDE 128

class System;

struct State
{
    Eigen::Vector3d m_position;
    Eigen::Vector3d m_velocity;
    Eigen::Vector4d p;
    Eigen::Vector4d pdot;
    State(
        Eigen::Vector3d m_position,
        Eigen::Vector3d m_velocity,
        Eigen::Vector4d p,
        Eigen::Vector4d pdot):
        m_position(m_position),
        m_velocity(m_velocity),
        p(p),
        pdot(pdot){}
};

class Body
{
private:
    GLuint texName;
    std::list<State> states;
    Body operator=(State const &state){
        m_position = state.m_position;
        m_velocity = state.m_velocity;
        p = state.p;
        pdot = state.pdot;
        return *this;
    }
    bool enabled;
public:
	double m_mass;
    double m_radius;
    Eigen::Vector3d m_position;
    Eigen::Vector3d m_velocity;
    Eigen::Vector4d p; // Euler parameters
    Eigen::Vector4d pdot; // time derivative
    Eigen::Vector4d rk_p;
    Eigen::Vector4d rk_pdot;
    Eigen::Vector4d rk_pddot; // second derivative of p
    Eigen::Matrix3d Jp; // Inertia tensor
    Eigen::Vector3d m_rk_position;
    Eigen::Vector3d m_rk_velocity;
    Eigen::Vector3d m_rk_acceleration;
    Eigen::Vector3d ex_position0;
    Eigen::Vector3d ex_position1;
    Eigen::Vector3d ex_velocity0;
    Eigen::Vector3d ex_velocity1;
    Eigen::Matrix3x4d m_kr;
    Eigen::Matrix3x4d m_kv;
    Eigen::Matrix4d k_p_dot;
    Eigen::Matrix4d k_p_ddot;
    Body(Body &b);
    Body(double mass,
         double radius,
         Eigen::Vector3d &position,
         Eigen::Vector3d &velocity,
         Eigen::Vector4d &p,
         Eigen::Vector4d &pdot,
         Eigen::Matrix3d &Jp,
         char* tex_image);

    void PushState(void);
    void PopState(void);
    void RestoreState(void);
    bool Enabled(void);
    void Enable(void);
    void Disable(void);

    void draw(glm::mat4 proj, glm::mat4 view);

    virtual void Prepare(void);
    virtual void Update(double dt);
    virtual void ForceAndTorque(
            double dt,
            Eigen::Vector3d &force,   // force in global coordinate space
            Eigen::Vector3d &torque); // torque in body coordinate space
    virtual double TimeStep(void);
    double CurvatureTimeStep(void);
    virtual Eigen::Vector3d rk_acceleration(
            Eigen::Vector3d const &rk_position,
            Eigen::Vector3d const &rk_velocity);

    void update_rotation( double dt );
};

class System
{
private:
    std::list<Body*> m_bodies;
public:
	void AddBody(Body* p_body);
    void rkIntegrate(double p_dt_total );
    void extrapolate(double t);
    void origin( Eigen::Vector3d r0 );
    void draw(glm::mat4 proj, glm::mat4 view);
    void PushState(void);
    void PopState(void);
    void RestoreState(void);
    Eigen::Vector3d accelerationAtPoint(Eigen::Vector3d p);
private:
    void rkPrepare(void);
    void rkUpdate(double dt);
    void rkAccelerations( double dt );
    void ortho_p_dot(Eigen::Vector4d const &p, Eigen::Vector4d &pdot);
	void rkPhase1Positions( void );
	void rkPhase2Positions( double p_dt );
	void rkPhase3Positions( double p_dt );
	void rkPhase4Positions( double p_dt );
	void rkPhase1Integrate( double p_dt );
	void rkPhase2Integrate( double p_dt );
	void rkPhase3Integrate( double p_dt );
	void rkPhase4Integrate( double p_dt );
    double rkTimeStep( void );
};

#endif
