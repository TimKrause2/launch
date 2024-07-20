#include "caams.hpp"
#include <GL/gl.h>
#include <list>

#define GM 6.67384e-11

#define SPHERE_N_LATITUDE 64
#define SPHERE_N_LONGITUDE 128

class Body
{
public:
    GLuint         texName;
    double         mass;
    double         radius;
    caams::matrix  Jp;
    Body(double mass, double radius, char* tex_file);
};

class State
{
public:
    Body* body;
    caams::matrix position;
    caams::matrix velocity;
    caams::matrix p; // rotational postion
    caams::matrix pdot; // rotational change in position with time
    caams::matrix a_thrust; // thrust acceleration
    caams::matrix r_thrust; // point of application of the thrust
    State(Body *body);
    State(State &state);
    void draw(void);
    caams::matrix rk_position;
    caams::matrix rk_acceleration;
    caams::matrix m_kr;
    caams::matrix m_kv;
    void IntegrateRotation(double dt);
};

class System
{
private:
    std::list<State> bodies;
public:
    void AddBody(Body *body);
    State& GetState(Body *body);
    void rkIntegrate(double dt);
    void draw(void);
    System Duplicate(void);
private:
    void rkAccelerations( void );
    void rkPhase1Positions( void );
    void rkPhase2Positions( double p_dt );
    void rkPhase3Positions( double p_dt );
    void rkPhase4Positions( double p_dt );
    void rkPhase1Integrate( double p_dt );
    void rkPhase2Integrate( double p_dt );
    void rkPhase3Integrate( double p_dt );
    void rkPhase4Integrate( double p_dt );
    void update_rotations( double dt );
};

