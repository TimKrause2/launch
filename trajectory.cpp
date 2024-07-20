#include <GL/glew.h>
#include "trajectory.h"

Trajectory::Trajectory(int Npoints) :
    Npoints(Npoints),
    positions(3,Npoints),
    velocities(3,Npoints)
{

}

void Trajectory::SetPoint(int i, Body *body)
{
    if(i>=Npoints || i<0)return;
    //positions.sub(~body->m_position,i+1,1);
    //velocities.sub(~body->m_velocity,i+1,1);
    positions.col(i) = body->m_position;
    velocities.col(i) = body->m_velocity;
}

void Trajectory::GetPoint(int i, Body *body){
    if(i>=Npoints || i<0)return;
//    body->m_position = ~positions.sub(1,3,i+1,1);
//    body->m_velocity = ~velocities.sub(1,3,i+1,1);
    body->m_position = positions.col(i);
    body->m_velocity = velocities.col(i);
}

void Trajectory::Draw(int N)
{
    //double* v = positions.data;
    glBegin(GL_LINE_STRIP);
    for(int i=0;i<N;i++){
        glVertex3dv(positions.col(i).data());
        //v+=3;
    }
    glEnd();
}


