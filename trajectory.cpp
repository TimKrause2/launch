#include "trajectory.h"
#include "shader.hpp"

#include <stdio.h>

GLuint Trajectory::program = 0;
int Trajectory::n_objects = 0;

Trajectory::Trajectory(int Npoints) :
    Npoints(Npoints),
    positions(3,Npoints),
    velocities(3,Npoints)
{
    LoadProgram();

    glGenVertexArrays(1, &vao);
    glBindVertexArray(vao);

    glGenBuffers(1, &vbuff);
    glBindBuffer(GL_ARRAY_BUFFER, vbuff);
    glBufferData(GL_ARRAY_BUFFER, sizeof(float)*3*Npoints, NULL, GL_DYNAMIC_DRAW);

    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, (void*)0);

    glBindVertexArray(0);
}

Trajectory::~Trajectory(void)
{
    UnLoadProgram();
    glDeleteVertexArrays(1, &vao);
    glDeleteBuffers(1, &vbuff);
}

void Trajectory::LoadProgram(void)
{
    if(n_objects==0){
        program = LoadShaders("trajectory.vtx.glsl", "trajectory.fgmt.glsl");
        mvp_id = glGetUniformLocation(program, "MVP");
        color_id = glGetUniformLocation(program, "color");
    }
    n_objects++;
}

void Trajectory::UnLoadProgram(void)
{
    if(n_objects){
        n_objects--;
        if(n_objects==0){
            glDeleteProgram(program);
        }
    }
}

void checkError() {
    GLenum err;
    while((err = glGetError()) != GL_NO_ERROR){
        printf("OpenGL Error:%d\n", err);
    }
}

void Trajectory::SetPoint(int i, Body *body)
{
    if(i>=Npoints || i<0)return;
    positions.col(i) = body->m_position;
    velocities.col(i) = body->m_velocity;
}

void Trajectory::GetPoint(int i, Body *body){
    if(i>=Npoints || i<0)return;
    body->m_position = positions.col(i);
    body->m_velocity = velocities.col(i);
}

void Trajectory::Draw(glm::mat4 mvp, glm::vec4 color, int N)
{
    if(N<2) return;

    // copy the position matrix to the vertex buffer
    glBindBuffer(GL_ARRAY_BUFFER, vbuff);
    glm::vec3 *data = (glm::vec3*)glMapBufferRange(
                GL_ARRAY_BUFFER,
                0, sizeof(glm::vec3)*N,
                GL_MAP_WRITE_BIT|GL_MAP_INVALIDATE_RANGE_BIT);
    for(int i=0;i<N;i++){
        Eigen::Vector3d ve = positions.col(i);
        data[i] = glm::vec3(ve(0), ve(1), ve(2));
    }
    glUnmapBuffer(GL_ARRAY_BUFFER);

    glUseProgram(program);
    glBindVertexArray(vao);
    glUniformMatrix4fv(mvp_id, 1, GL_FALSE, glm::value_ptr(mvp));
    glUniform4fv(color_id, 1, glm::value_ptr(color));

    glDrawArrays(GL_LINE_STRIP, 0, N);

    glBindVertexArray(0);
    glUseProgram(0);
}


