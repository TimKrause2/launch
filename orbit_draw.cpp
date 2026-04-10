#include "orbit_draw.h"
#include "orbit.h"
#include "shader.hpp"
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtx/transform.hpp>
#include <memory>
#include <cmath>

#define  N_VERTICES 1024

static bool initialized = false;
static GLuint program;
static GLuint vao;
static GLuint vbuff;
static GLuint mvp_id;
static GLuint color_id;

void orbit_draw_init(void)
{
    if(initialized) return;
    // load the program
    program = LoadShaders("trajectory.vtx.glsl", "trajectory.fgmt.glsl");
    mvp_id = glGetUniformLocation(program, "MVP");
    color_id = glGetUniformLocation(program, "color");

    // initialize the vertices to a unit radius circle in
    // the xy plane
    std::unique_ptr<glm::vec3[]> vertices
            = std::make_unique<glm::vec3[]>(N_VERTICES);
    for(int i=0;i<N_VERTICES;i++){
        float theta = (float)i*2.0f*M_PIf/N_VERTICES;
        vertices[i] = glm::vec3(cos(theta), sin(theta), 0.0f);
    }

    // create the vertex array object and initialize the vertex buffer
    glGenVertexArrays(1, &vao);
    glBindVertexArray(vao);

    glGenBuffers(1, &vbuff);
    glBindBuffer(GL_ARRAY_BUFFER, vbuff);
    glBufferData(GL_ARRAY_BUFFER, sizeof(glm::vec3)*N_VERTICES, vertices.get(), GL_STATIC_DRAW);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, (void*)0);

    glBindVertexArray(0);

    initialized = true;
}

void orbit_draw(glm::mat4 mvp,
                Body *p_body, glm::vec4 p_color,
                Body *s_body, glm::vec4 s_color)
{
    if(!initialized) return;

    // calculate the orbital elements
    Eigen::Vector3d r = p_body->m_position - s_body->m_position;
    Eigen::Vector3d v = p_body->m_velocity - s_body->m_velocity;
    OrbitalParams params;
    orbital_elements(r, v, p_body->m_mass, s_body->m_mass, params);

    // calculate the bary center
    Eigen::Vector3d bc =
            (p_body->m_position*p_body->m_mass + s_body->m_position*s_body->m_mass)/
            (p_body->m_mass + s_body->m_mass);

    // calculate the scale factors and center locations for the orbits
    double mass_pri = p_body->m_mass;
    double mass_sec = s_body->m_mass;
    double as = params.a*mass_pri/(mass_pri+mass_sec);
    double ap = params.a*mass_sec/(mass_pri+mass_sec);
    double alpha = sqrt(1.0-params.e*params.e);
    double bs = as*alpha;
    double bp = ap*alpha;
    double cs = as*params.e;
    double cp = ap*params.e;

    // initialize the scale matricis
    glm::vec3 scale_p(ap, bp, 1.0f);
    glm::mat4 Mscale_p = glm::scale(scale_p);
    glm::vec3 scale_s(as, bs, 1.0f);
    glm::mat4 Mscale_s = glm::scale(scale_s);

    // initialize the translation matricis
    glm::vec3 t_cp(cp, 0.0f, 0.0f);
    glm::mat4 Mt_cp = glm::translate(t_cp);
    glm::vec3 t_cs(-cs, 0.0f, 0.0f);
    glm::mat4 Mt_cs = glm::translate(t_cs);
    glm::vec3 t_bc(bc(0), bc(1), bc(2));
    glm::mat4 Mt_bc = glm::translate(t_bc);

    // initialize the orbit basis matrix
    Eigen::Matrix3d Mbasis3_e = orbit_basis(params);
    Eigen::Matrix4d Mbasis_e = Eigen::Matrix4d::Identity();
    Mbasis_e.block<3,3>(0,0) = Mbasis3_e;
    glm::mat4 Mbasis;
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            Mbasis[j][i] = Mbasis_e(i, j); // Note the [j][i] swap
        }
    }

    // initialize the mvp matricis
    glm::mat4 mvp_p = mvp*Mt_bc*Mbasis*Mt_cp*Mscale_p;
    glm::mat4 mvp_s = mvp*Mt_bc*Mbasis*Mt_cs*Mscale_s;

    // Use the program to draw the orbits
    glUseProgram(program);
    glBindVertexArray(vao);
    glUniformMatrix4fv(mvp_id, 1, GL_FALSE, glm::value_ptr(mvp_p));
    glUniform4fv(color_id, 1, glm::value_ptr(p_color));
    glDrawArrays(GL_LINE_LOOP, 0, N_VERTICES);

    glUniformMatrix4fv(mvp_id, 1, GL_FALSE, glm::value_ptr(mvp_s));
    glUniform4fv(color_id, 1, glm::value_ptr(s_color));
    glDrawArrays(GL_LINE_LOOP, 0, N_VERTICES);

    glUseProgram(0);
    glBindVertexArray(0);
}
