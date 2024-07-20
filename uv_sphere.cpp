#include "uv_sphere.h"
#include "shader.hpp"
#include <glm/gtc/type_ptr.hpp>

struct v_position
{
    float x,y,z;
};

struct v_uv
{
    float u,v;
};

struct triangle
{
    unsigned short v0, v1, v2;
};

bool initialized=false;
GLuint VertexArrayID;
GLuint ProgramID;
GLuint MatrixID;
GLuint TextureID;
int Ntri;

void uv_sphere_init(int Nlat, int Nlong)
{
    if(initialized)return;

    int Nvertices = (Nlat+1)*(Nlong+1);
    int Nfaces = Nlat*Nlong;
    Ntri = Nfaces*2;

    v_position *v_pos = new v_position[Nvertices];
    v_uv *v_uvs = new v_uv[Nvertices];

    v_position *v_pos_it = v_pos;
    v_uv *v_uv_it = v_uvs;

    for(int lat=0;lat<=Nlat;lat++){
        float phi = (float)lat*M_PI/Nlat;
        float v = (float)lat/Nlat;
        for(int lo=0;lo<=Nlong;lo++){
            float theta = (float)lo*2.0*M_PI/Nlong-M_PI;
            float u = (float)lo/Nlong;

            v_pos_it->x = cosf(theta)*sinf(phi);
            v_pos_it->y = sinf(theta)*sinf(phi);
            v_pos_it->z = cosf(phi);

            v_uv_it->u = u;
            v_uv_it->v = v;

            v_pos_it++;
            v_uv_it++;
        }
    }

    triangle *indices = new triangle[Ntri];
    triangle *indices_it = indices;

    for(int lat=0;lat<Nlat;lat++){
        unsigned short q_v0 = lat*(Nlong+1);
        unsigned short q_v1 = q_v0 + (Nlong+1);
        unsigned short q_v2 = q_v1 + 1;
        unsigned short q_v3 = q_v0 + 1;
        for(int lo=0;lo<Nlong;lo++){
            indices_it->v0 = q_v0;
            indices_it->v1 = q_v1;
            indices_it->v2 = q_v2;
            indices_it++;
            indices_it->v0 = q_v0;
            indices_it->v1 = q_v2;
            indices_it->v2 = q_v3;
            indices_it++;
            q_v0++;
            q_v1++;
            q_v2++;
            q_v3++;
        }
    }

    glGenVertexArrays(1, &VertexArrayID);
    glBindVertexArray(VertexArrayID);

    GLuint vertexbuffer;
    glGenBuffers(1, &vertexbuffer);
    glBindBuffer(GL_ARRAY_BUFFER, vertexbuffer);
    glBufferData(GL_ARRAY_BUFFER, sizeof(v_position)*Nvertices, v_pos, GL_STATIC_DRAW);

    GLuint uvbuffer;
    glGenBuffers(1, &uvbuffer);
    glBindBuffer(GL_ARRAY_BUFFER, uvbuffer);
    glBufferData(GL_ARRAY_BUFFER, sizeof(v_uv)*Nvertices, v_uvs, GL_STATIC_DRAW);

    GLuint elementbuffer;
    glGenBuffers(1, &elementbuffer);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, elementbuffer);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(triangle)*Ntri, indices, GL_STATIC_DRAW);

    glEnableVertexAttribArray(0);
    glBindBuffer(GL_ARRAY_BUFFER, vertexbuffer);
    glVertexAttribPointer(
        0,                  // attribute. No particular reason for 0, but must match the layout in the shader.
        3,                  // size
        GL_FLOAT,           // type
        GL_FALSE,           // normalized?
        0,                  // stride
        (void*)0            // array buffer offset
    );

    glEnableVertexAttribArray(1);
    glBindBuffer(GL_ARRAY_BUFFER, uvbuffer);
    glVertexAttribPointer(
        1,                                // attribute. No particular reason for 1, but must match the layout in the shader.
        2,                                // size : U+V => 2
        GL_FLOAT,                         // type
        GL_FALSE,                         // normalized?
        0,                                // stride
        (void*)0                          // array buffer offset
    );

    ProgramID = LoadShaders( "TransformVertexShader.vertexshader", "TextureFragmentShader.fragmentshader" );
    MatrixID = glGetUniformLocation(ProgramID, "MVP");
    TextureID  = glGetUniformLocation(ProgramID, "myTextureSampler");

    glBindVertexArray(0);

    delete [] v_pos;
    delete [] v_uvs;
    delete [] indices;

    initialized = true;
}

void uv_sphere_draw(glm::mat4 mvp, GLuint tex_buff)
{
    if(!initialized)return;

    glBindVertexArray(VertexArrayID);
    glUseProgram(ProgramID);
    glUniformMatrix4fv(MatrixID, 1, GL_FALSE, glm::value_ptr(mvp));

    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, tex_buff);
    glUniform1i(TextureID, 0);

    glDrawElements(
                GL_TRIANGLES,
                Ntri*3,
                GL_UNSIGNED_SHORT,
                (void*)0);

    glBindVertexArray(0);
    glUseProgram(0);
}


