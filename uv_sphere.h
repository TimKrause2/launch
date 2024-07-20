//#include <GL/glew.h>
//#define GL_GLEXT_PROTOTYPES 1
//#include <GL/gl.h>
//#include <GL/glext.h>
#include <GLES3/gl32.h>
#include <GLES3/gl3ext.h>

#include <glm/glm.hpp>

void uv_sphere_init(int Nlat, int Nlong);
void uv_sphere_draw(glm::mat4 mvp, GLuint tex_buff);

