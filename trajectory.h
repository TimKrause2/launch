#include <GLES3/gl32.h>
#include <GLES3/gl3ext.h>
#include "caams.hpp"
#include "gsim.h"
#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>

class Trajectory{
    int Npoints;
    Eigen::MatrixXd positions;
    Eigen::MatrixXd velocities;
    GLuint vao;
    GLuint vbuff;
    static GLuint program;
    static int n_objects;
    GLuint mvp_id;
    GLuint color_id;
    void LoadProgram(void);
    void UnLoadProgram(void);
public:
    Trajectory(int Npoints);
    ~Trajectory(void);
    void SetPoint(int i,Body *body);
    void GetPoint(int i,Body *body);
    void Draw(glm::mat4 mvp, glm::vec4 color, int Npoints);

};
