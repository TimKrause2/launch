#pragma once

#include <glm/glm.hpp>
#include "gsim.h"

void orbit_draw_init(void);
void orbit_draw(glm::mat4 mvp,
                Body *primary_body, glm::vec4 primary_color,
                Body *secondary_body, glm::vec4 secondary_color);
