#ifndef MXPHYS_FORCES_H
#define MXPHYS_FORCES_H

#include "types.h"

namespace mxphys {

struct impulse_point {
    vec2 origin;
    vec2 impulse;
    impulse_point() = delete;
    impulse_point(vec2 _origin, vec2 _impulse) : origin(_origin), impulse(_impulse) {}
};

struct contact_point {
    vec2 position;
    vec2 normal;
    body* body1;
    body* body2;
};

}

#endif
