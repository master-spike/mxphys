#ifndef MXPHYS_FORCES_H
#define MXPHYS_FORCES_H

#include <functional>
#include <vector>
#include <optional>
#include <utility>

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
    bool resolve();
    bool resolve(event::eventmanager& event_manager);
    contact_point() = delete;
    
    contact_point(vec2 _pos, vec2 _normal, body* _body1, body* _body2) :
        position(_pos), normal(_normal), body1(_body1), body2(_body2) {}
private:
    std::optional<impulse_point> calculate_impulse_point();
};

}

#endif
