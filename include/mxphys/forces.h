#ifndef MXPHYS_FORCES_H
#define MXPHYS_FORCES_H

#include <functional>
#include <vector>

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
    std::vector<std::function<void(contact_point*)>> pre_collision_callbacks;
    std::vector<std::function<void(contact_point*)>> post_collision_callbacks;
    bool resolve();
    contact_point() = delete;
    
    contact_point(vec2 _pos, vec2 _normal, body* _body1, body* _body2) :
        position(_pos), normal(_normal), body1(_body1), body2(_body2) {}
    
    contact_point(vec2 _pos, vec2 _normal, body* _body1, body* _body2,
        std::initializer_list<std::function<void(contact_point*)>> _pre_callbacks,
        std::initializer_list<std::function<void(contact_point*)>> _post_callbacks) :
        position(_pos), normal(_normal), body1(_body1), body2(_body2),
        pre_collision_callbacks(_pre_callbacks), post_collision_callbacks(_post_callbacks) {}
};

}

#endif
