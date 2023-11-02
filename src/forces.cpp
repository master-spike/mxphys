
#include "mxphys/forces.h"
#include "mxphys/body.h"

namespace mxphys {

bool contact_point::resolve() {

    for (auto& f : pre_collision_callbacks) {
        f(this);
    }

    double restitution = body1->getElasticity() * body2->getElasticity();
    auto rvel = body2->velocityAt(position) - body1->velocityAt(position);
    
    // objects moving away from eachother at this point
    if (rvel.dot(normal) <= 0) return false;
    
    auto rvel_fcomp = rvel.dot(normal) * normal;
    
    affine_2d p1 = body1->getPos();
    affine_2d p2 = body2->getPos();

    double d1 = normal.cross(position - p1.shift);
    double d2 = normal.cross(position - p2.shift);

    d1 *= d1;
    d2 *= d2;

    double t1 = 1.0 / body1->getMass();
    double t2 = d1  / body1->getMomentOfInertia();
    double t3 = 1.0 / body2->getMass();
    double t4 = d2  / body2->getMomentOfInertia();

    // do not collide objects that are completely unaccelerateable
    if (t1 + t2 + t3 + t4 == 0) return false;

    auto r_impulse = rvel_fcomp * ((1.0 + restitution) / (t1 + t2 + t3 + t4));

    impulse_point impulse_p1{position,  r_impulse};
    impulse_point impulse_p2{position, -r_impulse};

    body1->apply_impulse(impulse_p1);
    body2->apply_impulse(impulse_p2);

    for (auto& f : post_collision_callbacks) {
        f(this);
    }
    return true;
}

}
