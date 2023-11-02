
#include <cmath>
#include <memory>

#include "mxphys/event/eventmanager.h"
#include "mxphys/event/collision_event.h"
#include "mxphys/forces.h"
#include "mxphys/body.h"

namespace mxphys {

bool contact_point::resolve() {

    auto impulse_p1_opt = calculate_impulse_point();
    if (!impulse_p1_opt) return false;

    body1->apply_impulse(impulse_p1_opt.value());
    body2->apply_impulse(impulse_point{impulse_p1_opt.value().origin, -impulse_p1_opt.value().impulse});

    return true;
}

bool contact_point::resolve(event::eventmanager& event_manager) {
    auto impulse_p1_opt = calculate_impulse_point();
    if (!impulse_p1_opt) return false;

    auto ids = std::make_pair(body1->getID(), body2->getID());
    auto p_vels = std::make_pair(body1->getVelocity(), body2->getVelocity());
    auto p_afs = std::make_pair(body1->getAngularFreq(), body2->getAngularFreq());

    double mag = impulse_p1_opt.value().impulse.dot(impulse_p1_opt.value().impulse);
    mag = std::sqrt(mag);

    body1->apply_impulse(impulse_p1_opt.value());
    body2->apply_impulse(impulse_point{impulse_p1_opt.value().origin, -impulse_p1_opt.value().impulse});

    auto f_vels = std::make_pair(body1->getVelocity(), body2->getVelocity());
    auto f_afs = std::make_pair(body1->getAngularFreq(), body2->getAngularFreq());

    event_manager.add_event(
        event::event_type::COLLISION,
        std::make_unique<event::collision_event>(
            ids, std::make_pair(body1->getMass(), body2->getMass()),
            std::make_pair(body1->getMomentOfInertia(), body2->getMomentOfInertia()),
            p_vels, f_vels, p_afs, f_afs, mag, impulse_p1_opt.value().origin, impulse_p1_opt.value().impulse,
            body1->getElasticity() * body2->getElasticity()
        )
    );

    return true;
}

std::optional<impulse_point> contact_point::calculate_impulse_point()
{
    double restitution = body1->getElasticity() * body2->getElasticity();
    auto rvel = body2->velocityAt(position) - body1->velocityAt(position);
    
    // objects moving away from eachother at this point
    if (rvel.dot(normal) <= 0) return std::nullopt;
    
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
    if (t1 + t2 + t3 + t4 == 0) return std::nullopt;

    auto r_impulse = rvel_fcomp * ((1.0 + restitution) / (t1 + t2 + t3 + t4));

    return impulse_point{position, r_impulse};
}

}
