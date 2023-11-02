#ifndef MXPHYS_EVENT_COLLISION_EVENT_H
#define MXPHYS_EVENT_COLLISION_EVENT_H

#include <utility>

#include "mxphys/types.h"
#include "mxphys/forces.h"

namespace mxphys::event {

struct collision_event : public event_data {
    std::pair<uint64_t, uint64_t> ids;
    std::pair<double, double> masses;
    std::pair<double, double> moments_of_inertia;
    std::pair<vec2, vec2> vs_init;
    std::pair<vec2, vec2> vs_final;
    std::pair<double, double> afs_init;
    std::pair<double, double> afs_final;
    double impulse_mag;
    vec2 point;
    vec2 normal;
    double restitution;
    collision_event(
        auto _ids,
        auto _masses,
        auto _moments_of_inertia,
        auto _vs_init,
        auto _vs_final,
        auto _afs_init,
        auto _afs_final,
        auto _impulse_mag,
        auto _point,
        auto _normal,
        auto _rest
    ) : ids(_ids), masses(_masses), moments_of_inertia(_moments_of_inertia),
        vs_init(_vs_init), vs_final(_vs_final), afs_init(_afs_init),afs_final(_afs_final),
        impulse_mag(_impulse_mag), point(_point), normal(_normal), restitution(_rest) {}
};

}

#endif


