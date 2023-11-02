#include "mxphys/body.h"
#include "mxphys/forces.h"

namespace mxphys
{

void body::get_contact_points(body& other, std::vector<contact_point>& out_contact_points) {
    for (auto p : shape.getPoints()) {
        auto ptrue = position(p);
        auto pomap = other.position.inverse()(ptrue);
        auto normal = other.position.scale * other.shape.normalAt(pomap);
        if (normal == vec2{0.0, 0.0}) continue;
        out_contact_points.emplace_back(ptrue, normal, this, &other);
    }
    for (auto q : other.shape.getPoints()) {
        auto qtrue = other.position(q);
        auto qimap = position.inverse()(qtrue);
        auto normal = position.scale * shape.normalAt(qimap);
        if (normal == vec2{0.0, 0.0}) continue;
        out_contact_points.emplace_back(qtrue, normal, &other, this);
    }
}

void body::apply_impulse(impulse_point imp) {
    auto adj_origin = position.inverse()(imp.origin);
    auto adj_impulse = position.scale.inverse() * imp.impulse;
    angular_frequency += (1.0/moment_of_inertia) * adj_impulse.cross(-adj_origin);
    velocity += (1.0/mass) * imp.impulse;
}

uint64_t body::id_next = 0;
    
}
