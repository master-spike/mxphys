#include "mxphys/body.h"
#include "mxphys/forces.h"
#include <iostream>

namespace mxphys
{

void body::get_contact_points(body& other, std::vector<contact_point>& out_contact_points) {
    auto cn1 = sat_contact_normal(other);
    auto cn2 = other.sat_contact_normal(*this);
    if (!cn1 || !cn2) return;
    if (cn1->first <= cn2->first) {
        const vec2& norm = cn1->second;
        affine_2d map_to_own = position.inverse() * other.position;
        double threshold = norm.dot(
            *std::max_element(shape.getPoints().begin(), shape.getPoints().end(), [&norm](vec2 const& p1, vec2 const& p2) {
                return norm.dot(p1) < norm.dot(p2);
            })
        );
        std::for_each(other.shape.getPoints().begin(), other.shape.getPoints().end(), [&](const vec2& p) {
            if (norm.dot(map_to_own(p)) <= threshold)
                out_contact_points.emplace_back(other.position(p), position.scale * norm, &other, this);
        });
    }
    if (cn2->first <= cn1->first) {
        const vec2& norm = cn2->second;
        affine_2d map_to_other = other.position.inverse() * position;
        double threshold = norm.dot(
            *std::max_element(other.shape.getPoints().begin(), other.shape.getPoints().end(), [&norm](vec2 const& p1, vec2 const& p2) {
                return norm.dot(p1) < norm.dot(p2);
            })
        );
        std::for_each(shape.getPoints().begin(), shape.getPoints().end(), [&](const vec2& p) {
            if (norm.dot(map_to_other(p)) <= threshold)
                out_contact_points.emplace_back(position(p), other.position.scale * norm, this, &other);
        });
    }
}

void body::apply_impulse(impulse_point imp) {
    auto adj_origin = position.inverse()(imp.origin);
    auto adj_impulse = position.scale.inverse() * imp.impulse;
    angular_frequency += (1.0/moment_of_inertia) * adj_impulse.cross(-adj_origin);
    velocity += (1.0/mass) * imp.impulse;
}

bounding_box body::getBoundingBox() const
{
    return bounding_box{
        getTranslatedPoints()
    };
}

uint64_t body::id_next = 0;

std::optional<std::pair<double,vec2>> body::sat_contact_normal(const body& other) const {
    affine_2d map_to_own = position.inverse() * other.position;

    auto other_points = std::views::transform(other.shape.getPoints(), map_to_own);
    double min_depth = std::numeric_limits<double>::infinity();
    vec2 min_depth_norm = vec2::zerovec();

    for (auto it = shape.getPoints().begin(); it != shape.getPoints().end(); ++it) {
        vec2 normal = it + 1 == shape.getPoints().end() ? *shape.getPoints().begin() - *it: *(it+1) - *it;
        normal = mat2{0, 1, -1, 0} * normal.normalized();
        
        auto min_other = std::min_element(other_points.begin(), other_points.end(), [&normal](vec2 const& p1, vec2 const& p2){
            return normal.dot(p1) < normal.dot(p2);
        });
        auto max_this = std::max_element(shape.getPoints().begin(), shape.getPoints().end(), [&normal](vec2 const& p1, vec2 const& p2){
            return normal.dot(p1) < normal.dot(p2);
        });
        double depth = normal.dot(*max_this) - normal.dot(*min_other);
        
        // No collision!
        if (depth < 0) return std::nullopt;

        if (min_depth > depth) {
            min_depth = depth;
            min_depth_norm = normal;
        }
    }
    return std::make_pair(min_depth, min_depth_norm);
}

}
