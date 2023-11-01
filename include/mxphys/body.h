#ifndef MXPHYS_BODY_H
#define MXPHYS_BODY_H

#include <utility>
#include <vector>
#include <functional>
#include <limits>
#include <ranges>

#include "types.h"
#include "polygon.h"
#include "forces.h"

namespace mxphys {

class body {
private:
    double mass;
    double moment_of_inertia;
    vec2 velocity;
    double angular_frequency;
    double elasticity;

    // position of center of mass
    affine_2d position;
    convex_polygon shape;
    bool collision_iteration(body& other, const std::vector<std::pair<vec2,vec2>>& contact_points) {
        double c_elastic = elasticity * other.elasticity;
        for (auto const& [p, norm] : contact_points) {
            auto tang = mat2{0,-1.0,1.0,0} * norm;
            auto rvel = other.velocityAt(p) - velocityAt(p);
            if (rvel.dot(norm) <= 0) continue;
            auto rvel_fcomp = rvel.dot(norm) * norm;
            // min_impulse / m1 + min_impulse / m2 == rvel
            // min_impulse(1/m1 + 1/m2) == rvel
            // min_impulse == rvel(1 / (1/m1 + 1/m2))
            double d1 = norm.cross(p - position.shift);
            double d2 = norm.cross(p - other.position.shift);
            d1 = d1*d1;
            d2 = d2*d2;

            // this is the minimum impulse applied at the point to ensure that the points on each body will move away relative to eachother
            // on the next tick

            double t1 = mass > std::numeric_limits<double>::max() ? 0.0 : 1.0 / mass;
            double t2 = other.mass > std::numeric_limits<double>::infinity() ? 0.0 : 1.0 / other.mass;
            double t3 = moment_of_inertia > std::numeric_limits<double>::infinity() ? 0.0 : d1 / moment_of_inertia;
            double t4 = other.moment_of_inertia > std::numeric_limits<double>::infinity() ? 0.0 : d2 / other.moment_of_inertia;

            if (t1 + t2 + t3 + t4 == 0) return false;


            auto r_impulse = rvel_fcomp * ( (1.0 + c_elastic) / (t1 + t2 + t3 + t4));

            // max impulse is when pairwise elasticity is 1.0, and the relative velocity is negative
            // of the initial relative velocity, therefore max impuse is twice min impulse

            // interpolate by elasticity to find impulse
            impulse_point impulse_p1{p, r_impulse};
            impulse_point impulse_p2{p, -r_impulse};
            apply_impulse(impulse_p1);
            other.apply_impulse(impulse_p2);
            return true;
        }
        return false;
    }

public:
    body() = delete;
    body(const convex_polygon& _shape, const affine_2d& _pos, double _mass, double _moment_of_inertia, double _elasticity)
    : shape(_shape), position(_pos), mass(_mass), moment_of_inertia(_moment_of_inertia), elasticity(_elasticity),
      velocity{0.0, 0.0}, angular_frequency(0.0) {}

    void handle_collision(body& other) {
        std::vector<std::pair<vec2,vec2>> contact_points;
        for (auto p : shape.getPoints()) {
            auto ptrue = position(p);
            auto pomap = other.position.inverse()(ptrue);
            auto normal = other.position.scale * other.shape.normalAt(pomap);
            if (normal == vec2{0.0, 0.0}) continue;
            contact_points.emplace_back(ptrue, normal);
        }
        for (auto q : other.shape.getPoints()) {
            auto qtrue = other.position(q);
            auto qimap = position.inverse()(qtrue);
            auto normal = -position.scale * shape.normalAt(qimap);
            if (normal == vec2{0.0, 0.0}) continue;
            contact_points.emplace_back(qtrue, normal);
        }
        for (int i = 0; i < 4; ++i) {
            if (!collision_iteration(other, contact_points)) break;
        }
    }
    vec2 velocityAt(vec2 point) const {
        auto p = point - position.shift;
        return velocity + angular_frequency * (mat2{0,-1,1,0} * p);
    }
    void apply_impulse(impulse_point imp) {
        auto adj_origin = position.inverse()(imp.origin);
        auto adj_impulse = position.scale.inverse() * imp.impulse;
        angular_frequency += (1/moment_of_inertia) * adj_impulse.cross(-adj_origin);
        velocity += (1/mass) * imp.impulse;
    }
    void update(double delta) {
        position.scale = mat2::rot_mat(angular_frequency * delta) * position.scale;
        position.shift += velocity * delta;

    }
    const convex_polygon& getShape() const {
        return shape;
    }
    auto getTranslatedPoints() const {
        return std::views::transform(shape.getPoints(), [this](vec2 v) {return position(v);});
    }
    affine_2d getPos() const {
        return position;
    }
    void setVelocity(vec2 v) {
        velocity = v;
    }
    vec2 getVelocity() {
        return velocity;
    }
    double kineticEnergy() const {
        double translational = velocity == vec2::zerovec() ? 0.0 : velocity.dot(velocity) * mass;
        double rotational = angular_frequency == 0.0 ? 0.0 : angular_frequency * angular_frequency * moment_of_inertia;
        return 0.5 * (translational + rotational);
    }
};

}

#endif
