#ifndef MXPHYS_BODY_H
#define MXPHYS_BODY_H

#include <utility>
#include <vector>
#include <functional>
#include <ranges>
#include <optional>

#include "types.h"
#include "polygon.h"

namespace mxphys {

class body {
private:

    static uint64_t id_next;

    convex_polygon shape;
    // position of center of mass
    affine_2d position;


    double mass;
    double moment_of_inertia;
    double elasticity;

    vec2 velocity;
    double angular_frequency;

    const uint64_t id;

    std::optional<std::pair<double,vec2>> sat_contact_normal(const body& other) const;

public:
    body() = delete;
    body(const convex_polygon& _shape, const affine_2d& _pos, double _mass, double _moment_of_inertia, double _elasticity)
      : shape(_shape), position(_pos), mass(_mass), moment_of_inertia(_moment_of_inertia), elasticity(_elasticity),
        velocity{0.0, 0.0}, angular_frequency(0.0), id(id_next) {
        ++id_next;
    }

    void get_contact_points(body& other, std::vector<contact_point>& out_contact_points);
    vec2 velocityAt(vec2 point) const {
        auto p = point - position.shift;
        return velocity + angular_frequency * (mat2{0,-1,1,0} * p);
    }
    void apply_impulse(impulse_point imp);

    void update(double delta) {
        position.scale = mat2::rot_mat(angular_frequency * delta) * position.scale;
        position.shift += velocity * delta;

    }
    const convex_polygon& getShape() const {
        return shape;
    }

    void setShape(const convex_polygon& new_shape) {
        shape = new_shape;
    }

    auto getTranslatedPoints() const {
        return std::views::transform(shape.getPoints(), position);
    }

    affine_2d getPos() const {
        return position;
    }
    void setPos(affine_2d pos) {
        position = pos;
    }

    void setVelocity(vec2 v) {
        velocity = v;
    }
    vec2 getVelocity() const {
        return velocity;
    }

    double getAngularFreq() const {
        return angular_frequency;
    }
    void setAngularFreq(double af) {
        angular_frequency = af;
    }
    
    double getMass() const {
        return mass;
    }
    void setMass(double m) {
        mass = m;
    }

    double getMomentOfInertia() const {
        return moment_of_inertia;
    }
    void setMomentOfInertia(double i) {
        moment_of_inertia = i;
    }

    double getElasticity() const {
        return elasticity;
    }
    void setElasticity(double e) {
        elasticity = e;
    }

    uint64_t getID() const {
        return id;
    }

    bounding_box getBoundingBox() const;

    double kineticEnergy() const {
        double translational = velocity == vec2::zerovec() ? 0.0 : velocity.dot(velocity) * mass;
        double rotational = angular_frequency == 0.0 ? 0.0 : angular_frequency * angular_frequency * moment_of_inertia;
        return 0.5 * (translational + rotational);
    }
};

}

#endif
