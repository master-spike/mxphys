#ifndef MXPHYS_BODY_H
#define MXPHYS_BODY_H

#include <utility>
#include <vector>
#include <functional>
#include <ranges>

#include "types.h"
#include "polygon.h"

namespace mxphys {

class body {
private:

    convex_polygon shape;
    // position of center of mass
    affine_2d position;


    double mass;
    double moment_of_inertia;
    double elasticity;

    vec2 velocity;
    double angular_frequency;

public:
    body() = delete;
    body(const convex_polygon& _shape, const affine_2d& _pos, double _mass, double _moment_of_inertia, double _elasticity)
    : shape(_shape), position(_pos), mass(_mass), moment_of_inertia(_moment_of_inertia), elasticity(_elasticity),
      velocity{0.0, 0.0}, angular_frequency(0.0) {}

    void get_contact_points(body& other, std::vector<contact_point>& out_contact_points);
    void get_contact_points(body & other, std::vector<contact_point>& out_contact_points,
                            std::initializer_list<std::function<void(contact_point*)>> callbacks);
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
    auto getTranslatedPoints() const {
        return std::views::transform(shape.getPoints(), [this](vec2 v) {return position(v);});
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

    double kineticEnergy() const {
        double translational = velocity == vec2::zerovec() ? 0.0 : velocity.dot(velocity) * mass;
        double rotational = angular_frequency == 0.0 ? 0.0 : angular_frequency * angular_frequency * moment_of_inertia;
        return 0.5 * (translational + rotational);
    }
};

}

#endif
