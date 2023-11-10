#ifndef MXPHYS_ROPE_H
#define MXPHYS_ROPE_H

#include <functional>
#include <utility>

#include "types.h"
#include "forces.h"
#include "body.h"

namespace mxphys {

class elastic_rope {

    uint64_t body1_id;
    uint64_t body2_id;

    vec2 body1_attachment;
    vec2 body2_attachment;

    double resting_length;
    double f_tension;
    double f_compression;
    double shock_coeff;

public:

    elastic_rope(uint64_t _b1_id, uint64_t _b2_id, vec2 _b1_attach, vec2 _b2_attach,
        double _resting_len, double _f_tension, double _f_compression, double _shock_coeff)
        : body1_id(_b1_id), body2_id(_b2_id), body1_attachment(_b1_attach),
        body2_attachment(_b2_attach), resting_length(_resting_len),
        f_tension(_f_tension), f_compression(_f_compression), shock_coeff(_shock_coeff) {}

    void apply_impulses(double delta, std::function<mxphys::body*(uint64_t)> get_body_from_id) const {
        auto body1 = get_body_from_id(body1_id);
        auto body2 = get_body_from_id(body2_id);
        if (!body1 || !body2) return;
        vec2 p1 = body1->getPos()(body1_attachment);
        vec2 p2 = body2->getPos()(body2_attachment);
        if (p1 == p2) return;
        double extension = std::sqrt((p2-p1).dot(p2-p1)) - resting_length;
        if (extension > 0.0) extension *= f_tension;
        else if (extension < 0.0) extension *= f_compression;
        if (extension == 0) return;
        vec2 rvel = (p2 - p1).normalized().dot(
            body1->velocityAt(p1) - body2->velocityAt(p2)
        ) * (p2 - p1).normalized();
        vec2 impulse = ((p2 - p1).normalized() * extension - shock_coeff * rvel) * delta;
        body1->apply_impulse(impulse_point{p1, impulse});
        body2->apply_impulse(impulse_point{p2, -impulse});
    }

    std::pair<vec2, vec2> getRopePoints(std::function<mxphys::body*(uint64_t)> get_body_from_id) const {
        auto body1 = get_body_from_id(body1_id);
        auto body2 = get_body_from_id(body2_id);
        if (!body1 || !body2) {
            return std::make_pair(vec2::zerovec(), vec2::zerovec());
        }
        return std::make_pair(body1->getPos()(body1_attachment), body2->getPos()(body2_attachment));
    }

    double getRestingLength() const {
        return resting_length;
    }
    double getFTension() const {
        return f_tension;
    }
    double getFCompression() const {
        return f_compression;
    }

    double calculateScalarTension(std::function<const mxphys::body*(uint64_t)> get_body_from_id) const {
        auto body1 = get_body_from_id(body1_id);
        auto body2 = get_body_from_id(body2_id);
       if (!body1 || !body2) return 0;
        vec2 p1 = body1->getPos()(body1_attachment);
        vec2 p2 = body2->getPos()(body2_attachment);
        double extension = std::sqrt((p2-p1).dot(p2-p1)) - resting_length;
        if (extension > 0.0) extension *= f_tension;
        else if (extension < 0.0) extension *= f_compression;
        return extension;
    }

};

}


#endif
