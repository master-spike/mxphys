#include <vector>
#include <cmath>
#include <limits>
#include <iostream>

#include "mxphys/types.h"
#include "mxphys/polygon.h"
#include "mxphys/body.h"
#include "mxphys/forces.h"

using mxphys::vec2;
using mxphys::affine_2d;
using mxphys::mat2;
using mxphys::convex_polygon;
using mxphys::contact_point;

namespace collision_conservation_test_parameters {

bool cctp_verbose = false;

}

bool almost_equal(double a, double b) {
    return (0.99999999 <= a/b && a/b <= 1.0000001) || std::abs(a - b) <= std::numeric_limits<double>::epsilon() ;
}

std::function<void(contact_point*)> test_kinetic_energy_precallback(double* previous_KE) {
    return [=](contact_point* cp) {
        *previous_KE = cp->body1->kineticEnergy() + cp->body2->kineticEnergy();
    };
}

std::function<void(contact_point*)> test_linear_momentum_precallback(vec2* previous_LM) {
    return [=](contact_point* cp) {
        vec2 lm1 = (cp->body1->getVelocity() == vec2::zerovec()) ? vec2::zerovec() : cp->body1->getMass() * cp->body1->getVelocity();
        vec2 lm2 = (cp->body2->getVelocity() == vec2::zerovec()) ? vec2::zerovec() : cp->body2->getMass() * cp->body2->getVelocity();
        *previous_LM = lm1 + lm2;
    };
}


std::function<void(contact_point*)> test_kinetic_energy_postcallback(bool* out_var, double* previous_KE) {
    using namespace collision_conservation_test_parameters;
    return [=](contact_point* cp) {
        if (cctp_verbose) {
            std::cout << "Collision detected - kinetic energy " << *previous_KE << "->" << cp->body1->kineticEnergy() + cp->body2->kineticEnergy() << std::endl;
        }
        if (cp->body1->getElasticity() * cp->body2->getElasticity() < 1.0) {
            *out_var = cp->body1->kineticEnergy() + cp->body2->kineticEnergy() <= *previous_KE;
        }
        else if (cp->body1->getElasticity() * cp->body2->getElasticity() > 1.0) {
            *out_var = cp->body1->kineticEnergy() + cp->body2->kineticEnergy() >= *previous_KE;
        }
        else *out_var = almost_equal(cp->body1->kineticEnergy() + cp->body2->kineticEnergy(), *previous_KE);
    };
}

std::function<void(contact_point*)> test_linear_momentum_postcallback(bool* out_var, vec2* previous_momentum) {
    using namespace collision_conservation_test_parameters;
    return [=](contact_point* cp) {
        vec2 lm1 = (cp->body1->getVelocity() == vec2::zerovec()) ? vec2::zerovec() : cp->body1->getMass() * cp->body1->getVelocity();
        vec2 lm2 = (cp->body2->getVelocity() == vec2::zerovec()) ? vec2::zerovec() : cp->body2->getMass() * cp->body2->getVelocity();
        if (cctp_verbose) {
            std::cout << "Collision detected - linear momentum " << *previous_momentum << "->" << lm1 + lm2 << std::endl;
        }
        if (cp->body1->getMass() > std::numeric_limits<double>::max() || cp->body2->getMass() > std::numeric_limits<double>::max()) {
            // infinite mass objects are an exceptional case as they can act as infinite sinks / sources of momentum
            *out_var = true;
            return;
        }
        *out_var = almost_equal(lm1.x + lm2.x, previous_momentum->x) && almost_equal(lm1.y + lm2.y, previous_momentum->y);
    };
}

bool test_conservation(std::vector<mxphys::body>& bodies, size_t iterations, double delta) {
    bool test_ke = true;
    bool test_lm = true;
    double pre_ke = 0.0;
    vec2 pre_lm = vec2::zerovec();
    for (size_t i = 0; i < iterations; ++i) {
        std::vector<contact_point> contact_points;
        for (auto it = bodies.begin(); it != bodies.end(); ++it) {
            for (auto jt = it + 1; jt != bodies.end(); ++jt) {
                it->get_contact_points(*jt, contact_points, {
                    test_kinetic_energy_precallback(&pre_ke),
                    test_linear_momentum_precallback(&pre_lm)
                }, {
                    test_kinetic_energy_postcallback(&test_ke, &pre_ke),
                    test_linear_momentum_postcallback(&test_lm, &pre_lm)
                });
            }
        }

        bool unresolved = true;
        for (int i = 0; i < 5 && unresolved; ++i) {
            unresolved = false;
            for (contact_point& cp : contact_points) {
                bool coll = cp.resolve();
                unresolved |= coll;
                if (!(test_ke && test_lm)) {
                    if (!test_ke)
                        std::cout << "FAILED KINETIC ENERGY CONSERVATION" << std::endl;
                    if (!test_lm)
                        std::cout << "FAILED LINEAR MOMENTUM CONSERVATION" << std::endl;
                    std::cout << 
                    "masses: " << cp.body1->getMass() << "," << cp.body2->getMass() << std::endl <<
                    "surface normal: " << cp.normal << std::endl <<
                    "contact position: " << cp.position << std::endl <<
                    "" << std::endl;
                    return false;
                }
            }
        }

        for (mxphys::body& b : bodies) {
            b.update(delta);
        }

    }

    return true;
}

int main(int argc, char** argv) {

    constexpr double bb_l = -100.0;
    constexpr double bb_r = 100.0;
    constexpr double bb_u = -100.0;
    constexpr double bb_d = 100.0;

    for (int i = 0; i < argc; ++i) {
        std::cout << argv[i] << std::endl;
        if (std::string_view(argv[i]) == std::string_view("--verbose")) {
            std::cout << "Running in verbose mode" << std::endl;
            collision_conservation_test_parameters::cctp_verbose = true;
        }
    }

    auto regular_poly = [](int n, double r, mxphys::affine_2d pos){
        std::vector<vec2> points;
        vec2 v{r, 0};
        auto rmat = mxphys::mat2::rot_mat(M_PI * 2.0 / static_cast<double>(n));
        for (int i = 0; i < n; ++i){
            points.emplace_back(v);
            v = rmat * v;
        }
        return mxphys::body(
            mxphys::convex_polygon(points.begin(), points.end()),
            pos, M_PI*r*r, M_PI*r*r*r*r / 2.0, 1.0
        );
    };

    std::vector<mxphys::body> bodies;

    // create 81 sample polygons
    for (int i = 0; i < 9; ++i) {
        for (int j = 0; j < 9; ++j) {
            auto px = static_cast<double>(i + 1) / 10.0;
            auto py = static_cast<double>(j + 1) / 10.0; 
            bodies.emplace_back(regular_poly( (i + j) % 5 + 3, ((i + j) % 6 + 4) / 4.0,
                mxphys::affine_2d{
                    mxphys::mat2::identity(),
                    vec2{std::lerp(bb_l, bb_r, px), std::lerp(bb_u, bb_d, py)}
                }
            ));
            bodies.back().setVelocity(vec2{(0.5 - px) * 10.0, (0.5 - py) * 10.0});
        }
    }

    mxphys::convex_polygon xbound{
        vec2{bb_l, -1.0},
        vec2{bb_r, -1.0},
        vec2{bb_r, +1.0},
        vec2{bb_l, +1.0}
    };
    mxphys::convex_polygon ybound{
        vec2{-1.0, bb_u},
        vec2{+1.0, bb_u},
        vec2{+1.0, bb_d},
        vec2{-1.0, bb_d}
    };

    bodies.emplace_back(
        xbound,
        mxphys::affine_2d{mxphys::mat2::identity(), vec2{0.0, bb_u - 1.0}},
        std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity(), 1.0
    );
    bodies.emplace_back(
        xbound,
        mxphys::affine_2d{mxphys::mat2::identity(), vec2{0.0, bb_d + 0.99}},
        std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity(), 1.0
    );
    bodies.emplace_back(
        ybound,
        mxphys::affine_2d{mxphys::mat2::identity(), vec2{bb_l - 1.0, 0.0}},
        std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity(), 1.0
    );
    bodies.emplace_back(
        ybound,
        mxphys::affine_2d{mxphys::mat2::identity(), vec2{bb_r + 0.99, 0.0}},
        std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity(), 1.0
    );

    if (!test_conservation(bodies, 10000, 0.01)) return 1;
    else return 0;
}