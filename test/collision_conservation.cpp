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

bool almost_equal(double a, double b) {
    return (0.99999999 <= a/b && a/b <= 1.0000001) || std::abs(a - b) <= std::numeric_limits<double>::epsilon() ;
}

std::function<void(contact_point*)> test_kinetic_energy_callback(bool* out_var, double previous_KE) {
    return [=](contact_point* cp) {
        if (cp->body1->getElasticity() * cp->body2->getElasticity() < 1.0) {
            *out_var = cp->body1->kineticEnergy() + cp->body2->kineticEnergy() <= previous_KE;
        }
        else if (cp->body1->getElasticity() * cp->body2->getElasticity() > 1.0) {
            *out_var = cp->body1->kineticEnergy() + cp->body2->kineticEnergy() >= previous_KE;
        }
        else *out_var = almost_equal(cp->body1->kineticEnergy() + cp->body2->kineticEnergy(), previous_KE);
    };
}

bool test_kinetic_energy(std::vector<mxphys::body>& bodies, size_t iterations, double delta) {
    bool test_bool = true;
    for (size_t i = 0; i < iterations; ++i) {
        std::vector<contact_point> contact_points;
        for (auto it = bodies.begin(); it != bodies.end(); ++it) {
            for (auto jt = it + 1; jt != bodies.end(); ++jt) {
                double pke = it->kineticEnergy() + jt->kineticEnergy();
                it->get_contact_points(*jt, contact_points, {test_kinetic_energy_callback(&test_bool, pke)});
            }
        }

        bool unresolved = true;
        for (int i = 0; i < 5 && unresolved; ++i) {
            unresolved = false;
            for (contact_point& cp : contact_points) {
                unresolved |= cp.resolve();
                if (!test_bool) {
                    std::cout << "FAILED KINETIC ENERGY CONSERVATION" << std::endl;
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

int main() {

    constexpr double bb_l = -100.0;
    constexpr double bb_r = 100.0;
    constexpr double bb_u = -100.0;
    constexpr double bb_d = 100.0;

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
            bodies.back().setVelocity(vec2{(px - 0.5) * 7.0, (py - 0.5) * 10.0});
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
        mxphys::affine_2d{mxphys::mat2::identity(), vec2{0.0,bb_u - 1.0}},
        std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity(), 1.0
    );
    bodies.emplace_back(
        xbound,
        mxphys::affine_2d{mxphys::mat2::identity(), vec2{0.0,bb_d + 0.99}},
        std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity(), 1.0
    );
    bodies.emplace_back(
        ybound,
        mxphys::affine_2d{mxphys::mat2::identity(), vec2{bb_l -1.0, 0.0}},
        std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity(), 1.0
    );
    bodies.emplace_back(
        ybound,
        mxphys::affine_2d{mxphys::mat2::identity(), vec2{bb_r + 0.99, 0.0}},
        std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity(), 1.0
    );

    if (!test_kinetic_energy(bodies, 10000, 0.001)) return 1;
    else return 0;
}