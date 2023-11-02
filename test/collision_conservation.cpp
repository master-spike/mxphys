#include <vector>
#include <cmath>
#include <limits>
#include <iostream>

#include "mxphys/types.h"
#include "mxphys/polygon.h"
#include "mxphys/body.h"
#include "mxphys/forces.h"
#include "mxphys/event/eventmanager.h"
#include "mxphys/event/collision_event.h"
#include "mxphys/event/event_data.h"

using mxphys::vec2;
using mxphys::affine_2d;
using mxphys::mat2;
using mxphys::convex_polygon;
using mxphys::contact_point;

namespace collision_conservation_test_parameters {

bool cctp_verbose = false;

}

namespace collision_conservation_test_results {
    std::vector<mxphys::event::collision_event> failed_tests;
}

bool almost_equal(double a, double b) {
    return (0.99999999 <= a/b && a/b <= 1.0000001) || std::abs(a - b) <= std::numeric_limits<double>::epsilon() ;
}


void test_collision(std::unique_ptr<mxphys::event::event_data> const& edata) {
    mxphys::event::collision_event const& coll = static_cast<mxphys::event::collision_event const&>(*edata);
    using namespace collision_conservation_test_parameters;
    using namespace collision_conservation_test_results;
    
    auto calc_ke = [](std::pair<double, double> const& masses, std::pair<double, double> const& moments_of_inertia,
                      std::pair<vec2,vec2> const& velocities, std::pair<double, double> const& angular){
        double ke1l = (velocities.first == vec2{0.0, 0.0}) ? 0.0 : masses.first * velocities.first.dot(velocities.first);
        double ke2l = (velocities.second == vec2{0.0, 0.0}) ? 0.0 : masses.second * velocities.second.dot(velocities.second);
        double ke1a = (angular.first == 0.0) ? 0.0 : moments_of_inertia.first * angular.first * angular.first;
        double ke2a = (angular.second == 0.0) ? 0.0 : moments_of_inertia.second * angular.second * angular.second;
        return 0.5 * (ke1a + ke2a + ke1l + ke2l);
    };

    auto calc_momentum_l = [](std::pair<double, double> const& masses, std::pair<vec2,vec2> const& velocities) {
        if (std::max(masses.first, masses.second) > std::numeric_limits<double>::max())
            return vec2::zerovec();
        return masses.first * velocities.first + masses.second * velocities.second;
    };

    double ke_init = calc_ke(coll.masses, coll.moments_of_inertia, coll.vs_init, coll.afs_init);
    double ke_final = calc_ke(coll.masses, coll.moments_of_inertia, coll.vs_final, coll.afs_final);
    vec2 lm_init = calc_momentum_l(coll.masses, coll.vs_init);
    vec2 lm_final = calc_momentum_l(coll.masses, coll.vs_final);

    if (cctp_verbose) {
        std::cout << "Collision detected: IDs " << coll.ids.first << "," << coll.ids.second << std::endl
                  << "kinetic energy " << ke_init << "->" << ke_final << std::endl
                  << "linear momentum " << lm_init << "->" << lm_final << std::endl;
    }
    if ((coll.restitution < 1.0 && ke_init < ke_final) || 
        (coll.restitution > 1.0 && ke_init > ke_final) || 
        (coll.restitution == 1.0 && !almost_equal(ke_init, ke_final))) {
        failed_tests.emplace_back(coll);
    }
    else if (!almost_equal(lm_init.x, lm_final.x) ||
             !almost_equal(lm_init.y, lm_final.y)) {
        failed_tests.emplace_back(coll);
    }
}

bool test_conservation(std::vector<mxphys::body>& bodies, size_t iterations, double delta) {

    mxphys::event::eventmanager e_manager;
    e_manager.register_listener(mxphys::event::event_type::COLLISION, test_collision);
    for (size_t i = 0; i < iterations; ++i) {
        std::vector<contact_point> contact_points;
        for (auto it = bodies.begin(); it != bodies.end(); ++it) {
            for (auto jt = it + 1; jt != bodies.end(); ++jt) {
                it->get_contact_points(*jt, contact_points);
            }
        }

        bool unresolved = true;
        for (int i = 0; i < 5 && unresolved; ++i) {
            unresolved = false;
            for (contact_point& cp : contact_points) {
                unresolved |= cp.resolve(e_manager);
            }
        }

        e_manager.emit_all();

        for (mxphys::body& b : bodies) {
            b.update(delta);
        }

    }

    if (!collision_conservation_test_results::failed_tests.empty()) {
        std::cout << "FAILED TESTS" << std::endl;
        return false;
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