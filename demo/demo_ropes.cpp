#define _USE_MATH_DEFINES
#include <iostream>

#include <random>
#include <SDL.h>
#include <chrono>
#include <unordered_map>
#include <vector>

#include "mxphys/polygon.h"
#include "mxphys/body.h"
#include "mxphys/forces.h"
#include "mxphys/bv_heirarchy.h"
#include "mxphys/rope.h"

int main(int argc, char** argv) {

    // silence silly compiler warnings about these
    (void) argc;
    (void) argv;

    if (SDL_Init(SDL_INIT_EVERYTHING) != 0) {
        
        return 0;
    }

    constexpr int w_height = 960;
    constexpr int w_width = 960;
    SDL_Window* window = SDL_CreateWindow("mxphys demo", 10, 10, w_height, w_width, SDL_WINDOW_SHOWN);

    SDL_Renderer* renderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_ACCELERATED);

    using mxphys::vec2;

    constexpr double bb_l = -15.0;
    constexpr double bb_r = 15.0;
    constexpr double bb_u = -15.0;
    constexpr double bb_d = 15.0;

    auto world_to_scr = [bb_r, bb_l, bb_u, bb_d](mxphys::vec2 v) {
        double x_frac = (v.x - bb_l) / (bb_r - bb_l);
        int x = static_cast<int>(static_cast<double>(w_width) * x_frac);
        double y_frac = (v.y - bb_u) / (bb_d - bb_u);
        int y = static_cast<int>(static_cast<double>(w_height) * y_frac);
        return std::make_pair(x,y);
    };

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
    std::vector<mxphys::elastic_rope> ropes;

    // create 36 sample polygons
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 6; ++j) {
            auto px = static_cast<double>(i + 1) / 7;
            auto py = static_cast<double>(j + 1) / 7;
            bodies.emplace_back(regular_poly( (i + j) % 5 + 3, ((i + j) % 6 + 4) / 5.0,
                mxphys::affine_2d{
                    mxphys::mat2::identity(),
                    vec2{std::lerp(bb_l, bb_r, px), std::lerp(bb_u, bb_d, py)}
                }
            ));
            bodies.back().setVelocity(vec2{(px - 0.5) * 7.0, (py - 0.5) * 7.0});
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


    // infinite mass rectangles at the screen boundaries
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


    auto draw_body = [renderer, &world_to_scr](const mxphys::body& b) {
        auto view = b.getTranslatedPoints();
        auto itl = view.begin();
        auto itr = view.begin() + 1;
        while (itl != view.end()) {
            auto [x1, y1] = world_to_scr(*itl);
            if (itr == view.end()) itr = view.begin();
            auto [x2, y2] = world_to_scr(*itr);
            SDL_RenderDrawLine(renderer, x1, y1, x2, y2);
            ++itr; ++itl;
        }
    };



    auto draw_bounding_box = [renderer, &world_to_scr](const mxphys::bounding_box& bb) {
        auto [x1, y1] = world_to_scr(bb.bottom_left);
        auto [x2, y2] = world_to_scr(bb.top_right);
        SDL_RenderDrawLine(renderer, x1, y1, x2, y1);
        SDL_RenderDrawLine(renderer, x1, y1, x1, y2);
        SDL_RenderDrawLine(renderer, x2, y1, x2, y2);
        SDL_RenderDrawLine(renderer, x1, y2, x2, y2);
    };


    std::unordered_map<uint64_t, std::vector<mxphys::body>::iterator> id_to_body;
    for (auto it = bodies.begin(); it != bodies.end(); ++it) {
        id_to_body.emplace(it->getID(), it);
    }

    auto get_body = [&id_to_body](uint64_t id) -> mxphys::body* {
        auto it = id_to_body.find(id);
        if (it == id_to_body.end()) return nullptr;
        return &(*(it->second));
    };

    auto draw_rope = [renderer, &get_body, &world_to_scr](const mxphys::elastic_rope& rope) -> void {
        auto r_points = rope.getRopePoints(get_body);
        auto [x1, y1] = world_to_scr(r_points.first);
        auto [x2, y2] = world_to_scr(r_points.second);
        SDL_RenderDrawLine(renderer, x1, y1, x2, y2);
    };

    // create 10 sample ropes
    for (uint64_t i = 0; i < 18; ++i) {
        auto& body1 = bodies[i];
        auto& body2 = bodies[i + 18];
        vec2 diff = body1.getPos().shift - body2.getPos().shift;
        ropes.emplace_back(body1.getID(), body2.getID(), vec2::zerovec(), body2.getShape().getPoints().back(),std::sqrt(diff.dot(diff)), 50.0, 0.0, 1.0);
    }

    bool close = false;
    std::chrono::time_point t0 = std::chrono::steady_clock::now();
    std::size_t c = 0;
    while (!close) {
        std::chrono::time_point t1 = std::chrono::steady_clock::now();
        double delta = static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count());
        delta /= 1000000.0;
        t0 = std::chrono::steady_clock::now();
        std::string new_window_title = "mxphys demo - " + std::to_string(delta * 1000.0) + "ms | contact tests:" + std::to_string(c);
        delta = std::max(delta, 0.016);
        SDL_SetWindowTitle(window, new_window_title.c_str());

        mxphys::bounding_volume_heirarchy<uint64_t> bvh_by_id(
            bodies.begin(), bodies.end(),
            [](const mxphys::body & b) { return b.getID(); },
            [](const mxphys::body & b) { return b.getBoundingBox(); }
        );

        for (mxphys::elastic_rope& rope : ropes) {
            rope.apply_impulses(delta, get_body);
        }

        std::vector<mxphys::contact_point> contacts;
        c = 0;
        for (auto it = bodies.begin(); it < bodies.end(); ++it) {
            c += bvh_by_id.for_each_possible_colliding(
                it->getBoundingBox(),
                [&id_to_body, &contacts, &it](uint64_t other_id) {
                    if (it->getID() >= other_id) return;
                    auto jt = id_to_body.find(other_id);
                    if (jt == id_to_body.end()) return;
                    if (jt->second->getID() == other_id) it->get_contact_points(*(jt->second), contacts);
                }
            );
        }
        bool contacts_unresolved = true;
        
        int max_iters = 5;
        for (int i = 0; contacts_unresolved && i < max_iters; ++i) {
            contacts_unresolved = false;
            for (auto& cp : contacts) {
                contacts_unresolved |= cp.resolve();
            }
        }

        for (mxphys::body& b : bodies) {
            b.update(delta);
        }
        
        SDL_Event event;
        while(SDL_PollEvent(&event)) {
            if (event.type == SDL_QUIT) close = true;
        }

        if (delta < 0.001) SDL_Delay(1);

        SDL_SetRenderDrawColor(renderer, 0,0,0,255);
        SDL_RenderClear(renderer);

        SDL_SetRenderDrawColor(renderer, 0, 0, 180, 127);
        for (auto& bb : bvh_by_id.getBoundingBoxes()) {
            draw_bounding_box(bb);
        }

        SDL_SetRenderDrawColor(renderer, 180,180,0,255);
        for (const mxphys::elastic_rope& r : ropes) {
            double tension = r.calculateScalarTension(get_body);
            if (tension < 0.0) {
                SDL_SetRenderDrawColor(renderer, 180,255,0,255);
            }
            else {
                SDL_SetRenderDrawColor(renderer, 255,180,0,255);
            }
            draw_rope(r);
            SDL_SetRenderDrawColor(renderer, 180,180,0,255);
        }
        


        SDL_SetRenderDrawColor(renderer, 255, 255, 255, 255);
        for (const mxphys::body& b : bodies) {
            draw_body(b);
        }

        SDL_SetRenderDrawColor(renderer, 0,0,0,255);

        SDL_RenderPresent(renderer);
    }

    
    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);

    return 0;

}