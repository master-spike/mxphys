#ifndef MXPHYS_POLYGON_H
#define MXPHYS_POLYGON_H

#include <vector>
#include <span>
#include <algorithm>
#include <optional>
#include <utility>
#include <functional>
#include "types.h"

namespace mxphys {

class convex_polygon {
private:
    // counter clockwise
    std::vector<vec2> points;

    // turns points vector into convex polygon containing all points.
    // called by constructor once. Points vector is immutable after construction.
    void convexify();
    
public:
    // counter clockwise
    std::span<const vec2> getPoints() const {
        return std::span{points.cbegin(), points.cend()};
    }
    convex_polygon() = delete;
    template<class InputIt>
    convex_polygon(InputIt first, InputIt last) : points(first, last) {
        if (points.size() <= 3) return;
        convexify();
    }

    convex_polygon(std::initializer_list<vec2> _il_points) : points(_il_points) {
        if (points.size() <= 3) return;
        convexify();
    }

    /*
    std::optional<std::pair<convex_polygon, convex_polygon>> split(vec2 l1, vec2 l2) const {
        if (l1 == l2 || points.size() <= 1) return std::nullopt;
        auto above_line = [l1,l2](const vec2& p) {
            if (l1.x == l2.x) return p.x > l1.x;
            double t = (p.x - l1.x) / (l2.x - l1.x);
            return std::lerp(l1.y, l2.y, t) > p.y; 
        };
        
        std::vector<vec2> new_points(points);
        auto it = above_line(new_points[0]) ? std::find_if(new_points.begin(), new_points.end(), std::not1(above_line))
                                            : std::find_if(new_points.begin(), new_points.end(), above_line);
        
        if (it == new_points.end()) return std::nullopt; // line doesn't intersect polygon
        std::rotate(new_points.begin(), it, new_points.end());
        it = new_points.begin();
        it = above_line(*it) ? std::find_if(new_points.begin(), new_points.end(), std::not1(above_line))
                             : std::find_if(new_points.begin(), new_points.end(), above_line);
        
    }
    */
    vec2 normalAt(vec2 v) const {
        double min_dist = +INFINITY;
        vec2 best_normal{0,0};
        for (auto it = points.cbegin(); it != points.cend(); ++it) {
            auto const& p1 = *it;
            auto const& p2 = it + 1 == points.cend() ? *points.cbegin() : *(it+1);
            vec2 v_ = v - p1;
            vec2 w_ = p2 - p1;
            // alignw_ is (v,0) where v > 0;
            mat2 align{w_.x, w_.y, -w_.y, w_.x};
            if ((align * v_).y < 0) return vec2{0,0};
            if ((align * v_).y < min_dist) {
                min_dist = (align * v_).y;
                best_normal = vec2{w_.y, -w_.x};
            }
        } 
        return best_normal.normalized();
    }
};

}

#endif
