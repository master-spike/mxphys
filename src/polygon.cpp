#include "mxphys/polygon.h"

namespace mxphys {
void convex_polygon::convexify() {

    std::sort(points.begin(), points.end(), [](const vec2& lhs, const vec2& rhs) {
        if (lhs.x < rhs.x) return true;
        if (lhs.x == rhs.x) return lhs.y < rhs.y;
        return false;
    });

        vec2 p_bl = *points.begin();
    vec2 p_tr = *points.rbegin();
    if (p_bl.x == p_tr.x) {
        points = {p_bl, p_tr};
        return;
    }
    auto above_line = [p_bl, p_tr](const vec2& p) -> bool {
        double t = (p.x - p_bl.x) / (p_tr.x - p_bl.x);
        return std::lerp(p_bl.y, p_tr.y, t) >= p.y; 
    };
    auto gradient = [](const vec2& p1, const vec2& p2) {
        return (p2.y - p1.y) / (p2.x - p1.x);
    };
    auto part = std::stable_partition(points.begin(), points.end(), above_line);
    auto l_end = part;
    if (std::distance(points.begin(), l_end) > 2) {
        for (auto it = points.begin() + 2; it < l_end; ++it) {
            while (it >= points.begin() + 2) {
                auto const& v1 = *(it-2);
                auto const& v2 = *(it-1);
                auto const& v3 = *it;
                if (gradient(v1, v2) < gradient(v1, v3)) break;
                --it;
                std::rotate(it, it+1, l_end);
                 --l_end;
            }
        }
    }
    auto u_end = points.end();
    if (std::distance(part, u_end) > 2){ 
        for (auto it = part + 2; it < u_end; ++it) {
            while (it >= part + 2) {
                auto const& v1 = *(it-2);
                auto const& v2 = *(it-1);
                auto const& v3 = *it;
                if ((v1.x == v2.x) | ((gradient(v1, v2) > gradient(v1, v3)))) break;
                --it;
                std::rotate(it, it+1, u_end);
                --u_end;
            }
        }
    }
    auto end = l_end;
    for (; part < u_end; ++part, ++end) {
        std::swap(*end, *part);
    }
    std::ranges::reverse(l_end, end);
    points.erase(end, points.end());
    points.shrink_to_fit();
    
}

}
