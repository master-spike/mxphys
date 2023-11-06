#ifndef MXPHYS_BV_HEIRARCHY_H
#define MXPHYS_BV_HEIRARCHY_H

#include <memory>
#include <utility>
#include <vector>
#include <functional>
#include <optional>
#include <type_traits>
#include <cstdint>
#include <algorithm>
#include <numeric>
#include <unordered_set>
#include <iterator>

#include "types.h"

namespace mxphys {

template<typename T>
class bounding_volume_heirarchy {


    struct bvh_node {
        
        bounding_box m_BB;
        
        std::optional<T> m_Value;
        std::optional<std::size_t> m_left;
        std::optional<std::size_t> m_right;

        bvh_node(bounding_box _bb, std::optional<T> _val, std::optional<std::size_t> _l, std::optional<std::size_t> _r)
        : m_BB(_bb), m_Value(_val), m_left(_l), m_right(_r) {

        }

    };

    std::vector<bvh_node> nodes;

    std::optional<std::size_t> constructBVHNode(std::vector<std::pair<T, bounding_box>>::iterator ts_l,
                                 std::vector<std::pair<T, bounding_box>>::iterator ts_r) {
        if (std::distance(ts_l, ts_r) == 0) {
            return std::nullopt;
        }
        std::size_t s = nodes.size();
        if (std::distance(ts_l, ts_r) == 1) {
            nodes.emplace_back(ts_l->second, ts_l->first, std::nullopt, std::nullopt);
            return s;
        }

        bounding_box node_bb = ts_l->second;
        node_bb = std::accumulate(ts_l, ts_r, node_bb, [](const bounding_box& lhs, const std::pair<T, bounding_box>& rhs) {
            return lhs.bb_union(rhs.second);
        });

        auto bb_centroidx2 = [](const bounding_box& bb) -> vec2 {
            return vec2{bb.bottom_left.x + bb.top_right.x, bb.bottom_left.y + bb.top_right.y};
        };
        auto cmp_x = [&bb_centroidx2](const auto& lhs, const auto& rhs) {
            return bb_centroidx2(lhs.second).x < bb_centroidx2(rhs.second).x;
        };
        auto cmp_y = [&bb_centroidx2](const auto& lhs, const auto& rhs) {
            return bb_centroidx2(lhs.second).y < bb_centroidx2(rhs.second).y;
        };

        auto midpoint = ts_l;
        std::advance(midpoint, std::distance(ts_l, ts_r) / 2);


        std::nth_element(ts_l, midpoint, ts_r, cmp_x);
        double x_part_line = midpoint->second.bottom_left.x;
        std::size_t x_overlap = std::count_if(ts_l, ts_r, [x_part_line](const std::pair<T, bounding_box>& pair) {
            return pair.second.top_right.x >= x_part_line;
        }) + std::count_if(ts_l, ts_r, [x_part_line](const std::pair<T, bounding_box>& pair) {
            return pair.second.bottom_left.x <= x_part_line;
        }) - std::distance(ts_l, ts_r);

        // find overlap along y-axis
        std::nth_element(ts_l, midpoint, ts_r, cmp_y);
        double y_part_line = midpoint->second.bottom_left.y;
        std::size_t y_overlap = std::count_if(ts_l, ts_r, [y_part_line](const std::pair<T, bounding_box>& pair) {
            return pair.second.top_right.y >= y_part_line;
        }) + std::count_if(ts_l, ts_r, [y_part_line](const std::pair<T, bounding_box>& pair) {
            return pair.second.bottom_left.y <= y_part_line;
        }) - std::distance(ts_l, ts_r);

        if (x_overlap < y_overlap) {
            std::nth_element(ts_l, midpoint, ts_r, cmp_x);
        }

        nodes.emplace_back(node_bb, std::nullopt, 0, 0);
        nodes[s].m_left = constructBVHNode(ts_l, midpoint);
        nodes[s].m_right = constructBVHNode(midpoint, ts_r);

        return s;
    }

    std::size_t traverse(const bounding_box& bb, std::function<void(const T&)> func, std::size_t idx) const {
        if (!bb.intersects(nodes[idx].m_BB)) return 0;
        std::size_t c = 0;
        if (nodes[idx].m_Value.has_value()) {
            func(nodes[idx].m_Value.value());
            c = 1;
        }
        if (nodes[idx].m_left.has_value()) {
            c += traverse(bb,func, nodes[idx].m_left.value());
        }
        if (nodes[idx].m_right.has_value()) {
            c += traverse(bb,func, nodes[idx].m_right.value());
        }
        return c;
    }

    void get_bounding_boxes(std::size_t idx, auto output_it) const {
        *output_it = nodes[idx].m_BB;
        ++output_it;
        if (nodes[idx].m_left.has_value()) {
            get_bounding_boxes(nodes[idx].m_left.value(), output_it);
        }
        if (nodes[idx].m_right.has_value()) {
            get_bounding_boxes(nodes[idx].m_right.value(), output_it);
        }
    }

public:
    bounding_volume_heirarchy()  = delete;

    template<class InputIt>
    bounding_volume_heirarchy(InputIt it_elems, InputIt it_end,
                              auto t_from_u,
                              auto bb_from_u)
    {
        static_assert(std::is_trivial_v<T>, "T is not a trivial type");
        std::vector<std::pair<T, bounding_box>> ts;
        std::transform(it_elems, it_end, std::back_inserter(ts), [&t_from_u, &bb_from_u](auto const& u) {
            return std::make_pair(t_from_u(u), bb_from_u(u));
        });
        nodes.reserve(ts.size() * 2);
        constructBVHNode(ts.begin(), ts.end());
    }
    
    std::size_t for_each_possible_colliding(const bounding_box& bb, std::function<void(const T&)> func) const {
        if (!nodes.empty()) return traverse(bb, func, 0);
        return 0;
    }

    std::vector<bounding_box> getBoundingBoxes() const {
        std::vector<bounding_box> boxes;
        if (!nodes.empty()) get_bounding_boxes(0, std::back_inserter(boxes));
        return boxes;
    }
};

}

#endif
