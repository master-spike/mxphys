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

#include "types.h"

namespace mxphys {

template<typename T>
class bounding_volume_heirarchy {
    struct bvh_node {
        bounding_box m_BB;
        std::optional<T> m_Value;
        std::optional<
            std::pair<std::unique_ptr<bvh_node>, std::unique_ptr<bvh_node>>
        > m_Children;
        bvh_node(std::vector<std::pair<T, bounding_box>>::iterator ts_l, std::vector<std::pair<T, bounding_box>>::iterator ts_r) : m_BB{vec2::zerovec(), vec2::zerovec()}, m_Value{std::nullopt}, m_Children{std::nullopt}{
            if (std::distance(ts_l, ts_r) == 0) {
                return;
            }
            m_BB = ts_l->second;
            if (std::distance(ts_l, ts_r) == 1) {
                m_Value = ts_l->first;
                return;
            }
            m_BB = std::accumulate(ts_l, ts_r, m_BB, [](const bounding_box& lhs, const std::pair<T, bounding_box>& rhs) {
                return lhs.bb_union(rhs.second);
            });

            // we don't need to divide by 2 as we are simply using it to compare
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
            m_Value = midpoint->first;
            m_Children = std::make_pair(
                std::make_unique<bvh_node>(ts_l, midpoint),
                std::make_unique<bvh_node>(midpoint+1, ts_r)
            );
        }
        std::size_t traverse(const bounding_box& bb, std::function<void(const T&)> func) const {
            if (!bb.intersects(m_BB)) return 0;
            std::size_t c = 0;
            if (m_Value.has_value()) {
                func(m_Value.value());
                c = 1;
            }
            if (m_Children.has_value()) {
                c += m_Children.value().first->traverse(bb,func);
                c += m_Children.value().second->traverse(bb,func);
            }
            return c;
        }

        void getBoundingBoxes(auto out) const {
            *out = m_BB;
            ++out;
            if (m_Children.has_value()) {
                m_Children.value().first->getBoundingBoxes(out);
                m_Children.value().second->getBoundingBoxes(out);
            }
        }

        int countElems(std::unordered_set<T>& seen) const {
            int c = 0;
            if (m_Value.has_value()) {
                ++c;
                seen.emplace(m_Value().value());
            }
            if (m_Children.has_value()) {
                c += m_Children.value().first->countElems(seen);
                c += m_Children.value().first->countElems(seen);
            }
            return c;
        }
    };

    std::optional<bvh_node> top_node;

public:
    bounding_volume_heirarchy()  = delete;
    bounding_volume_heirarchy(auto it_elems, auto it_end,
                              auto t_from_u,
                              auto bb_from_u) : top_node(std::nullopt)
    {
        static_assert(std::is_trivial_v<T>, "T is not a trivial type");
        std::vector<std::pair<T, bounding_box>> ts;
        std::transform(it_elems, it_end, std::back_inserter(ts), [&t_from_u, &bb_from_u](auto const& u) {
            return std::make_pair(t_from_u(u), bb_from_u(u));
        });
        top_node = bvh_node{ts.begin(), ts.end()};
    }
    
    std::size_t for_each_possible_colliding(const bounding_box& bb, std::function<void(const T&)> func) {
        if (top_node.has_value()) return top_node->traverse(bb, func);
        return 0;
    }

    std::vector<bounding_box> getBoundingBoxes() const {
        std::vector<bounding_box> boxes;
        if (top_node.has_value()) top_node->getBoundingBoxes(std::back_inserter(boxes));
        return boxes;
    }
};

}

#endif
