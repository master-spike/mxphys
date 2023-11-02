#ifndef MXPHYS_TYPES_H
#define MXPHYS_TYPES_H

#include <ostream>
#include <utility>
#include <cmath>

namespace mxphys {

class body;
class convex_polygon;
struct contact_point;
struct impulse_point;
namespace event {
    class eventmanager;
}

struct vec2 {
    double x;
    double y;
    
    vec2() = delete;
    constexpr vec2(double _x, double _y) : x(_x), y(_y) {}

    constexpr vec2& operator+=(const vec2& o) {
        x += o.x; y += o.y;
        return *this;
    }
    constexpr friend vec2 operator+(vec2 v, const vec2& o) {
        v += o;
        return v;
    }
    constexpr vec2& operator-=(const vec2& o) {
        x -= o.x, y -= o.y;
        return *this;
    }
    constexpr friend vec2 operator-(vec2 v, const vec2& o) {
        v -= o;
        return v;
    }
    constexpr double dot(const vec2& o) const {
        return x * o.x + y * o.y;
    }
    constexpr double cross(const vec2& o) const {
        return x * o.y - y * o.x;
    }

    constexpr vec2 normalized() {
        return 1.0 / std::sqrt(this->dot(*this)) * *this;
    }

    constexpr vec2& operator*=(double f) {
        x *= f;
        y *= f;
        return *this;
    }

    constexpr friend vec2 operator*(vec2 v, double f) {
        v *= f;
        return v;
    }

    constexpr friend vec2 operator*(double f, const vec2& v) {
        return v * f;
    }

    constexpr friend vec2 operator-(vec2 v) {
        v.x = -v.x;
        v.y = -v.y;
        return v;
    }

    constexpr static vec2 zerovec(){
        return vec2{0.0,0.0};
    }

    constexpr bool operator==(const vec2& v) const {
        return (v.x == x) & (v.y == y);
    }

    friend std::ostream& operator<<(std::ostream& o, const vec2& v) {
        o << "[" << v.x << "," << v.y << "]";
        return o;
    }

};


struct mat2 {
    double a;
    double b;
    double c;
    double d;
    
    mat2() = delete;
    constexpr mat2(double _a, double _b, double _c, double _d)
    : a(_a), b(_b), c(_c), d(_d) {}
    
    constexpr mat2& operator*=(const mat2& o) {
        double ta = a * o.a + b * o.c;
        b = a * o.b + b * o.d;
        double tc = c * o.a + d * o.c;
        d = c * o.b + d * o.d;
        a = ta; c = tc;
        return *this;
    }
    constexpr friend mat2 operator*(mat2 m, const mat2& o) {
        m *= o;
        return m;
    }
    constexpr vec2 operator*(const vec2& o) const {
        return vec2{
            a * o.x + b * o.y,
            c * o.x + d * o.y
        };
    }
    constexpr mat2& operator+=(const mat2& o) {
        a += o.a;
        b += o.b;
        c += o.c;
        d += o.d;
        return *this;
    }
    constexpr friend mat2 operator+(mat2 m, const mat2& o) {
        m += o;
        return m;
    }
    constexpr mat2& operator-=(const mat2& o) {
        a -= o.a;
        b -= o.b;
        c -= o.c;
        d -= o.d;
        return *this;
    }
    constexpr friend mat2 operator-(mat2 m, const mat2& o) {
        m -= o;
        return m;
    }
    constexpr mat2& operator*=(double f) {
        a *= f; b *= f; c *= f; d *= f;
        return *this;
    }
    constexpr friend mat2 operator*(mat2 m, double f) {
        m *= f;
        return m;
    }
    constexpr friend mat2 operator*(double f, const mat2& m) {
        return m * f;
    }
    constexpr friend mat2 operator-(mat2 m) {
        m.a = -m.a;
        m.b = -m.b;
        m.c = -m.c;
        m.d = -m.d;
        return m;
    }
    static constexpr mat2 identity() {
        return mat2{1.0,0.0,0.0,1.0};
    }
    constexpr mat2 inverse() const {
        double f = 1.0 / (a*d - b*c);
        return mat2{
            d * f,
            -b * f,
            -c * f,
            a * f,
        };
    }
    constexpr mat2 transpose() const {
        return mat2{
            a, c, b, d
        };
    }
    constexpr static mat2 scale_mat(double sf) {
        return identity() * sf;
    }

    // constexpr support not until c++26
    static mat2 rot_mat(double theta) {
        return mat2{
            std::cos(theta),
            -std::sin(theta),
            std::sin(theta),
            std::cos(theta)
        };
    }
};

struct affine_2d {
    mat2 scale;
    vec2 shift;
    affine_2d() = delete;
    constexpr affine_2d(mat2 _sc, vec2 _sh) : scale(_sc), shift(_sh) {}
    
    constexpr vec2 operator()(vec2 v) const {
        v = scale*v + shift;
        return v;
    }
    constexpr friend affine_2d operator*(const mat2& m, affine_2d aff) {
        aff.scale = m * aff.scale;
        aff.shift = m * aff.shift;
        return aff;
    }
    constexpr friend affine_2d operator+(const vec2& v, affine_2d aff) {
        aff.shift += v;
        return aff;
    }
    constexpr friend affine_2d operator-(const vec2& v, affine_2d aff) {
        aff.shift -= v;
        return aff;
    }
    constexpr friend affine_2d operator*(const affine_2d& lhs, affine_2d rhs) {
        rhs = lhs.scale * rhs;
        rhs = lhs.shift + rhs;
        return rhs;
    }
    
    constexpr affine_2d inverse() const {
        mat2 matinv = scale.inverse();
        return affine_2d{
            matinv,
            matinv * -shift
        };
    }
    static constexpr affine_2d identity() {
        return affine_2d{
            mat2::identity(),
            vec2::zerovec()
        };
    }
    
};

}


#endif
