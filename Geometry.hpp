#ifndef GEOMETRY_HPP
#define GEOMETRY_HPP

namespace Geometry{
    using f80 = long double;
    const f80 PI = acosl(-1.l);
    const f80 noVal = std::nanl("invalid value");
    constexpr f80 INF = 1e20l;
    constexpr f80 EPS = 1e-12l;
    int sign(f80 x){
        if(!std::isfinite(x) || fabsl(x) <= EPS) return 0;
        return x > 0 ? 1 : -1;
    }

    struct Point{
        f80 x, y;
        Point() = default;
        Point(f80 x, f80 y) : x(x), y(y) {}
        Point(const Point&) = default;

        friend std::ostream& operator<<(std::ostream&out, const Point&p){
            return out << '(' << p.x << ',' << p.y << ')';
        }
        friend bool operator<(const Point&a, const Point&b){
            return std::make_pair(a.x, a.y) < std::make_pair(b.x, b.y);
        }
        friend Point operator+(const Point&a, const Point&b){
            return {a.x + b.x, a.y + b.y};
        }
        friend Point operator-(const Point&a, const Point&b){
            return {a.x - b.x, a.y - b.y};
        }
        friend f80 operator*(const Point&a, const Point&b){
            return dot(a, b);
        }
        friend Point operator*(const Point&p, f80 r){
            return {r * p.x, r * p.y};
        }
        friend Point operator*(f80 r, const Point&p){
            return p * r;
        }
        friend Point operator/(const Point&p, f80 r){
            return {p.x / r, p.y / r};
        }
        friend bool operator==(const Point& a, const Point&b){
            return sign(a.x - b.x) == 0 && sign(a.y - b.y) == 0;
        }

        f80 square() const{
            return x * x + y * y;
        }
        f80 norm() const{
            return sqrtl(square());
        }
        Point unit() const{
            return (*this) * (1.L / norm());
        }
        friend f80 dot(const Point&a, const Point&b){
            return a.x * b.x + a.y * b.y;
        }
        friend f80 cross(const Point&a, const Point&b){
            return a.x * b.y - a.y * b.x;
        }
        friend f80 dist(const Point&a, const Point&b){
            return (a - b).norm();
        }

        bool between(const Point&a, const Point&b) const{
            auto& p = *this;
            return sign(dot(p - a, p - b)) <= 0
                   && sign(cross(p - b, a - b)) == 0;
        }
        f80 theta() const{
            return fabsl(x) + fabsl(y) > EPS ? atan2l(y, x) : -INF;
        }
        Point rotate(f80 theta) const{
            f80 cos_theta = cosl(theta), sin_theta = sinl(theta);
            return Point{
                    x * cos_theta - y * sin_theta,
                    x * sin_theta + y * cos_theta
            };
        }
        bool cmp_theta(const Point&a, const Point&b) const{
            auto& c = *this;
            auto ac = a - c, bc = b - c;
            auto cr = cross(ac, bc);
            if(cr == 0 && dot(ac, bc) >= 0){
                return ac.square() < bc.square();
            }else{
                auto xa = ac.y > 0 || (ac.y == 0 && ac.x < 0),
                        xb = bc.y > 0 || (bc.y == 0 && bc.x < 0);
                return xa == xb ? cr > 0 : xa < xb;
            }
        }
        bool cmp_angular(const Point&a, const Point&b) const;
    } noPoint = {noVal, noVal};
    using Points = std::vector<Point>;

    struct Line : public std::array<Point, 2>{
        using std::array<Point, 2>::array;
        Line() = default;
        Line(const Point&a, const Point&b) : array({a, b}) {}

        explicit operator Point()const{
            auto& L = *this;
            return L[1] - L[0];
        }
        int onLeft(const Point&p) const{
            auto& L = *this;
            return sign(cross(Point(L), p - L[0]));
        }
        friend std::ostream& operator<<(std::ostream&out, const Line&l){
            return out << '{' << l[0] << ',' << l[1] << '}';
        }
        friend Point intersection(const Line&A, const Line&B){
            f80 S1 = cross(B[1] - A[0], Point(A)), S2 = cross(B[0]-A[0], Point(A));
            if(fabsl(S1-S2) < EPS) return noPoint;
            return Point{
                (S1 * B[0].x - S2 * B[1].x) / (S1 - S2),
                (S1 * B[0].y - S2 * B[1].y) / (S1 - S2)
            };
        }
        friend bool parallel(const Line&A, const Line&B){
            return sign(cross(Point(A), Point(B))) == 0;
        }
        friend bool isCross_ray(const Line&A, const Line&B){
            if(parallel(A, B)) return A[1].between(A[0], B[0]) || B[1].between(A[0], B[0]);
            int sgn = sign(cross(Point(A), Point(B)));
            return sgn * B.onLeft(A[0]) >= 0 && sgn * A.onLeft(B[0]) <= 0;
        }
        friend bool isCross_seg(const Line&A, const Line&B){
            if(parallel(A, B)) return
                A[0].between(B[0], B[1]) || A[1].between(B[0], B[1]) ||
                B[0].between(A[0], A[1]) || B[1].between(A[0], A[1]);
            return A.onLeft(B[0]) * A.onLeft(B[1]) <= 0 &&
                   B.onLeft(A[0]) * B.onLeft(A[1]) <= 0;
            // not strict
            // return A.onLeft(B[0]) * A.onLeft(B[1]) < 0 &&
            //        B.onLeft(A[0]) * B.onLeft(A[1]) < 0;
        }
        friend bool isCross_line_to_ray(const Line&A, const Line&B){
            if(parallel(A, B)) return parallel(A, Line{A[0], B[0]});
            auto L = Line{B[0], Point(A)};
            return L.onLeft(B[1]) * L.onLeft(A[0]) > 0;
        }
        friend bool isCross_line_to_seg(const Line&A, const Line&B){
            return A.onLeft(B[0]) * A.onLeft(B[1]) <= 0;
        }
        friend bool isCross_ray_to_seg(const Line&A, const Line&B){
            return isCross_line_to_seg(A, B) &&
                   sign(dot(B[0] - A[0], Point(A))) >= 0 &&
                   sign(dot(B[1] - A[0], Point(A))) >= 0;
        }
        Point projection(const Point&p)const{
            const Line&L = *this;
            return L[0] + Point(L) * (dot(Point(L), p - L[0]) / Point(L).square());
        }
        Point symmetry(const Point&p)const{
            return projection(p) * 2 - p;
        }
        friend f80 dist(const Line&L, const Point&p){
            f80 A = L[0].y - L[1].y,
                B = L[1].x - L[0].x,
                C = - L[0].x * A - L[0].y * B;
            return (A * p.x + B * p.y + C) / sqrtl(A * A + B * B);
        }
    } noLine = {noPoint, noPoint};
    using Lines = std::vector<Line>;

    bool Point::cmp_angular(const Geometry::Point &a, const Geometry::Point &b) const{
        auto&p = *this;
        auto sgn = Line{p, a}.onLeft(b);
        return sgn == 0 ? (a - p).square() < (b - p).square() : sgn > 0;
    }

    struct Circle{
        Point c{};
        f80 r = 0.;
        Circle() = default;
        explicit Circle(const Line&L) : c(L[0]), r(Point(L).norm()) {}
        Circle(const Point&c, f80 r) : c(c), r(r) {}
        explicit Circle(const std::array<Point, 3>&triangle){
            Point p1 = (triangle[0] + triangle[1]) * 0.5;
            Point p2 = (triangle[1] + triangle[2]) * 0.5;
            this->c = intersection(
                Line{p1, p1 + (triangle[0] - triangle[1]).rotate(PI / 2.l)},
                Line{p2, p2 + (triangle[1] - triangle[2]).rotate(PI / 2.l)}
            );
            this->r = dist(this->c, triangle[1]);
        }

        explicit operator Point() const{
            return c;
        }
        f80 area() const{
            return PI * r * r;
        }
        bool inside(const Point&p) const{
            return sign(dist(p, c) - r) <= 0;
        }
        friend std::ostream& operator<<(std::ostream&out, const Circle&C){
            return out << '{' << C.c << ',' << C.r << '}';
        }
        Point inversion(const Point&p) const{
            return c + (p - c).unit() * (r * r / dist(p, c));
        }
        Circle inversion(const Line&l) const{
            Point p = l.projection(c);
            f80 d = (p - c).norm(), R = r * r / (2.l * d);
            Point dv = (p - c).unit() * R;
            return Circle{c + dv, R};
        }
        std::variant<Line, Circle> inversion(const Circle&C2) const{
            f80 d = dist(C2.c, c);
            if(sign(d - C2.r) == 0){
                return Line{intersection(*this, C2)};
            }else{
                f80 nr = ((1.l / (d - C2.r)) - (1.l / (d + C2.r))) * r * r / 2.;
                return Circle{c + (C2.c - c).unit() * (nr / C2.r), nr};
            }
        }
        friend Line intersection(const Circle&C1, const Circle&C2){
            f80 d = (C2.c - C1.c).norm();
            if(C1.r + C2.r < d || fabsl(C1.r - C2.r) > d) return noLine;
            f80 dt = acosl((C1.r * C1.r + d * d - C2.r * C2.r) / (2.l * d * C1.r));
            return Line{
                C1.c + ((C2.c - C1.c).unit() * C1.r).rotate(-dt),
                C1.c + ((C2.c - C1.c).unit() * C1.r).rotate(dt)
            };
        }
        friend Line intersection(const Circle&C, const Line&L){
            f80 d = dist(L, C.c);
            if(C.r < d) return noLine;
            return intersection(C, C.inversion(L));
        }
        std::array<Point, 2> tangency_line(const Point&p) const{
            f80 d = dist(p, c);
            if(d + EPS < r) return noLine;
            if(fabsl(d - r) < EPS) return Line{
                        p + (c - p).rotate(PI/2.l) * EPS,
                        p + (c - p).rotate(-PI/2.l) * EPS
                };
            return intersection(*this, Circle{(p + c) / 2.l, (p - c).norm() / 2.l});
        }
        //    std::vector<Line> tangency_line(const Circle&o) const{}
    } noCircle = {noPoint, noVal};
    using Circles = std::vector<Circle>;

    struct Angle : protected Point{
        using Point::Point, Point::theta;
        Angle() = default;
        Angle(const Angle&t) = default;
        explicit Angle(const Point&p) : Angle(p.x, p.y) {}
        Angle(const Angle&a, const Angle&b) : Angle(b - a) {}
        Angle(f80 _x, f80 _y) : Point(Point(_x, _y).unit()) {}
        explicit Angle(f80 theta) : Angle(cosl(theta), sinl(theta)) {}

        friend std::ostream& operator<<(std::ostream&out, const Angle&t){
            return out << '{' << (Point)t << ',' << t.theta() << '}';
        }
        friend Angle operator+(const Angle&a, const Angle&b){
            return Angle{a.x * b.x - a.y * b.y, a.x * b.y + a.y * b.x};
        }
        friend Angle operator-(const Angle&a, const Angle&b){
            return a + Angle(b.x, -b.y);
        }
        bool operator<(const Angle&t) const{
            return Point{}.cmp_theta(*this, t);
        }
    };
}

#endif