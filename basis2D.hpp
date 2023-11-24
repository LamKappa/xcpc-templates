#ifndef BASIS2D
#define BASIS2D
// 二维点集

namespace Basis2D{
    using f80 = long double;
    constexpr f80 PI = acosl(-1.L);
    constexpr f80 INF = 1e8;
    constexpr f80 EPS = 1e-8;
    int sign(f80 x){
        return std::abs(x) <= EPS ? 0 : (x > 0 ? 1 : -1);
    }
    struct Point{
        double x,y;
        Point(double _x=0, double _y=0):x(_x),y(_y){}

        Point& operator=(const Point&o){
            x = o.x; y = o.y;
            return *this;
        }
        friend std::ostream& operator<<(std::ostream&out, const Point&p){
            return out<<'('<<p.x<<','<<p.y<<')';
        }
        friend bool operator<(const Point& a, const Point&b){
            return std::abs(a.x-b.x)<EPS?a.y<b.y:a.x<b.x;
        }
        friend Point operator+(const Point& a, const Point&b){
            return {a.x+b.x, a.y+b.y};
        }
        friend Point operator-(const Point&p){
            return {-p.x, -p.y};
        }
        friend Point operator-(const Point& a, const Point&b){
            return a + (-b);
        }
        friend double operator*(const Point& a, const Point&b){
            return a.dot(b);
        }
        friend Point operator*(const Point& p, double r){
            return {r*p.x, r*p.y};
        }
        friend bool operator==(const Point& a, const Point&b){
            return a.between(b,b);
        }
        double norm()const{
            return std::hypotl(x, y);
        }
        double dot(const Point&o)const{
            return x*o.x + y*o.y;
        }
        double cross(const Point&o)const{
            return x*o.y - y*o.x;
        }
        bool between(Point a,Point b)const{
            const Point&p = *this;
            return sign((p-a).dot(p-b)) <= 0 && sign((p-b).cross(a-b)) == 0;
        }
        double theta()const{
            return std::abs(x)+std::abs(y)>EPS?atan2l(y,x):-INF;
        }
        Point rotate(double theta)const{
            return Point{
                x*cosl(theta) - y*sinl(theta),
                x*sinl(theta) + y*cosl(theta)
            };
        }
    } noPoint = {NAN, NAN};
    using Points = std::vector<Point>;

    double norm(const Point&p){
        return p.norm();
    }
    Point unit(const Point&p){
        return p * (1.L/p.norm());
    }
    double dot(const Point&a, const Point&b){
        return a.dot(b);
    }
    double cross(const Point&a, const Point&b){
        return a.cross(b);
    }
    double theta(const Point&p){
        return p.theta();
    }
    Point rotate(const Point&p, double theta){
        return p.rotate(theta);
    }
    double dist(const Point&a, const Point&b){
        return (a-b).norm();
    }
    double tan(const Point&p){
        return p.y / p.x;
    }

    struct Line : public std::array<Point,2>{
        using std::array<Point,2>::array;
        Line(Point a, Point b):array({a,b}){}

        operator Point()const{
            return (*this)[1] - (*this)[0];
        }
        int onleft(const Point&p)const{
            return sign(cross(*this, p - (*this)[0]));
        }
        friend std::ostream& operator<<(std::ostream&out, const Line&l){
            return out<<l[0]<<"->"<<l[1];
        }
        Point intersection(const Line&o)const{
            double S1 = cross(o[1]-(*this)[0], *this);
            double S2 = cross(o[0]-(*this)[0], *this);
            if(std::abs(S1-S2) < EPS) return noPoint;
            return Point{
                (S1 * o[0].x - S2 * o[1].x) / (S1 - S2),
                (S1 * o[0].y - S2 * o[1].y) / (S1 - S2)
            };
        }
        bool parallel(const Line&o)const{
            return sign(cross(*this, o)) == 0;
        }
        bool co_line(const Line&o)const{
            if(!parallel(o)) return false;
            auto&l = *this;
            return l[0].between(o[0], o[1]) || l[1].between(o[0], o[1]) ||
                o[0].between(l[0], l[1]) || o[1].between(l[0], l[1]);
        }
        bool intersection_ray(const Line&o)const{
            if(parallel(o)) return co_line(o);
            const auto&l = *this;
            int sgn = sign(cross(l, o));
            return sign(o.onleft(l[0])) == sgn && sign(l.onleft(o[0])) == -sgn;
        }
        bool intersection_strict(const Line&o)const{
            if(parallel(o)) return co_line(o);
            const auto&l = *this;
            if(l.onleft(o[0])*l.onleft(o[1]) > 0 ||
                o.onleft(l[0])*o.onleft(l[1]) > 0) {
                return false;
            }
            return true;
        }
        bool intersection_nostrict(const Line&o)const{
            if(parallel(o)) return co_line(o);
            const auto&l = *this;
            if(l.onleft(o[0])*l.onleft(o[1]) >= 0 ||
                o.onleft(l[0])*o.onleft(l[1]) >= 0) {
                return false;
            }
            return true;
        }
        Point projection(const Point&p)const{
            const Line&l = *this;
            return l[0] + (Point)l * (dot(l, p - l[0]) / dot(l, l));
        }
    } noLine = {noPoint, noPoint};
    using Lines = std::vector<Line>;

    double dist(const Point&p, const Line&l){
        double A = l[0].y - l[1].y;
        double B = l[1].x - l[0].x;
        double C = -l[0].x*A + -l[0].y*B;
        return (A*p.x + B*p.y + C) / std::sqrt(A*A + B*B);
    }
    double dist(const Line&l, const Point&p){
        return dist(p, l);
    }
    Point intersection(const Line&a, const Line&b){
        return a.intersection(b);
    }

    struct Circle{
        Point c;
        double r = 0.;
        Circle(Line l={}):c(l[0]),r(norm(l)){}
        Circle(Point c, double r):c(c),r(r){}
        Circle(const std::array<Point, 3>&triangle){
            Point p1 = (triangle[0] + triangle[1]) * 0.5;
            Point p2 = (triangle[1] + triangle[2]) * 0.5;
            this->c = Basis2D::intersection(
                {p1, p1 + rotate(triangle[0] - triangle[1], PI / 2.)},
                {p2, p2 + rotate(triangle[1] - triangle[2], PI / 2.)}
            );
            this->r = dist(this->c, triangle[1]);
        }

        operator Point()const{
            return c;
        }
        double area(){
            return PI*r*r;
        }
        friend std::ostream& operator<<(std::ostream&out, const Circle&c){
            return out<<c.c<<" radius:"<<c.r;
        }
        Point inversion(const Point&p)const{
            return c + unit(p - c) * (r*r / dist(p, c));
        }
        Circle inversion(const Line&l)const{
            Point p = l.projection(c);
            double d = norm(p-c);
            double R = r*r / (2*d);
            Point dv = unit(p-c) * R;
            return Circle{c + dv, R};
        }
        Line intersection(const Circle&o)const{
            double d = (o.c-c).norm();
            if(r + o.r < d || std::abs(r - o.r) > d) return noLine;
            double dt = acosl((r*r + d*d - o.r*o.r) / (2*d*r));
            return Line{
                c + rotate(unit(o.c-c)*r, -dt),
                c + rotate(unit(o.c-c)*r, dt)
            };
        }
        Line intersection(const Line&l)const{
            double d = dist(l, c);
            if(r < d) return noLine;
            return intersection(inversion(l));
        }
        Line tangency_line(const Point&p){
            return intersection(Circle{(p+c)*0.5, norm(p-c)});
        }
        std::vector<Line> tangency_line(const Circle&o){

        }
        std::pair<Line, Circle> inversion(const Circle&);
    } noCircle = {noPoint, NAN};
    using Circles = std::vector<Circle>;

    std::pair<Line, Circle> Circle::inversion(const Circle&o){
        double d = (o.c-c).norm();
        if(d == o.c){
            return {intersection(o), noCircle};
        }else{
            double nr = ((1./(d-o.r)) - (1./(d+o.r))) * r * r / 2.;
            return {noLine, {c + unit(o.c-c)*(nr/o.r), nr}};
        }
    }
}

#endif