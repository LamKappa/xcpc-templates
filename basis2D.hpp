#ifndef BASIS2D
#define BASIS2D
// 二维点集

namespace Basis2D{
    constexpr double PI = acosl(-1.);
    constexpr double INF = 1e20;
    constexpr double EPS = 1e-12;
    struct Point;
    using Points = std::vector<Point>;
    using Line = std::array<Point,2>;
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
        bool operator<(const Point&o)const{
            return std::abs(x-o.x)<EPS?y<o.y:x<o.x;
        }
        Point operator+(const Point&o)const{
            return {x+o.x,y+o.y};
        }
        Point operator-()const{
            return {-x,-y};
        }
        Point operator-(const Point&o)const{
            return (*this) + (-o);
        }
        double operator*(const Point&o)const{
            return dot(o);
        }
        Point operator*(double r)const{
            return {r*x, r*y};
        }
        bool operator==(const Point&o)const{
            return between(o,o);
        }
        double norm()const{
            return std::hypotl(x, y);
        }
        double dot(const Point&o)const{
            return x*o.x+y*o.y;
        }
        double cross(const Point&o)const{
            return x*o.y-y*o.x;
        }
        bool between(Point a,Point b)const{
            if(std::abs((a.x-x)*(b.y-y)-(b.x-x)*(a.y-y))>EPS)return false;
            if(a.x>b.x)std::swap(a.x,b.x);
            if(a.y>b.y)std::swap(a.y,b.y);
            return a.x-EPS<=x&&x<=b.x+EPS &&
                a.y-EPS<=y&&y<=b.y+EPS;
        }
        double theta()const{
            return std::abs(x)+std::abs(y)>EPS?atan2(y,x):-INF;
        }
        Point rotate(double theta){
            return Point{
                x*cos(theta) - y*sin(theta),
                x*sin(theta) + y*cos(theta)
            };
        }
        double dist(const Point&o){
            return ((*this)-o).norm();
        }
        double directed_dist(const Line&o){
            double A = o[0].y - o[1].y;
            double B = o[1].x - o[0].x;
            double C = -o[0].x*A + -o[0].y*B;
            return (A*x + B*y + C) / std::sqrt(A*A + B*B);
        }
        double dist(const Line&o){
            return std::abs(directed_dist(o));
        }
    };
    Point intersection(const Line&a, const Line&b){
        double S1 = (b[1]-a[0]).cross(a[1]-a[0]);
        double S2 = (b[0]-a[0]).cross(a[1]-a[0]);
        return Point{
            (S1 * b[0].x - S2 * b[1].x) / (S1 - S2),
            (S1 * b[0].y - S2 * b[1].y) / (S1 - S2)
        };
    }
}

#endif