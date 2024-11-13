# 2D Data Structure

```cpp
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
```

# Basic 2D Algorithm

```cpp
#ifndef ALGORITHM2D
#define ALGORITHM2D
// 平面点集算法

// depends basis2D.hpp

namespace Algorithm2D{
    using namespace Basis2D;
    namespace Check_inside{
        template<typename ForwardIt>
        bool __check_inside_poly(const Point&p,
                const ForwardIt&begin,
                const ForwardIt&end){
            static std::mt19937 eng(std::random_device{}());
            static std::uniform_real_distribution<> dis(EPS,INF);
            double K = dis(eng);
            double B = p.y - K*p.x;

            int cnt = 0;
            for(auto i=begin;i!=end;i++){
                auto j = std::next(i)==end?begin:std::next(i);
                Point pi{*i}, pj{*j}, pk;

                double k = (pi.y-pj.y) / (pi.x-pj.x);
                double b = pi.y - k*pi.x;

                pk.x = (b-B) / (K-k);
                if(!finite(k)) pk.x = pi.x;
                pk.y = K*pk.x + B;

                if(pk.x>=p.x && pk.between(pi,pj))cnt++;
            }
            return cnt&1;
        }
        template<typename ForwardIt>
        bool __check_inside_chull(const Point&p,
                const ForwardIt&begin,
                const ForwardIt&end){
            for(auto i=begin;i!=end;i++){
                auto j = std::next(i)==end?begin:std::next(i);
                Point pi{*i}, pj{*j};
                if(cross(pi-p, pi-pj)) return false;
            }
            return true;
        }
        template<typename Iterable>
        bool solve(const Point&p, const Iterable&points, bool isChull = false){
            if(points.size()<3)return false;
            auto check_inside = __check_inside_poly<typename Iterable::iterator>;
            if(isChull) check_inside = __check_inside_chull<typename Iterable::iterator>;
            return check_inside(p, points.begin(), points.end());
        }
    }
    namespace Closet_pair{
        typedef std::pair<int,int> resultType;
        template<typename RandomAccessable>
        std::pair<resultType,double> __closet_pair_rec(const RandomAccessable&points){
            resultType res;
            double resv = INF;
            std::vector<int> sorted(points.size());
            std::iota(sorted.begin(),sorted.end(),0);
            std::sort(sorted.begin(),sorted.end(),[&points](auto a,auto b){
                return points[a] < points[b];
            });
            int BF_n = std::max<int>(64, std::log2(points.size()));
            std::function<void(int,int)> rec = [&](int l,int r){
                if(r-l<=BF_n){
                    std::sort(sorted.begin()+l,sorted.begin()+r,[&points](auto a,auto b){
                        return points[a].y < points[b].y;
                    });
                    for(int i=l;i<r;i++){
                        while(points[sorted[l]].y+resv<points[sorted[i]].y) l++;
                        for(int j=l;j<i;j++){
                            double tmp = (points[sorted[i]]-points[sorted[j]]).norm();
                            if(tmp<resv)resv=tmp,res={sorted[i],sorted[j]};
                        }
                    }
                }else{
                    int mid = (l+r) / 2;
                    auto&middot = points[sorted[mid]];
                    rec(l, mid); rec(mid, r);
                    std::inplace_merge(sorted.begin()+l,sorted.begin()+mid,sorted.begin()+r,
                        [&points](auto a,auto b){
                            return points[a].y < points[b].y;
                    });

                    std::array<std::deque<int>, 2> st;
                    for(int i=l;i<r;i++){
                        auto&dot = points[sorted[i]];
                        if(std::abs(dot.x-middot.x)>resv)continue;
                        bool lr = dot < middot;
                        while(!st[lr].empty()&&points[st[lr].front()].y+resv<dot.y) st[lr].pop_front();
                        for(auto it=st[lr].rbegin();it!=st[lr].rend();it++){
                            double tmp = (dot-points[*it]).norm();
                            if(tmp<resv)resv=tmp,res={sorted[i],*it};
                        }
                        st[!lr].push_back(sorted[i]);
                    }
                }
            };
            rec(0, points.size());
            return {res,resv};
        }
        template<typename RandomAccessable>
        std::pair<resultType,double> __closet_pair_multiset(const RandomAccessable&points){
            resultType res;
            double resv = INF;
            std::vector<int> sortx(points.size());
            std::iota(sortx.begin(),sortx.end(),0);
            std::sort(sortx.begin(),sortx.end(),[&points](auto a,auto b){
                return points[a] < points[b];
            });
            auto cmpy = [&points](auto a,auto b){
                return points[a].y==points[b].y?points[a].x<points[b].x:points[a].y<points[b].y;
            };
            auto sorty = sortx;
            std::sort(sorty.begin(),sorty.end(),cmpy);
            std::multiset<int,decltype(cmpy)> s(cmpy);
            for(int i=0,l=0;i<points.size();i++){
                auto&dot = points[sortx[i]];
                while(l<i&&points[l].x+resv<dot.x) s.erase(l++);
                auto it = s.lower_bound(*lower_bound(sorty.begin(),sorty.end(),dot.y-resv,
                    [&points](auto a,double b){
                        return points[a].y<b;
                }));
                while(it!=s.end()&&points[*it].y<=dot.y+resv){
                    double tmp = (dot-points[*it]).norm();
                    if(tmp<resv)resv=tmp,res={sortx[i],*it};
                    it++;
                }
                s.insert(sortx[i]);
            }
            return {res,resv};
        }
        template<typename RandomAccessable>
        std::pair<resultType,double> __solve(const RandomAccessable&points){
            return __closet_pair_rec(points);
        }
        template<typename RandomAccessable>
        resultType solve_pair(const RandomAccessable&points){
            return __solve(points).first;
        }
        template<typename RandomAccessable>
        double solve_dist(const RandomAccessable&points){
            return __solve(points).second;
        }
    }
    namespace MinCoverage{
        template<typename Iterable>
        Points solve_rectangle(const Iterable&_points){
            Points poly = Polygon::Chull::solve(_points);
            if(poly.size()<=2) return poly;
            int i=0, j=1;
            int k1=0, k2=0, k3=0;
            double ans = 1e18;
            Points ans_p(4);
            auto update = [&](){
                double H = poly[k2].dist({poly[i],poly[j]});
                double W1 = poly[k1].dist({poly[i],poly[i]+(poly[j]-poly[i]).rotate(PI/2)});
                double W2 = poly[k3].dist({poly[i],poly[i]+(poly[j]-poly[i]).rotate(PI/2)});
                double S = H*(W1+W2);
                if(S < ans){
                    Point D = unit(poly[j] - poly[i]);
                    ans_p[0] = poly[i]; ans_p[0] = ans_p[0] + D*(-W1);
                    ans_p[1] = poly[i]; ans_p[1] = ans_p[1] + D*(W2);
                    D = Point{-D.y, D.x};
                    ans_p[2] = ans_p[1]; ans_p[2] = ans_p[2] + D*(H);
                    ans_p[3] = ans_p[0]; ans_p[3] = ans_p[3] + D*(H);

                    ans = S;
                }
            };
            while((poly[j]-poly[i]).cross(poly[k1]-poly[i]) <= (poly[j]-poly[i]).cross(poly[(k1+1)%n]-poly[i])){
                k1 = (k1+1)%n;
            }
            for(; i<n; i++, j=(j+1)%n){
                while((poly[j]-poly[i]).rotate(PI/2).cross(poly[k1]-poly[i]) <= (poly[j]-poly[i]).rotate(PI/2).cross(poly[(k1+1)%n]-poly[i])){
                    k1 = (k1+1)%n;
                }
                while((poly[j]-poly[i]).cross(poly[k2]-poly[i]) <= (poly[j]-poly[i]).cross(poly[(k2+1)%n]-poly[i])){
                    k2 = (k2+1)%n;
                }
                while((poly[j]-poly[i]).rotate(PI/2).cross(poly[k3]-poly[i]) >= (poly[j]-poly[i]).rotate(PI/2).cross(poly[(k3+1)%n]-poly[i])){
                    k3 = (k3+1)%n;
                }
                update();
            }
            return ans_p;
        }
        template<typename Iterable>
        Circle solve_circle(const Iterable&_points){
            Points points(_points.begin(), _points.end());
            std::shuffle(points.begin(), points.end(), std::random_device{});
            Circle c;
            for(auto i=points.begin(); i!=points.end(); i++){
                if(dist(c, *i) <= c.r) continue;
                c = {*i, 0.};
                for(auto j=points.begin(); j!=i; j++){
                    if(dist(c, *j) <= c.r) continue;
                    c.c = (*i + *j) * 0.5;
                    c.r = dist(c, *j);
                    for(auto k=points.begin(); k!=j; k++){
                        if(dist(c, *k) <= c.r) continue;
                        c = {{*i, *j, *k}};
                    }
                }
            }
            return c;
        }
    }
    namespace HalfPlane{
        template<typename LineContainer>
        Points overlap(const LineContainer&_lines){
            Lines lines(_lines.begin(), _lines.end());
            std::sort(lines.begin(), lines.end(), [](auto a,auto b){
                return theta(a) < theta(b);
            });
            if(lines.empty() || norm(lines[0]) < EPS) return {};

            std::deque<Point> res;
            std::deque<Line> dq;
            for(auto l : lines){
                while(!res.empty() && cross(l, res.back() - l[0]) < -EPS){
                    res.pop_back(); dq.pop_back();
                }
                while(!res.empty() && cross(l, res.front() - l[0]) < -EPS){
                    res.pop_front(); dq.pop_front();
                }
                dq.push_back(l);
                if(dq.size()>1){
                    if(std::abs(cross(dq.back(), dq[dq.size()-2])) < EPS){
                        Line bk = dq.back(); dq.pop_back();
                        if(cross(bk, dq.back()[0] - bk[0]) < -EPS){
                            if(dot(bk, dq.back()) < -EPS) return {};
                            dq.back() = l;
                        }
                        if(!res.empty()) res.pop_back();
                    }
                    if(dq.size()>1) res.push_back(intersection(dq.back(), dq[dq.size()-2]));
                }
            }
            while(!res.empty() && cross(dq.front(), res.back() - dq.front()[0]) < -EPS){
                res.pop_back(); dq.pop_back();
            }
            if(dq.size()>1) res.push_front(intersection(dq.back(), dq.front()));
            assert(res.size() >= 3 && cross(dq.front(), dq.back()) <= 0);

            return Points(res.begin(), res.end());
        }
    }
}

#endif
```

# Polygon&Convex Hull Algorithm

```cpp
#ifndef POLYGON
#define POLYGON
// 平面凸包

// depends basis2D.hpp

namespace Polygon{
    using namespace Basis2D;
    namespace Chull{
        typedef Points __Chull;
        template<typename PointsContainer>
        __Chull __Graham_chull(const PointsContainer&_points){
            Points points(_points.begin(), _points.end());
            Point c = *std::min_element(points.begin(), points.end());
            std::sort(points.begin(), points.end(), [&c](auto a,auto b){
                auto ac = a - c, bc = b - c;
                auto cr = cross(ac, bc);
                if(cr == 0 && dot(ac, bc) >= 0){
                    return ac.norm() < bc.norm();
                }else{
                    auto xa = ac.y > 0 || (ac.y == 0 && ac.x < 0),
                        xb = bc.y > 0 || (bc.y == 0 && bc.x < 0);
                    return xa == xb ? cr > 0 : xa < xb;
                }
            });
            __Chull H;
            for(int i=0;i<points.size();i++){
                while(H.size()>1&&cross(H.back()-H[H.size()-2], points[i]-H.back())<=0){
                    H.pop_back();
                }
                H.push_back(points[i]);
            }
            return H;
        }
        template<typename PointsContainer>
        __Chull __Andrew_chull(const PointsContainer&_points){
            Points points(_points.begin(), _points.end());
            std::sort(points.begin(),points.end());
            std::vector<int> I;
            std::vector<bool> used(points.size(), true);
            for(int i=0;i<points.size();i++){
                while(I.size()>1&&cross(points[I.back()]-points[I[I.size()-2]], points[i]-points[I.back()])<=0){
                    used[I.back()] = false; I.pop_back();
                }
                I.push_back(i);
            }
            used[0] = false;
            std::size_t s = I.size();
            for(int i=points.size()-2;i>=0;i--){
                if(used[i]) continue;
                while(I.size()>s&&cross(points[I.back()]-points[I[I.size()-2]], points[i]-points[I.back()])<=0){
                    I.pop_back();
                }
                I.push_back(i);
            }
            __Chull H(I.size()-1);
            for(int i=0;i<I.size()-1;i++) H[i] = points[I[i]];
            return H;
        }
        template<typename Iterable>
        __Chull solve(const Iterable&points){
            if(points.size()<3)return points;
            return __Andrew_chull(points);
        }
    }
    namespace Get_diameter{
        typedef std::pair<int,int> resultType;
        template<typename RandomAccessable>
        std::pair<resultType,double> __get_diameter_rotating_calipers(const RandomAccessable&points){
            std::size_t n = points.size();
            if(n<3)return {{0,n-1},(points.back()-points.front()).norm()};
            resultType res;
            double resv = 0.;
            int j = 2;
            for(int i=0;i<n;i++){
                while(std::abs(cross(points[(i+1)%n]-points[i], points[j]-points[(i+1)%n])) <=
                    std::abs(cross(points[(i+1)%n]-points[i], points[(j+1)%n]-points[(i+1)%n]))) j = (j+1)%n;
                auto tmp = std::max<std::pair<double,int>>({(points[i]-points[j]).norm(),i},
                    {(points[(i+1)%n]-points[j]).norm(),(i+1)%n});
                if(tmp.first>resv)resv=tmp.first,res={tmp.second,j};
            }
            return {res,resv};
        }
        template<typename RandomAccessable>
        std::pair<resultType,double> __solve(const RandomAccessable&points, bool isChull){
            auto get_diameter = __get_diameter_rotating_calipers;
            if(!isChull)return get_diameter(Chull::solve(points));
            return get_diameter(points);
        }
        template<typename RandomAccessable>
        resultType solve_pair(const RandomAccessable&points, bool isChull = false){
            return __solve(points, isChull).first;
        }
        template<typename RandomAccessable>
        double solve_dist(const RandomAccessable&points, bool isChull = false){
            return __solve(points, isChull).second;
        }
    }
    namespace Area{
        template<typename Iterable>
        double __solve(const Iterable&points){
            if(points.begin() == points.end()) return 0.;
            double res = 0;
            Point o = *points.begin();
            for(auto it=points.begin(); std::next(it)!=points.end(); it++){
                Point vec = (*std::next(it) - *it);
                res += cross(*it - o, vec) / 2.;
            }
            return res;
        }
        template<typename Iterable>
        double solve(const Iterable&points, bool isChull = false){
            if(!isChull) return __solve(Chull::solve(points));
            else return __solve(points);
        }
    }
}

#endif
```

# Optimize&Search Algorithm

```cpp
#ifndef ALGORITHM_OPT_SEARCH
#define ALGORITHM_OPT_SEARCH
// 优化查找通用算法

// depends random.hpp

namespace Search{
    constexpr double PI = acosl(-1.);
    constexpr double INF = 1e20;
    constexpr double EPS = 1e-12;
    namespace Binary{
        template<typename T>
        T solve_max(T L, T R, const std::function<bool(T)>&check){
            T l = L, r = R;
            bool spj = std::is_integral<T>::value;
            while(l+EPS<r){
                T mid = (l+r) / 2;
                if(check(mid)) l = mid + spj;
                else r = mid;
            }
            return (l+r)/2;
        }
        template<typename T>
        T solve_min(T L, T R, const std::function<bool(T)>&check){
            return solve_max(L, R, [&check](T x)->bool{return !check(x);})
                + std::is_integral<T>::value;
        }
    }
    namespace Ternary{
        template<typename T>
        T solve_max(T L, T R, const std::function<double(T)>&calc){
            T l = L, r = R;
            bool spj = std::is_integral<T>::value;
            while(l+EPS<r){
                T midl = l + (r-l) / 3;
                T midr = r - (r-l) / 3;
                if(calc(midl)<calc(midr)) l = midl + spj;
                else r = midr - spj;
            }
            return (l+r)/2;
        }
        template<typename T>
        T solve_min(T L, T R, const std::function<double(T)>&calc){
            return solve_max(L, R, [&calc](T x)->double{return -calc(x);});
        }
    }
    namespace Simulated_annealing{
        using namespace Random;
        std::pair<double,double> solve_max2D(const std::function<double(double)>&func,
            double startx=0., double T=1e4, double d=0.992, double T0=1e-3, int retry=1e2){
            double x = startx;
            double res = x, resv = func(x);
            while(T>T0){
                double nxtx = x + T*(int_uniform(0,1)*2-1);
                double tmp = func(nxtx);
                double delta = resv - tmp;
                if(std::exp(-delta/T) > double_uniform(0,1)) x = nxtx;
                if(delta<0) resv=tmp,res=x;
                T *= d;
            }
            while(retry--){
                double nxtx = x + T*(int_uniform(0,1)*2-1);
                double tmp = func(nxtx);
                if(tmp>resv) resv=tmp,res=x=nxtx;
            }
            return {res,resv};
        }
        std::pair<double,double> solve_min2D(const std::function<double(double)>&func,
            double startx=0., double T=1e4, double d=0.992, double T0=1e-3, int retry=1e2){
            return solve_max2D([&func](double x){return -func(x);},
                    startx, T, d, T0, retry);
        }
        double solve_max2D_pos(const std::function<double(double)>&func,
            double startx=0., double T=1e4, double d=0.992, double T0=1e-3, int retry=1e2){
            return solve_max2D(func, startx, T, d, T0, retry).first;
        }
        double solve_max2D_value(const std::function<double(double)>&func,
            double startx=0., double T=1e4, double d=0.992, double T0=1e-3, int retry=1e2){
            return solve_max2D(func, startx, T, d, T0, retry).second;
        }
        double solve_min2D_pos(const std::function<double(double)>&func,
            double startx=0., double T=1e4, double d=0.992, double T0=1e-3, int retry=1e2){
            return solve_min2D(func, startx, T, d, T0, retry).first;
        }
        double solve_min2D_value(const std::function<double(double)>&func,
            double startx=0., double T=1e4, double d=0.992, double T0=1e-3, int retry=1e2){
            return solve_min2D(func, startx, T, d, T0, retry).second;
        }
        std::pair<std::pair<double,double>,double> solve_max3D(const std::function<double(double,double)>&func,
            double startx=0., double starty=0., double T=1e4, double d=0.992, double T0=1e-3, int retry=1e2){
            typedef std::pair<double,double> resultType;
            double x = startx, y = starty;
            resultType res = {x,y};
            double resv = func(x,y);
            while(T>T0){
                double theta = double_uniform(-PI,PI);
                double nxtx = x + T*cos(theta);
                double nxty = y + T*sin(theta);
                double tmp = func(nxtx,nxty);
                double delta = resv - tmp;
                if(std::exp(-delta/T) > double_uniform(0,1)) x = nxtx, y = nxty;
                if(delta<0) resv=tmp,res={x,y};
                T *= d;
            }
            x = res.first; y = res.second;
            while(retry--){
                double theta = double_uniform(-PI,PI);
                double nxtx = x + T0*cos(theta);
                double nxty = y + T0*sin(theta);
                double tmp = func(nxtx,nxty);
                if(tmp>resv) resv=tmp,res={x=nxtx,y=nxty};
            }
            return {res,resv};
        }
        std::pair<std::pair<double,double>,double> solve_min3D(const std::function<double(double,double)>&func,
            double startx=0., double starty=0., double T=1e4, double d=0.992, double T0=1e-3, int retry=1e2){
            return solve_max3D([&func](double x,double y){return -func(x,y);},
                    startx, starty, T, d, T0, retry);
        }
        std::pair<double,double> solve_max3D_pos(const std::function<double(double,double)>&func,
            double startx=0., double starty=0., double T=1e4, double d=0.992, double T0=1e-3, int retry=1e2){
            return solve_max3D(func, startx, starty, T, d, T0, retry).first;
        }
        double solve_max3D_value(const std::function<double(double,double)>&func,
            double startx=0., double starty=0., double T=1e4, double d=0.992, double T0=1e-3, int retry=1e2){
            return solve_max3D(func, startx, starty, T, d, T0, retry).second;
        }
        std::pair<double,double> solve_min3D_pos(const std::function<double(double,double)>&func,
            double startx=0., double starty=0., double T=1e4, double d=0.992, double T0=1e-3, int retry=1e2){
            return solve_min3D(func, startx, starty, T, d, T0, retry).first;
        }
        double solve_min3D_value(const std::function<double(double,double)>&func,
            double startx=0., double starty=0., double T=1e4, double d=0.992, double T0=1e-3, int retry=1e2){
            return solve_min3D(func, startx, starty, T, d, T0, retry).second;
        }
    }
}

#endif
```

# Fenwick Tree Data Structure

```cpp
#ifndef FENWICK_TREE
#define FENWICK_TREE
// 树状数组

int lowbit(int x){return x&-x;}

template<typename T=int>
struct BIT{
    int n;
    std::vector<T> node;
    BIT(int n=0){init(n);}
    BIT(const std::vector<T>&initarr){
        init(initarr);
    }
    void init(const std::vector<T>&initarr){
        init(initarr.size());
        for(int i=1; i<=n; i++){
            node[i] += initarr[i-1];
            if(i+lowbit(i)<=n)node[i+lowbit(i)] += node[i];
        }
    }
    void init(int n){
        this->n = n;
        node.assign(n+1,T());
    }
    void add(int k, T v){
        if(!k)return;
        for(;k<=n;k+=lowbit(k))node[k]+=v;
    }
    T ask(int k){
        T ret = 0;
        for(;k>0;k-=lowbit(k))ret+=node[k];
        return ret;
    }
    T ask(int l, int r){
        if(l>r)std::swap(l,r);
        return ask(r) - ask(l-1);
    }
    int kth(T k){
        int x = 0;
        for(int i=1<<std::__lg(n); i; i>>=2){
            if(x+i<=n && k>=node[x+i-1]){
                x += i;
                k -= node[x-1];
            }
        }
        return x;
    }
};

template<typename T=int>
struct BIT_range{
    int n;
    BIT<T> di,dn;
    BIT_range(int n=0):di(n),dn(n){}
    BIT_range(const std::vector<T>&initarr){
        init(initarr);
    }
    void init(const std::vector<T>&initarr){
        init(initarr.size());
        for(int i=1; i<=n; i++){
            di.node[i] += i*initarr[i-1];
            if(i<n) di.node[i+1] -= (i+1)*initarr[i-1];
            if(i+lowbit(i)<=n) di.node[i+lowbit(i)] += di.node[i];
            dn.node[i] += initarr[i-1];
            if(i<n) dn.node[i+1] -= initarr[i-1];
            if(i+lowbit(i)<=n) dn.node[i+lowbit(i)] += dn.node[i];
        }
    }
    void init(int n){
        this->n = n;
        di.init(n+1); dn.init(n+1);
    }
    void __add(int k, int v){
        di.add(k,k*v);
        dn.add(k,v);
    }
    void add(int l, int r, T v){
        if(l>r) std::swap(l,r);
        __add(l,v); __add(r+1,-v);
    }
    void add(int k,T v){add(k,k,v);}
    T __ask(int k){
        return (k+1)*dn.__ask(k) - di.__ask(k);
    }
    T ask(int l, int r){
        if(l>r) std::swap(l,r);
        return __ask(r) - __ask(l-1);
    }
};

#endif
```

# Disjoint Set Union Data Structure

```cpp
#ifndef DISJOINT_SET_UNION
#define DISJOINT_SET_UNION
// 并查集

struct DSU {
    std::vector<int> fa, sz, st;
    DSU() {}
    DSU(int n) {init(n);}
    void init(int n) {
        fa.resize(n);
        std::iota(fa.begin(), fa.end(), 0);
        sz.assign(n, 1);
        st.clear();
    }
    int find(int x) {
        while(x!=fa[x]) {
            x = fa[x] = fa[fa[x]];
        }
        return x;
    }
    bool same(int x, int y) {
        return find(x) == find(y);
    }
    bool merge(int x, int y) {
        x = find(x); y = find(y);
        if(x==y) { return false; }
        if(sz[x]<sz[y]) { std::swap(x,y); }
        sz[x] += sz[y];
        st.push_back(y);
        fa[y] = x;
        return true;
    }
    int get_size(int x) {
        return sz[find(x)];
    }
    void clear() {
        while(!st.empty()) {
            int x = st.back(), p = find(x);
            fa[p] = p; sz[p] = 1;
            fa[x] = x; sz[x] = 1;
            st.pop_back();
        }
    }
};

#endif
```

# Edge Biconnected Component Algorithm

```cpp
#ifndef EDGE_BICONNECTED_COMPONENT
#define EDGE_BICONNECTED_COMPONENT

std::set<std::pair<int, int>> E;

struct EBCC{
    int n, cur, cnt;
    std::vector<std::vector<int>> adj;
    std::vector<int> stk;
    std::vector<int> dfn, low, bel;
    
    EBCC() = default;
    EBCC(int n){init(n);}
    
    void init(int n){
        this->n = n;
        adj.assign(n, {});
        dfn.assign(n, -1);
        low.resize(n);
        bel.assign(n, -1);
        stk.clear();
        cur = cnt = 0;
    }
    
    void add_edge(int u, int v){
        adj[u].push_back(v);
        adj[v].push_back(u);
    }
    
    void dfs(int u, int fa){
        dfn[u] = low[u] = cur++;
        stk.push_back(u);
        
        for (auto v : adj[u]){
            if(v == fa){
                fa = -1;
                continue;
            }
            if(dfn[v] == -1){
                E.emplace(u, v);
                dfs(v, u);
                low[u] = std::min(low[u], low[v]);
            }else if(bel[v] == -1 && dfn[v] < dfn[u]){
                E.emplace(u, v);
                low[u] = std::min(low[u], dfn[v]);
            }
        }
        
        if(dfn[u] == low[u]){
            for(int v=-1; v!=u; stk.pop_back()){
                v = stk.back();
                bel[v] = cnt;
            }
            cnt++;
        }
    }
    
    std::vector<int> work(){
        for(int i=0; i<n; i++){
            if(dfn[i] != -1) continue;
            dfs(i, -1);
        }
        return bel;
    }
};

#endif
```

# Euler Loop Algorithm

```cpp
#ifndef EULER_LOOP
#define EULER_LOOP

struct EulerLoop{
    std::vector<std::unordered_map<int, int>> adj;

    EulerLoop(int n){
        adj.assign(n, {});
    }

    void add_edge(int u, int v){
        adj[u][v]++;
    }

    std::vector<int> find(int rt){
        std::vector<int> ans;
        // auto adj = this->adj;
        auto dfs = [&](auto&&dfs, int u)->void{
            while(!adj[u].empty()){
                auto v = adj[u].begin()->first;
                if(0 == --adj[u][v]){
                    adj[u].erase(v);
                }
                dfs(dfs, v);
            }
            ans.push_back(u);
        };
        dfs(dfs, rt);
        reverse(ans.begin(), ans.end());
        return ans;
    }
    
    bool empty()const{
        bool fl = true;
        for(const auto&x : adj){
            fl &= x.empty();
        }
        return fl;
    }
};

#endif
```

# Fraction Data Structure

```cpp
#ifndef FRACTION
#define FRACTION

using i64 = long long;
using i128 = __int128;
using f80 = long double;
const i64 INF = -1ULL >> 15;

std::ostream& operator<<(std::ostream&out, i128 x){
    static int st[33];
    int top = 0;
    do{st[top++] = x%10;}while(x/=10);
    while(top>0) out<<st[--top];
    return out;
}

std::array<std::pair<i64,i64>, 2> fraction(f80 x, i64 M = INF, i64 N = INF){
    i64 m = 1, n = 1, lm=0, ln=1, rm=1, rn=0;
    while(m<=M && n<=N){
        int k = 0, fl = x*n > m;
        if(fl) lm = m, ln = n;
        else rm = m, rn = n;
        while(k>=0){
            if(fl) m = lm + (rm<<k), n = ln + (rn<<k);
            else m = (lm<<k) + rm, n = (ln<<k) + rn;
            if(m>M || n>N) break;
            if((x*n > m) == fl){
                if(fl) lm = m, ln = n;
                else rm = m, rn = n;
                k++;
            }else k--;
        }
        m=lm+rm; n=ln+rn;
    }
    return {std::make_pair(lm, ln), std::make_pair(rm, rn)};
}
struct Fraction{
    i64 a = 0, b = 1;
    i128 c = 0;

    Fraction(f80 x){
        auto[p1, p2] = fraction(std::abs(x));
        if(p1.first > p2.first) std::swap(p1, p2);
        if(p1.second == 0) std::swap(p1, p2);
        *this = Fraction((x<0?-1:1) * p1.first, p1.second);
    }
    Fraction(i128 a, i128 b, i128 c = 0){
        this->c = c + a / b;
        a %= b;
        constexpr i64 MAX_I64 = -1ULL>>1;
        auto mx = std::max(std::abs(a), std::abs(b));
        int l = 0, r = 64;
        if(mx <= MAX_I64) r = 0;
        while(l < r){
            int mid = (l+r) / 2;
            if((mx >> mid) <= MAX_I64){
                r = mid;
            }else l = mid + 1;
        }
        pack(a>>l, b>>l);
    }
    void pack(i64 a, i64 b){
        bool neg = (a<0) ^ (b<0);
        a = std::abs(a), b = std::abs(b);
        i64 g = std::__gcd(a, b);
        if(neg) a = -a;
        this->a = a / g;
        this->b = b / g;
    }

    operator f80()const{
        return c + (f80)a / b;
    }
    friend std::ostream& operator<<(std::ostream&out, const Fraction&f){
        return out<<"("<<f.c<<","<<f.a<<"/"<<f.b<<")";
    }
    Fraction operator+(const Fraction&o)const{
        i64 g = std::__gcd(b, o.b);
        return Fraction{a*(o.b/g) + o.a*(b/g), b*(o.b/g), c + o.c};
    }
    Fraction operator-(const Fraction&o)const{
        return *this + (-o);
    }
    Fraction operator-()const{
        return Fraction{-a, b, -c};
    }
    Fraction operator*(const Fraction&o)const{
        i64 g1 = std::__gcd(a, o.b);
        i64 g2 = std::__gcd(b, o.a);
        return Fraction{(a/g1) * (o.a/g2), (b/g2) * (o.b/g1), c * o.c} +
            Fraction{c*o.a, o.b} + Fraction{o.c*a, b};
    }
    Fraction operator/(const Fraction&o)const{
        return *this * o.inv();
    }
    Fraction inv()const{
        return Fraction{b, c*b + a};
    }

    Fraction operator+(const f80&o)const{
        return *this + Fraction(o);
    }
    Fraction operator-(const f80&o)const{
        return *this - Fraction(o);
    }
    Fraction operator*(const f80&o)const{
        return *this * Fraction(o);
    }
    Fraction operator/(const f80&o)const{
        return *this / Fraction(o);
    }
};

#endif
```

## Stern-Brocot树

用于解决：已知实数$r$求一组足够接近（分子不超过$N$分母不超过$M$）的最简有理数$\frac{q}{p}$其中$\gcd(p,q)=1$

初始两个分数$\frac{0}{1}$和$\frac{1}{0}$，然后每次对于相邻的两个分数$\frac{m}{n}$和$\frac{m'}{n'}$，把$\frac{m+m'}{n+n'}$插入到它们中间可以证明这棵树具有二叉搜索性质且遍历了所有非负最简有理数，于是对于某个实数$r$只需如此寻找。

# Heavy-Light Decomposition Data Structure

```cpp
#ifndef HLD
#define HLD
// 重链剖分

// depends SegTree.hpp

template<class Info, class Tag>
struct HLD{
    int n;
    std::vector<int> sz, top, dpt, fa, in, out, seq;
    std::vector<std::vector<int>> adj;
    SegTree<Info, Tag> tree;
    
    HLD(int n=0){init(n);}
    void init(int n){
        this->n = n;
        sz.resize(n); top.resize(n); dpt.resize(n);
        fa.resize(n); in.resize(n); out.resize(n);
        seq.resize(n); adj.assign(n,{});
    }
    void init_val(const std::vector<Info>&initarr){
        std::vector<Info> mapping(n);
        for(int i=0;i<n;i++){
            mapping[in[i]] = initarr[i];
        }
        tree.init(mapping);
    }
    void add_edge(int u,int v){
        adj[u].push_back(v);
        adj[v].push_back(u);
    }
    void work(int rt=0){
        top[rt] = rt; fa[rt] = -1;
        int cur = dpt[rt] = 0;
        std::function<void(int)> dfs1 = [&](int u){
            sz[u] = 1;
            for(auto v : adj[u]){
                if(fa[u] == v) continue;
                fa[v] = u;
                dpt[v] = dpt[u] + 1;
                dfs1(v);
                sz[u] += sz[v];
            }
            if(adj[u].empty()) return;
            std::swap(*std::max_element(adj[u].begin(), adj[u].end(), [&](int x,int y){
                return fa[u]==x || sz[x]<sz[y];
            }), adj[u].front());
        };
        std::function<void(int)> dfs2 = [&](int u){
            in[u] = cur++;
            seq[in[u]] = u;
            for(auto v : adj[u]){
                if(fa[u] == v) continue;
                top[v] = (v==adj[u][0] ? top[u] : v);
                dfs2(v);
            }
            out[u] = cur;
        };
        dfs1(rt); dfs2(rt);
    }

    void apply(int u, int v, const Tag&t){
        while(top[u]!=top[v]){
            if(dpt[top[u]] < dpt[top[v]]) std::swap(u, v);
            tree.apply(in[top[u]], in[u], t);
            u = fa[top[u]];
        }
        tree.apply(in[u], in[v], t);
    }

    void apply_sub(int rt, const Tag&t){
        tree.apply(in[rt], out[rt]-1, t);
    }
    
    Info ask(int u, int v){
        Info res = Info();
        while(top[u]!=top[v]){
            if(dpt[top[u]] < dpt[top[v]]) std::swap(u, v);
            res = res + tree.ask(in[top[u]], in[u]);
            u = fa[top[u]];
        }
        return res + tree.ask(in[u], in[v]);
    }

    Info ask_sub(int rt){
        return tree.ask(in[rt], out[rt]-1);
    }

    int lca(int u, int v){
        while(top[u] != top[v]){
            if(dpt[top[u]] < dpt[top[v]]) std::swap(u, v);
            u = fa[top[u]];
        }
        return dpt[u]<dpt[v] ? u : v;
    }
    
    int dist(int u, int v){
        return dpt[u] + dpt[v] - dpt[lca(u, v)]*2;
    }

    int jump(int u, int k){
        if(dpt[u] < k) return -1;
        int d = dpt[u]-k;
        while(dpt[top[u]] > d){
            u = fa[top[u]];
        }
        return seq[in[u] - dpt[u]+d];
    }
    
    bool check_anc(int u, int v){
        return in[u]<=in[v] && in[v]<out[u];
    }
    
    int rooted_fa(int u, int rt){
        if(u == rt) return -1;
        if(!check_anc(u, rt)) return fa[u];
        auto it = std::upper_bound(adj[u].begin(), adj[u].end(), rt, [&](int x, int y){
            return in[x]<in[y];
        }) - 1;
        return *it;
    }
    
    int rooted_size(int u, int rt){
        if(u == rt) return n;
        if(!check_anc(u,rt)) return sz[u];
        return n - sz[rooted_fa(u, rt)];
    }
    
    int rooted_lca(int a, int b, int rt){
        return lca(a, b) ^ lca(b, rt) ^ lca(rt, a);
    }
};

#endif
```

# Linear Planning Algorithm

```cpp
#ifndef LINEAR_PROGRAMMING
#define LINEAR_PROGRAMMING
// 线性规划 支持无解、无穷大判定 常数大 缓存命中差

enum Relation {
    LESS_EQ,
    EQUAL,
    GREATER_EQ
};

struct LP {
    constexpr static double EPS = 1e-5;
    std::size_t n = 0, m = 0;
    std::vector<int> base, varmp, baseX;
    std::vector<double> b, c;
    std::vector<std::vector<double>> A;
    
    bool exist = true;
    
    // stdandard: Ax <= b maximize cx
    // assert xi >= 0
    // initialize: Ajx + x{n+j} == bj
    
    void add_restrict(std::vector<double> aj, Relation r, double bj) {
        n = std::max(n, aj.size());
        if(r != GREATER_EQ) {
            A.push_back(aj); b.push_back(bj); m++;
        }
        if(r != LESS_EQ) {
            A.push_back(aj); b.push_back(EPS-bj); m++;
            for(auto&x : A.back()) { x = -x; }
        }
    }
    
    void pivot(int x, int j) {
        baseX[j] = x;
        double Ajx = A[j][x], cx = c[varmp[x]];
        b[j] /= Ajx;
        c[n + m + 1] += cx * b[j];
        for(int i=0; i<n+1; i++) {
            c[varmp[i]] -= cx * (A[j][i] /= Ajx);
        }
        c[base[j]] -= cx / Ajx;
        A[j][x] = 1. / Ajx;
        
        for(int jj=0; jj<m; jj++) {
            if(jj == j) { continue; }
            if(std::abs(A[jj][x]) < EPS) { continue; }
            for(int i=0; i<n+1; i++) {
                if(i!=x) { A[jj][i] -= A[jj][x] * A[j][i]; }
            }
            b[jj] -= A[jj][x] * b[j];
            A[jj][x] /= -Ajx;
        }
        std::swap(varmp[x], base[j]);
    }
    
    void simplex() {
        for(;;) {
            int max_i = -1, min_j = -1;;
            for(int i=0; i<n+1; i++) {
                if(c[varmp[i]] < EPS) { continue; }
                if(max_i < 0 || c[varmp[max_i]] < c[varmp[i]]) {
                    max_i = i;
                }
            }
            if(max_i < 0) { return; }
            for(int j=0; j<m; j++) {
                if(A[j][max_i] < EPS) { continue; }
                if(min_j < 0 ||
                        b[j]*A[min_j][max_i] < b[min_j]*A[j][max_i]) {
                    min_j = j;
                }
            }
            if(min_j < 0) { return (void)(exist = false); }
            pivot(max_i, min_j);
        }
    }
    
    void initialize() {
        base.resize(m);
        baseX.resize(m);
        varmp.resize(n+1);
        std::iota(varmp.begin(),varmp.end(),0);
        varmp[n] = n + m;
        int min_j = 0;
        for(int j=0; j<m; j++) {
            A[j].resize(n+1, 0.);
            base[j] = n + j;
            baseX[j] = n;
            
            if(b[j] < b[min_j]) {
                min_j = j;
            }
        }
        if(b[min_j] < 0) {
            for(int j=0; j<m; j++) {
                A[j][n] = -1.;
            }
            decltype(c) _c; std::swap(_c, c);
            c.assign(n+m+2, 0.); c[n + m] = -1.;
            pivot(n, min_j); simplex();
            if(std::abs(c[n + m] + 1) > EPS) {
                return (void)(exist = false);
            }
            std::swap(c, _c);
            for(int j=0; j<m; j++) {
                for(int i=0; i<n+1; i++) {
                    if(varmp[i] == n + m) {
                        A[j][i] = 0.;
                    }
                    c[varmp[i]] -= c[base[j]] * A[j][i];
                }
                c[n + m + 1] += c[base[j]] * b[j];
                c[base[j]] = 0.;
            }
            for(int j=0, i; j<m; j++) {
                if(base[j] != n + m) { continue; }
                for(i=0; i<n+1 && std::abs(A[j][i])<EPS; i++);
                pivot(i, j);
                for(i=0; i<n+1 && varmp[i] != n+m; i++);
                A[j][i] = 0.;
                break;
            }
        }
    }
    
    double maximize(std::vector<double> _c) {
        exist = true;
        c = _c; c.resize(n+m+2, 0.);
        c[n + m + 1] = c[n]; c[n] = 0.;
        initialize();
        if(!exist) { return NAN; }
        simplex();
        if(!exist) { return 1. / 0.; }
        return c[n + m + 1];
    }
    
    double minimize(std::vector<double> _c) {
        for(auto&x : _c) { x = -x; }
        return -maximize(_c);
    }
    
};

#endif

// e.g.

// LP lp1;
// lp1.add_restrict({1.,1.,3.}, LESS_EQ, 30.);
// lp1.add_restrict({2.,2.,5.}, LESS_EQ, 24.);
// lp1.add_restrict({4.,1.,2.}, LESS_EQ, 36.);
// cout<<lp1.maximize({3.,1.,2.})<<endl;    // 28

// LP lp2;
// lp2.add_restrict({2.,-1.}, LESS_EQ, 2.);
// lp2.add_restrict({1.,-5.}, LESS_EQ, -4.);
// cout<<lp2.maximize({2.,-1.})<<endl;      // 2

// LP lp3;
// lp3.add_restrict({1,-1}, LESS_EQ, -2);
// cout<<lp3.maximize({1})<<endl;              // inf

// LP lp4;
// lp4.add_restrict({1}, LESS_EQ, -1);
// lp4.add_restrict({1}, GREATER_EQ, 1);
// cout<<lp4.maximize({1})<<endl;              // nan
```

# Mathematics Algorithm

```cpp
#ifndef MATH
#define MATH
// 数学

namespace Math {
    using i64 = long long;
    using i128 = __int128;
    const i64 INF = -1ULL >> 3;
    const long double EPS = 1e-16;

    namespace Sieve {
        std::vector<i64> minp, primes, phi, mu;
        
        void init(int N) {
            int last_n = std::max<int>(2, minp.size());
            if(primes.empty()){
                minp.assign(N+1, 0);
                phi.assign(N+1, 0);
                mu.assign(N+1, 0); mu[1] = 1;
                primes.clear();
            }
            
            for(int n=last_n; n<=N; n++) {
                if(0==minp[n]) {
                    minp[n] = n;
                    phi[n]  = n-1;
                    mu[n]   = -1;
                    primes.push_back(n);
                }
                for(int p : primes) {
                    if(n * p > N) { break; }
                    minp[n * p] = p;
                    if(p==minp[n]) {
                        phi[n * p]  = phi[n] * p;
                        mu[n * p]   = 0;
                        break;
                    }
                    phi[n * p]  = phi[n] * (p-1);
                    mu[n * p]   = mu[n] * mu[p];
                }
            }
        }
        
        std::vector<std::array<int,3>> FAC_3;
        std::vector<std::vector<int>> linear_gcd_lst;
        void linear_gcd_init(int N) {
            FAC_3.resize(N+1);
            minp.assign(N+1, 0);
            primes.clear();
            
            for(int n=2; n<=N; n++) {
                if(0==minp[n]) {
                    FAC_3[n]= {1, 1, n};
                    minp[n] = n;
                    primes.push_back(n);
                }
                for(int p : primes) {
                    if(n * p > N) { break; }
                    FAC_3[n * p] = FAC_3[n]; FAC_3[n * p][0] *= p;
                    std::sort(FAC_3[n * p].begin(), FAC_3[n * p].end());
                    minp[n * p] = p;
                    if(p==minp[n]) { break; }
                }
            }
            int sqN = std::sqrt(N);
            linear_gcd_lst.resize(sqN+1);
            for(int i=1; i<=sqN; i++) {
                linear_gcd_lst[i].resize(i+1);
                linear_gcd_lst[i][0] = i;
                for(int j=1; j<=i; j++) {
                    linear_gcd_lst[i][j] = linear_gcd_lst[j][i % j];
                }
            }
        }
        int linear_gcd(int a, int b) {
            if(a > b) { std::swap(a, b); }
            if(a <= 1) { return a ? 1 : b; }
            int res = 1;
            for(auto k : FAC_3[a]) {
                int tmp = 1;
                if(k > linear_gcd_lst.size()) {
                    if(b % k == 0) { tmp = k; }
                } else { tmp = linear_gcd_lst[k][b % k]; }
                b /= tmp;
                res *= tmp;
            }
            return res;
        }
    }
    
    namespace Basic {
        i64 gcd(i64 a, i64 b) {
            i64 neg = ((a < 0) ^ (b < 0)) ? -1 : 1;
            a = std::abs(a); b = std::abs(b);
            if((a&-a) < (b&-b)) { std::swap(a, b); }
            if(b == 0) { return a; }
            int bz = __builtin_ctzll(b);
            b >>= bz;
            while(a) {
                a >>= __builtin_ctzll(a);
                i64 tmp = b - a;
                if(a < b) { b = a; }
                a = std::abs(tmp);
            }
            return neg * (b << bz);
        }
        // a*x + b*y = gcd(a, b) returns min positive x
        std::array<i64, 3> exgcd(i64 a, i64 b){
            auto __exgcd = [](auto&&__exgcd, auto a, auto b, auto&x, auto&y)->i64{
                if(b==0){
                    x = 1; y = 0; return a;
                }
                auto d = __exgcd(__exgcd, b, a%b, x, y);
                std::tie(x, y) = std::make_pair(y, x - a / b * y);
                return d;
            };
            i64 x, y;
            auto d = __exgcd(__exgcd, a, b, x, y);
            auto tx = b / d;
            auto ty = a / d;
            x = x + tx * (i64)ceil((1. - x) / tx);
            y = (d - a * x) / b;
            return {x, y, d};
        }
        i64 qmul(i64 a, i64 b, i64 mod) {
            i64 res = a*b - mod*(i64)(1.L/mod*a*b);
            return res - mod*(res>=mod) + mod*(res<0);
        }
        i64 qpow(i64 a, int b){
            i64 res = 1;
            while(b) {
                if(b & 1) res *= a;
                a *= a;
                b >>= 1;
            }
            return res;
        }
        i64 qpow(i64 a, i64 b, i64 mod, i64 phi_mod = 0) {
            a %= mod;
            if(phi_mod > 0) b = std::min(b, b % phi_mod + phi_mod);
            i64 res = 1ll;
            while(b) {
                if(b & 1) { res = qmul(res, a, mod); }
                a = qmul(a, a, mod);
                b >>= 1;
            }
            return res;
        }
        bool miller_rabin(i64 x) {
            constexpr std::array<i64, 7> test_p = {2, 325, 9375, 28178, 450775, 9780504, 1795265022ll};
            if(x<=3 || x%2 == 0) { return x==2 || x==3; }
            if(x%6 != 1 && x%6 != 5) { return false; }
            i64 k = x-1;
            int r = __builtin_ctzll(k), j;
            k >>= r;
            for(auto p : test_p) {
                if(p % x == 0) { continue; }
                i64 v = qpow(p,k,x);
                if(v == 1) { continue; }
                for(j=1; j<=r; j++, v=qmul(v,v,x)) if(v==x-1) { break; }
                if(j > r) { return false; }
            }
            return true;
        }
        i64 pollard_rho(i64 num) {
            if(num % 2 == 0) { return 2; }
            i64 val = 1ll, s = 0, t = 0;
            static std::mt19937 eng(std::random_device{}());
            i64 c = std::uniform_int_distribution<i64>(1ll, num-1)(eng);
            auto func = [](i64 x,i64 c,i64 mod) {
                return (qmul(x,x,mod) + c) % mod;
            };
            for(int goal=1; ; goal<<=1, s=t, val=1ll) {
                for(int steps=1; steps<=goal; ++steps) {
                    t = func(t, c, num);
                    val = qmul(val, std::abs(t-s), num);
                    if(steps % 127 == 0) {
                        i64 gcd_ = gcd(val, num);
                        if(gcd_ > 1) { return gcd_; }
                    }
                }
                i64 gcd_ = gcd(val, num);
                if(gcd_ > 1) { return gcd_; }
            }
            assert(false);
        }
        i64 max_prime_factor(i64 num) {
            i64 max_factor = 1;
            auto dfs = [&max_factor](auto dfs, i64 num) {
                if(num <= max_factor || num < 2) { return; }
                if(miller_rabin(num)) {
                    return (void)(max_factor = std::max(max_factor, num));
                }
                i64 factor_ = pollard_rho(num);
                while(factor_ >= num) { factor_ = pollard_rho(num); }
                while(num % factor_ == 0) { num /= factor_; }
                dfs(dfs, num); dfs(dfs, factor_);
            };
            dfs(dfs, num);
            return max_factor;
        }
        std::vector<std::pair<i64, int>> factorize(i64 x){
            std::vector<std::pair<i64, int>> fac;
            while(x > 1){
                i64 p = (Sieve::minp.size() > x) ? Sieve::minp[x] : max_prime_factor(x), cnt = 0;
                for(auto pp=p, k=1ll; x%pp==0 || ((pp=p,k=1ll) && x%p==0); pp*=pp, k<<=1) x/=pp, cnt+=k;
                fac.emplace_back(p, cnt);
            }
            return fac;
        }
        i64 phi(i64 x, const std::vector<std::pair<i64, int>>&fac={}) {
            if(Sieve::phi.size() > x) return Sieve::phi[x];
            i64 res = x;
            for(auto [p, _] : fac.empty() ? factorize(x) : fac){
                res = res / p * (p-1);
            }
            return res;
        }
        i64 inv(i64 a, i64 mod) {
            auto[x, y, d] = exgcd(a, mod);
            if(d > 1) return -1;
            return x;
            // return qpow(a, phi(mod) - 1, mod);
        }
        std::array<i64, 2> primitive_root(i64 x) {
            // x is 2 || 4 || prime^k || 2*prime^k
            auto fac_x = factorize(x);
            if(fac_x.size() > 1){
                if(fac_x[0].first > fac_x[1].first) std::swap(fac_x[0], fac_x[1]);
                if(fac_x.size() != 2 || fac_x[0].second != 1){
                    return {-1};
                }
            }else if(x % 2 == 0){
                if(x <= 4) return {x-1, x/2};
                return {-1};
            }
            auto phi_x = phi(x, fac_x);
            auto fac_phi_x = factorize(phi_x);
            for(i64 g=2; g<x; g++){
                bool fl = true;
                if(gcd(g, x) > 1) continue;
                for(auto [fac, _] : fac_phi_x){
                    if(qpow(g, phi_x/fac, x) == 1){
                        fl = false;
                        break;
                    }
                }
                if(fl) return {g, phi_x};
            }
            return {-1};
        }
        // x == ai (mod mi)
        std::array<i64, 2> exCRT(const std::vector<std::pair<i64, i64>>&eq){
            if(eq.empty()) return {-1};
            auto[a0 ,m0] = eq.front(); a0 %= m0;
            for(int i=1; i<eq.size(); i++){
                auto[ai, mi] = eq[i]; ai %= mi;
                // (ai - a0) == x * m0 - y * mi
                auto[x, y, d] = exgcd(m0, mi);
                if((ai - a0) % d != 0) return {-1};
                mi = m0 / d * mi;
                x = qmul((ai - a0) / d, x, mi);
                a0 = (a0 + qmul(x, m0, mi)) % mi;
                m0 = mi;
            }
            return {a0, m0};
        }
        // ax == b (mod m)
        i64 liEu(i64 a, i64 b, i64 m){
            a %= m; b %= m;
            if(a == 0) return b == 0 ? 0 : -1;
            auto d = gcd(a, m);
            if(b % d != 0) return -1;
            // x = x0 + i * n / d
            return (b / d) * inv(a / d, m / d) % (m / d);
        }
        // c * a^x == b (mod m) exBSGS
        i64 exBSGS(i64 a, i64 b, i64 m, i64 c = 1) {
            a %= m; b %= m; c %= m;
            if(b == c || m == 1) return 0;
            // a^k/D * a^(x-k) == b/D (mod m/D)
            i64 k = 0, D = 1, akD = 1;
            for(i64 d; (d=gcd(a, m/D)) > 1; ){
                if(b/D % d != 0) return -1;
                k++; D *= d; akD = qmul(akD, (a / d), m);
                if(qmul(c, qmul(akD, D, m), m) == b){
                    return k;
                }
            }
            m = m / D; b = b / D;
            a %= m; b %= m; c = qmul(c, akD, m);
            
            i64 sqm = std::ceil(std::sqrt(m));

            std::unordered_map<i64, i64> baby_steps;
            for(i64 i=0, ax=1; i<sqm; i++){
                i64 baby_step = qmul(ax, b, m);
                baby_steps[baby_step] = i;
                ax = qmul(ax, a, m);
            }
            i64 x = -1;
            for(i64 j=1, asqm=qpow(a, sqm, m), ax=asqm; j<=sqm; j++){
                i64 giant_step = qmul(ax, c, m);
                if(baby_steps.count(giant_step)){
                    x = j * sqm - baby_steps[giant_step];
                    break;
                }
                ax = qmul(ax, asqm, m);
            }
            if(x == -1) return -1;
            return x + k;
        }
        // g^x == {b0, b1, b2, ...} (mod m)
        std::pair<i64, std::vector<i64>> logs(const std::vector<i64>&b_vec, i64 m, i64 g=-1, i64 phi_m=-1){
            // assert(gcd(bi, m) == 1)
            if(g == -1){
                std::tie(g, phi_m) = std::tuple_cat(primitive_root(m));
            }else if(phi_m == -1){
                phi_m = phi(m);
            }

            auto n = b_vec.size();
            i64 B = std::sqrt(phi_m * n);
            i64 BB = (phi_m + B - 1) / B;

            std::unordered_map<i64, i64> baby_steps;
            for(i64 i=0, gx=1; i<B; i++){
                baby_steps[gx] = i;
                gx = qmul(gx, g, m);
            }
            std::vector<i64> x(n, -1);
            for(int i=0; i<n; i++){
                auto b = b_vec[i];
                if(b == 1 || m == 1) { x[i] = 0; continue; }
                for(i64 j=1, gB=qpow(g, B, m), gx=qmul(gB, inv(b, m), m); j<=BB; j++){
                    if(baby_steps.count(gx)){
                        x[i] = j * B - baby_steps[gx];
                        break;
                    }
                    gx = qmul(gx, gB, m);
                }
            }
            return {g, x};
        }
        // x^a == b (mod m)
        // TODO: b==0 || m not prime
        std::vector<i64> LOG_BSGS(i64 a, i64 b, i64 m) {
            // g^(a*c) == b (mod m)
            // ac == t (mod phi(m))
            auto[g, phi_m] = primitive_root(m);
            if(g == -1) return {};

            auto t = exBSGS(g, b, m);
            if(t == -1) return {};
            auto c = liEu(a, t, phi_m);
            if(c == -1) return {};
            auto delta = phi_m / gcd(a, phi_m);
            std::vector<i64> ans;
            for(auto cur=c%delta; cur<phi_m; cur+=delta){
                ans.push_back(qpow(g, cur, m));
            }
            return ans;
        }
        std::array<std::pair<i64,i64>, 2> fraction(long double x, i64 M = INF, i64 N = INF) {
            i64 m = 1, n = 1, lm=0, ln=1, rm=1, rn=0;
            while(m<=M && n<=N) {
                int k = 0, fl = x*n > m;
                if(fl) { lm = m, ln = n; }
                else { rm = m, rn = n; }
                while(k>=0) {
                    if(fl) { m = lm + (rm<<k), n = ln + (rn<<k); }
                    else { m = (lm<<k) + rm, n = (ln<<k) + rn; }
                    if(m>M || n>N) { break; }
                    if((x*n > m) == fl) {
                        if(fl) { lm = m, ln = n; }
                        else { rm = m, rn = n; }
                        k++;
                    } else { k--; }
                }
                m=lm+rm; n=ln+rn;
            }
            return {std::make_pair(lm, ln), std::make_pair(rm, rn)};
        }
        
        struct Lucas{
            i64 mod;
            std::vector<i64> fac, inv, inv_fac;
            Lucas(i64 _mod) : mod(_mod){
                fac.resize(mod); 
                inv.resize(mod);
                inv_fac.resize(mod);
                fac[1] = inv[1] = inv_fac[1] = 1;
                for(int i=2; i<mod; i++){
                    fac[i] = qmul(fac[i-1], i, mod);
                    inv[i] = qmul(mod - mod / i, inv[mod % i], mod);
                    inv_fac[i] = qmul(inv_fac[i-1], inv[i], mod);
                }
            }
            i64 C(i64 n, i64 m)const{
                if(n <= m || m <= 0){
                    return n == m || m == 0;
                }
                if(n < mod){
                    return qmul(qmul(fac[n], inv_fac[m], mod), inv_fac[n - m], mod);
                }else{
                    return qmul(C(n % mod, m % mod), C(n / mod, m / mod), mod);
                }
            }
            i64 operator()(i64 n, i64 m)const{
                return C(n, m);
            }
        };

        struct exLucas{
            i64 mod;
            std::vector<std::pair<i64, int>> fac;
            std::vector<i64> pk;
            std::vector<std::vector<i64>> prod;
            exLucas(i64 _mod) : mod(_mod), fac(factorize(mod)){
                pk.resize(fac.size());
                prod.resize(fac.size());
                for(int i=0; i<fac.size(); i++){
                    auto[p, k] = fac[i];
                    pk[i] = qpow(p, k);
                    prod[i].resize(pk[i]); prod[i][0] = 1;
                    for(int j=0; j<pk[i]; j+=p){
                        if(j) prod[i][j] = prod[i][j-1];
                        for(int J=1; J<p; J++){
                            prod[i][j + J] = qmul(prod[i][j + J - 1], j + J, pk[i]);
                        }
                    }
                }
            }
            i64 operator()(i64 n, i64 m)const{
                if(n <= m || m <= 0){
                    return n == m || m == 0;
                }
                std::vector<std::pair<i64, i64>> eq(fac.size());
                for(int i=0; i<fac.size(); i++){
                    auto[p, k] = fac[i];
                    auto exponent = [&](i64 n)->i64{
                        i64 res = 0;
                        while(n /= p) res += n;
                        return res;
                    };
                    auto product = [&](i64 n)->i64{
                        i64 res = 1;
                        do{
                            if((n / pk[i]) % 2) res = pk[i] - res;
                            res = qmul(res, prod[i][n % pk[i]], pk[i]);
                        }while(n /= p);
                        return res;
                    };
                    auto e = exponent(n) - exponent(m) - exponent(n - m);
                    if(e >= k){
                        eq[i] = {0, pk[i]};
                    }else{
                        i64 c = qpow(p, e, pk[i], pk[i]-pk[i]/k);
                        c = qmul(c, product(n), pk[i]);
                        c = qmul(c, inv(product(m), pk[i]), pk[i]);
                        c = qmul(c, inv(product(n - m), pk[i]), pk[i]);
                        eq[i] = {c, pk[i]};
                    }
                }
                return exCRT(eq)[0];
            }
        };

        template<std::size_t B>
        struct BSGS{
            i64 m, g, cyc_g;
            std::unordered_map<i64, i64> baby_steps;
            BSGS(i64 _m, i64 _g=-1, i64 _cyc_g=-1) : m(_m), g(_g), cyc_g(_cyc_g){
                // g is generator of m with cyclic of cyc_g
                if(g == -1){
                    std::tie(g, cyc_g) = std::tuple_cat(primitive_root(m));
                }else if(cyc_g == -1){
                    cyc_g = phi(m);
                }

                for(i64 i=0, gx=1; i<B; i++){
                    baby_steps[gx] = i;
                    gx = qmul(gx, g, m);
                }
            }
            // g^x == b (mod m)
            i64 operator()(i64 b)const{
                if(b == 1 || m == 1) return 0;
                i64 BB = (cyc_g + B - 1) / B;

                for(i64 j=1, gB=qpow(g, B, m), gx=qmul(gB, inv(b, m), m); j<=BB; j++){
                    if(baby_steps.count(gx)){
                        return j * B - (baby_steps.find(gx)->second);
                    }
                    gx = qmul(gx, gB, m);
                }
                return -1;
            }
            // a^x == b (mod m)
            // xa*x == xb (mod cyc_g)
            i64 operator()(i64 a, i64 b)const{
                auto&bsgs = *this;
                return liEu(bsgs(a), bsgs(b), cyc_g);
            }
        };

        // based on Index Calculus
        template<std::size_t B>
        struct IndexCalculus{
            i64 mod, g, phi_m;
            std::vector<i64> primes, x_primes;
            IndexCalculus(i64 _mod, i64 _g=-1, i64 _phi_m=-1) : mod(_mod), g(_g), phi_m(_phi_m){
                if(g == -1){
                    std::tie(g, phi_m) = std::tuple_cat(primitive_root(mod));
                }else if(phi_m == -1){
                    phi_m = phi(mod);
                }

                Sieve::init(B);
                auto pn = std::upper_bound(Sieve::primes.begin(), Sieve::primes.end(), B) - Sieve::primes.begin() - 1;
                primes = {Sieve::primes.begin(), Sieve::primes.begin() + pn + 1};
                x_primes = logs(primes, mod, g, phi_m).second;
            }
            // g^x == b (mod m);
            i64 operator()(i64 b)const{
                static std::mt19937 eng{std::random_device{}()};
                for(i64 y = 0; ; y = eng() % phi_m){
                    auto z = qmul(b, qpow(g, y, mod, phi_m), mod);
                    auto max_p = max_prime_factor(z);
                    if(max_p > B) continue;
                    auto fac = factorize(z / max_p);
                    sort(fac.begin(), fac.end());
                    if(fac.empty() || fac.back().first != max_p){
                        if(max_p > 1) fac.emplace_back(max_p, 1);
                    }else fac.back().second++;
                    i64 ans = 0;
                    for(int i=0, j=0; i<fac.size(); i++){
                        auto[p, k] = fac[i];
                        while(primes[j] < p) j++;
                        ans = (ans + qmul(k, x_primes[j], phi_m)) % phi_m;
                    }
                    return (ans + phi_m - y) % phi_m;
                }
            }
        };

        // based on Pohlig–Hellman
        template<std::size_t B>
        struct PohligHellman{
            i64 m, g, phi_m;
            std::vector<std::pair<i64, int>> fac;
            std::vector<i64> pk, gi;
            std::vector<BSGS<B>> bsgs;
            PohligHellman(i64 _m, i64 _g=-1, i64 _phi_m=-1) : m(_m), g(_g), phi_m(_phi_m){
                if(g == -1){
                    std::tie(g, phi_m) = std::tuple_cat(primitive_root(m));
                }else if(phi_m == -1){
                    phi_m = phi(m);
                }

                fac = factorize(phi_m);
                pk.resize(fac.size());
                gi.resize(fac.size());
                bsgs.reserve(fac.size());
                for(int i=0; i<fac.size(); i++){
                    auto[p, k] = fac[i];
                    pk[i] = qpow(p, k);
                    gi[i] = qpow(g, phi_m / pk[i], m);
                    bsgs.emplace_back(m, gi[i], pk[i]);
                }
            }
            // g^x == h (mod m);
            i64 operator()(i64 h)const{
                std::vector<std::pair<i64, i64>> eq;
                for(int i=0; i<fac.size(); i++){
                    auto[p, k] = fac[i];
                    auto hi = qpow(h, phi_m / pk[i], m);
                    auto calc = [&](i64 gi, i64 hi)->i64{
                        i64 x = 0, pe = 1, ga = qpow(gi, pk[i]/p, m);
                        for(int e=0; e<k; e++, pe*=p){
                            auto he = qpow(qmul(hi, inv(qpow(gi, x, m), m), m), pk[i]/pe/p, m);
                            x = (x + qmul(pe, bsgs[i](ga, he), pk[i])) % pk[i];
                        }
                        return x;
                    };
                    auto xi = calc(gi[i], hi);
                    eq.emplace_back(xi, pk[i]);
                }
                return exCRT(eq)[0];
            }
            i64 operator()(i64 a, i64 b)const{
                auto&bsgs = *this;
                return liEu(bsgs(a), bsgs(b), phi_m);
            }
        };
    }
}

#endif
```

## including

+ Linear Sieve(Euler Sieve)
+ Eratosthenes Sieve
+ Linear Preprocessing GCD Algorithm
+ Greatest Common Divisor(Euclid's Algorithm、Euclidean Algorithm、Binary GCD Algorithm)
+ Extended Euclidean Algorithm
+ Binary Exponentiation
+ Miller Rabin Primality Test
+ Pollard Rho Algorithm
+ Factorization Algorithm
+ Fermat's Little Theorem
+ Primitive Root Algorithm
+ Chinese Remainder Theorem(extended CRT Algorithm)
+ Linear Congruence Equation Algorithm
+ Baby-step Giant-step Algorithm(extended BSGS Algorithm、Amortized BSGS Data Structure)
+ Stern Brocot Tree Algorithm
+ Lucas's Theorem(Wilson's Theorem、Kummer Theory、extended Lucas Algorithm)
+ Discrete Logarithms Algorithm(Index Calculus Algorithm、Pohlig-Hellman Algorithm)

# Matrix of fixed size Data Structure

```cpp
#ifndef MATRIX
#define MATRIX
// 矩阵

template<std::size_t N, std::size_t M, typename I=long long>
struct Matrix {
    I val[N][M];
    
    static Matrix IE;
    
    I* operator[](int i) {
        return val[i];
    }
    const I* operator[](int i) const {
        return val[i];
    }
    Matrix<N,M,I> operator+(const Matrix<N,M,I>&o)const {
        Matrix<N,M,I> res;
        for(int i=0; i<N; i++) {
            for(int j=0; j<M; j++) {
                res[i][j] = val[i][j] + o[i][j];
            }
        }
        return res;
    }
    template<std::size_t K>
    Matrix<N,K,I> operator*(const Matrix<M,K,I>&o)const {
        Matrix<N,K,I> res;
        for(int i=0; i<N; i++) {
            for(int k=0; k<K; k++) {
                res[i][k] = I();
                for(int j=0; j<M; j++) {
                    res[i][k] += val[i][j] * o[j][k];
                }
            }
        }
        return res;
    }
    friend std::ostream& operator<<(std::ostream&out, const Matrix<N,M,I>&x) {
        for(int i=0; i<N; i++) {
            for(int j=0; j<M; j++) {
                out<<x[i][j]<<" \n"[j+1==M];
            }
        }
        return out;
    }
};

template<std::size_t N, std::size_t M, typename I>
Matrix<N,M,I> Matrix<N,M,I>::IE = []() {
    Matrix<N,M,I> IE;
    for(int i=0; i<N; i++) {
        for(int j=0; j<M; j++) {
            IE[i][j] = i==j ? I(1) : I(0);
        }
    }
    return IE;
}();

#endif
```

# Mo's Algorithm

```cpp
#ifndef MO_ALG
#define MO_ALG

struct Info{
    int ans = 0;
    std::vector<int> info, cnt;
    Info(int n) : info(n), cnt(1e6+1, 0){}
};

struct MO : public Info{
    int N = 1, M = 1, B;
    int l, r, t, id;

    struct Query{
        int l, r, t, id;
    };
    std::vector<Query> query;

    struct Change{
        int x, v;
    };
    std::vector<Change> change;


    MO(int n) : Info(n){};

    void add_query(Query&&q){
        q.t = change.size()-1;
        q.id = query.size();
        query.emplace_back(q);
    }

    void add_change(Change&&c){
        change.emplace_back(c);
    }

    void init(){
        B = std::pow(std::pow(info.size(), 2) * (change.size()+1) / (query.size()+1), 1./3.);
        std::sort(query.begin(), query.end(), [this](auto&a,auto&b){
            return std::make_tuple(a.l/B, a.r/B, a.t/B) < 
                std::make_tuple(b.l/B, b.r/B, b.t/B);
        });
        r = t = -1, l = id = 0;
    }

    void redo(int i){
        if(0==cnt[info[i]]++){
            ans++;
        }
    }

    void undo(int i){
        if(0==--cnt[info[i]]){
            ans--;
        }
    }

    void update(int i){
        auto&[x, v] = change[i];
        if(l<=x && x<=r) undo(x);
        std::swap(v, info[x]);
        if(l<=x && x<=r) redo(x);
    }

    Query next(){
        auto x = query[id++];
        while(l < x.l) undo(l++);
        while(l > x.l) redo(--l);
        while(r < x.r) redo(++r);
        while(r > x.r) undo(r--);
        while(t < x.t) update(++t);
        while(t > x.t) update(t--);
        return x;
    }

    bool finish(){
        return id == query.size();
    }
};

#endif
```

# Modular Arithmetic

```cpp
#ifndef MODINT_CPP
#define MODINT_CPP
// 模数循环群

using i64 = long long;
 
template<class T>
constexpr T qpow(T a,i64 b){
    if(0==a)return a;
    T res = 1;
    while(b){
        if(b&1) res *= a;
        a *= a;
        b >>= 1;
    }
    return res;
}
constexpr i64 mul(i64 a,i64 b,i64 p){
    i64 res = a*b - p*i64(1.L*a*b/p);
    res %= p;
    if(res<0) res += p;
    return res;
}
template<i64 P=0LL>
struct MLong{
    i64 x;
    constexpr MLong():x(0){}
    constexpr MLong(i64 x):x(norm(x%getMod())){}
    
    static i64 Mod;
    constexpr static i64 getMod(){
        return P>0?P:Mod;
    }
    constexpr static void setMod(i64 Mod_){
        Mod = Mod_;
    }
    constexpr i64 norm(i64 x)const{
        if(x<0){
            x += getMod();
        }
        if(x>=getMod()){
            x -= getMod();
        }
        return x;
    }
    constexpr i64 val()const{
        return x;
    }
    explicit constexpr operator i64()const{
        return x;
    }
    constexpr MLong operator-()const{
        MLong res;
        res.x = norm(getMod() - x);
        return res;
    }
    constexpr MLong inv()const{
        assert(x != 0);
        return qpow(*this,getMod() - 2);
    }
    constexpr MLong &operator*=(MLong rhs)&{
        x = mul(x,rhs.x,getMod());
        return *this;
    }
    constexpr MLong &operator+=(MLong rhs)&{
        x = norm(x + rhs.x);
        return *this;
    }
    constexpr MLong &operator-=(MLong rhs)&{
        x = norm(x - rhs.x);
        return *this;
    }
    constexpr MLong &operator/=(MLong rhs)&{
        return *this *= rhs.inv();
    }
    friend constexpr MLong operator*(MLong lhs,MLong rhs){
        MLong res = lhs;
        res *= rhs;
        return res;
    }
    friend constexpr MLong operator+(MLong lhs,MLong rhs){
        MLong res = lhs;
        res += rhs;
        return res;
    }
    friend constexpr MLong operator-(MLong lhs,MLong rhs){
        MLong res = lhs;
        res -= rhs;
        return res;
    }
    friend constexpr MLong operator/(MLong lhs,MLong rhs){
        MLong res = lhs;
        res /= rhs;
        return res;
    }
    friend constexpr istream &operator>>(istream &is,MLong &a){
        i64 v;
        is >> v;
        a = MLong(v);
        return is;
    }
    friend constexpr ostream &operator<<(ostream &os,const MLong &a){
        return os << a.val();
    }
    friend constexpr bool operator==(MLong lhs,MLong rhs){
        return lhs.val() == rhs.val();
    }
    friend constexpr bool operator!=(MLong lhs,MLong rhs){
        return !(lhs == rhs);
    }
    friend constexpr bool operator<(MLong lhs,MLong rhs){
        return lhs.val() < rhs.val();
    }
    friend constexpr bool operator<=(MLong lhs,MLong rhs){
        return lhs < rhs || lhs == rhs;
    }
    friend constexpr bool operator>(MLong lhs,MLong rhs){
        return !(lhs <= rhs);
    }
    friend constexpr bool operator>=(MLong lhs,MLong rhs){
        return !(lhs < rhs);
    }
};
 
template<>
i64 MLong<0LL>::Mod = 1;
 
template<int P=0>
struct MInt{
    int x;
    constexpr MInt():x(0){}
    constexpr MInt(i64 x):x(norm(x%getMod())){}
    
    static int Mod;
    constexpr static int getMod(){
        return P>0?P:Mod;
    }
    constexpr static void setMod(int Mod_){
        Mod = Mod_;
    }
    constexpr int norm(int x)const{
        if(x<0){
            x += getMod();
        }
        if(x>=getMod()){
            x -= getMod();
        }
        return x;
    }
    constexpr int val()const{
        return x;
    }
    explicit constexpr operator int()const{
        return x;
    }
    constexpr MInt operator-()const{
        MInt res;
        res.x = norm(getMod() - x);
        return res;
    }
    constexpr MInt inv()const{
        assert(x != 0);
        return qpow(*this,getMod() - 2);
    }
    constexpr MInt &operator*=(MInt rhs)&{
        x = 1LL * x * rhs.x % getMod();
        return *this;
    }
    constexpr MInt &operator+=(MInt rhs)&{
        x = norm(x + rhs.x);
        return *this;
    }
    constexpr MInt &operator-=(MInt rhs)&{
        x = norm(x - rhs.x);
        return *this;
    }
    constexpr MInt &operator/=(MInt rhs)&{
        return *this *= rhs.inv();
    }
    friend constexpr MInt operator*(MInt lhs,MInt rhs){
        MInt res = lhs;
        res *= rhs;
        return res;
    }
    friend constexpr MInt operator+(MInt lhs,MInt rhs){
        MInt res = lhs;
        res += rhs;
        return res;
    }
    friend constexpr MInt operator-(MInt lhs,MInt rhs){
        MInt res = lhs;
        res -= rhs;
        return res;
    }
    friend constexpr MInt operator/(MInt lhs,MInt rhs){
        MInt res = lhs;
        res /= rhs;
        return res;
    }
    friend constexpr istream &operator>>(istream &is,MInt &a){
        i64 v;
        is >> v;
        a = MInt(v);
        return is;
    }
    friend constexpr ostream &operator<<(ostream &os,const MInt &a){
        return os << a.val();
    }
    friend constexpr bool operator==(MInt lhs,MInt rhs){
        return lhs.val() == rhs.val();
    }
    friend constexpr bool operator!=(MInt lhs,MInt rhs){
        return lhs.val() != rhs.val();
    }
    friend constexpr bool operator<(MInt lhs,MInt rhs){
        return lhs.val() < rhs.val();
    }
    friend constexpr bool operator<=(MInt lhs,MInt rhs){
        return lhs < rhs || lhs == rhs;
    }
    friend constexpr bool operator>(MInt lhs,MInt rhs){
        return !(lhs <= rhs);
    }
    friend constexpr bool operator>=(MInt lhs,MInt rhs){
        return !(lhs < rhs);
    }
};
 
template<>
int MInt<0>::Mod = 1;
 
template<int V,int P>
constexpr MInt<P> CInv = MInt<P>(V).inv();
 
constexpr int P = 1000000007;
using Z = MInt<P>;

#endif
```

# Policy-Based Data Structures

```cpp
#ifndef STL_PBDS
#define STL_PBDS

#define _EXT_CODECVT_SPECIALIZATIONS_H 1
#define _EXT_ENC_FILEBUF_H 1
#include <bits/extc++.h>
// using namespace __gnu_cxx;
// using namespace __gnu_pbds;

// 平衡树
template<typename K_T,
    typename Cmp = std::less<K_T>,
    typename V_T = __gnu_pbds::null_type
>
using order_set = __gnu_pbds::tree<
    K_T,                                                // key_type
    V_T,                                                // value_type
    Cmp,                                                // comparator
    __gnu_pbds::rb_tree_tag,                            // tag
    __gnu_pbds::tree_order_statistics_node_update       // policy
>;
template<typename K_T,
    typename V_T = __gnu_pbds::null_type,
    typename Cmp = std::less<K_T>
>
using order_map = order_set<K_T, Cmp, V_T>;
template<typename K_T, typename V_T, typename Cmp>
V_T& operator%(order_map<K_T,V_T,Cmp>&mp, const K_T&x){
    return mp.lower_bound(x)->second;
}
/*
tag:
rb_tree_tag                                 红黑树
splay_tree_tag                              Slpay
ov_tree_tag                                 向量树

Itr ::point_iterator

std::pair<point_iterator, bool> insert(T)   插入
bool erase(T/Itr)                           删除元素/迭代器
int order_of_key(T)                         返回排名
Itr find_by_order(T)                        排名对应元素
Itr lower_bound(T)/upper_bound(T)           二分查找
void join(order_set)                        合并
void split(T,order_set)                     保留<=,其余分离覆盖到order_set中
bool empty()                                判空
size_t size()                               大小
Itr begin()/end()                           首尾迭代器
*/

/***************/

// 堆
template<typename T,
    typename Cmp = std::greater<T>
>
using heap = __gnu_pbds::priority_queue<
    T,                                                  // type
    Cmp,                                                // comparator
    __gnu_pbds::pairing_heap_tag                        // tag
>;
/*
tag:
pairing_heap_tag        配对堆
thin_heap_tag           斐波那契堆
binomial_heap_tag       二项堆
binary_heap_tag         二叉堆

Itr ::point_iterator    可以指定为nullptr

usage:
Itr push(T)             入堆
void pop()              出堆
T top()                 堆顶
void modify(Itr, T)     修改
void join(heap)         合并,清空heap
bool empty()            判空
size_t size()           大小
void clear()            清空
*/

/***************/

// 哈希表
const int RANDOM = std::chrono::high_resolution_clock::now().time_since_epoch().count();
template<typename K_T>
struct Chash{
    static std::hash<K_T> hash;
    int operator()(K_T x)const{return hash(x)^RANDOM;}
};
template<typename K_T, typename V_T, typename Hash = Chash<K_T>>
using hash_table = __gnu_pbds::cc_hash_table<K_T, V_T, Hash>;

/*
tag:
cc_hash_table           拉链法
gp_hash_table           二次探测法

V_T& operaotr[](K_T)    映射
*/


// 字典树
using trie = __gnu_pbds::trie<
    std::string,                                    //
    __gnu_pbds::null_type,                          //
    __gnu_pbds::trie_string_access_traits<>,        //
    __gnu_pbds::pat_trie_tag,                       // tag
    __gnu_pbds::trie_prefix_search_node_update      // policy
>;
/*
Itr insert(string)                          插入
void erase(string)                          删除
void join(trie)                             合并trie
std::pair<Itr, Itr> prefix_range(string)    前缀遍历[beign,end)
*/

#endif
```


# Randomization

```cpp
#ifndef RANDOM
#define RANDOM
// 随机

namespace Random{
    std::mt19937 eng(std::random_device{}());
    double double_uniform(double l,double r){
        std::uniform_real_distribution<> dis(l,r);
        return dis(eng);
    }
    int int_uniform(int l,int r){
        std::uniform_int_distribution<> dis(l,r);
        return dis(eng);
    }
};

#endif
```

# Strongly Connected Component Algorithm

```cpp
#ifndef STRONGLY_CONNECTED_COMPONENT
#define STRONGLY_CONNECTED_COMPONENT

struct SCC {
    int n, cur, cnt;
    std::vector<std::vector<int>> adj;
    std::vector<int> stk;
    std::vector<int> dfn, low, bel;
    
    SCC() = default;
    SCC(int n) {init(n);}
    
    void init(int n) {
        this->n = n;
        adj.assign(n, {});
        dfn.assign(n, -1);
        low.resize(n);
        bel.assign(n, -1);
        stk.clear();
        cur = cnt = 0;
    }
    
    void add_edge(int u, int v) {
        adj[u].push_back(v);
    }
    
    void dfs(int u) {
        dfn[u] = low[u] = cur++;
        stk.push_back(u);
        
        for(auto v : adj[u]) {
            if(dfn[v] == -1) {
                dfs(v);
                low[u] = std::min(low[u], low[v]);
            } else if(bel[v] == -1) {
                low[u] = std::min(low[u], dfn[v]);
            }
        }
        
        if(dfn[u] == low[u]) {
            for(int v=-1; v!=u; stk.pop_back()) {
                v = stk.back();
                bel[v] = cnt;
            }
            cnt++;
        }
    }
    
    std::vector<int> work() {
        for(int i=0; i<n; i++) {
            if(dfn[i] != -1) { continue; }
            dfs(i);
        }
        return bel;
    }
};

#endif
```

# ZKW Segment Tree Data Structure

```cpp
#ifndef SEGTREE_ZKW_HPP
#define SEGTREE_ZKW_HPP

struct Info{
    int g = 0;

    Info operator+(const Info&o){
        return {__gcd(g, o.g)};
    }
};

template<class Info>
struct SegTreeZKW{
    int n, N;
    std::vector<Info> info;

    SegTreeZKW() = default;
    template<typename ...Args>
    explicit SegTreeZKW(Args&&... args){ init(std::forward<Args>(args)...); }

    void init(int n){
        this->n = n;
        this->N = 1 << (__lg(n) + 1);
        info.assign(N << 1, Info{});
    }

    void init(const std::vector<Info>&arr){
        init(arr.size());
        std::copy(arr.begin(), arr.end(), info.begin() + N);
        for(int x=N - 1; x; x--){
            info[x] = info[x * 2] + info[x * 2 + 1];
        }
    }

    void apply(int x, const Info&v){
        info[x += N] = v;
        while(x /= 2){
            info[x] = info[x * 2] + info[x * 2 + 1];
        }
    }

    Info ask(int l, int r){
        if(l > r) return {};
        l += N, r += N;
        Info infol = info[l], infor = l == r ? Info{} : info[r];
        for(; l + 1 < r; l /= 2, r /= 2){
            if(~l & 1) infol = infol + info[l ^ 1];
            if(r & 1) infor = info[r ^ 1] + infor;
        }
        return infol + infor;
    }
};

#endif
```

# ZKW Segment Tree Data Structure with Lazy Tag

```cpp
#ifndef SEGTREE_ZKW_TAG_HPP
#define SEGTREE_ZKW_TAG_HPP

struct Tag{
    char ch{-1};

    void operator+=(const Tag&t){
        if(t.ch < 0) return;
        ch = t.ch;
    }
};

struct Info{
    std::array<int,26> cnt {};
    int size = 0;

    Info operator+(const Info&o){
        Info res;
        for(int i=0;i<26;i++)res.cnt[i] = cnt[i] + o.cnt[i];
        res.size = size + o.size;
        return res;
    }
    void operator+=(const Tag&t){
        if(t.ch < 0) return;
        cnt.fill(0);
        cnt[t.ch - 'a'] = size;
    }
};

template<class Info, class Tag>
struct SegTreeZKW{
    int n, lg, N;
    std::vector<Info> info;
    std::vector<Tag> tag;

    SegTreeZKW() = default;
    template<typename ...Args>
    explicit SegTreeZKW(Args&&... args){ init(std::forward<Args>(args)...); }

    void init(int n){
        this->n = n;
        this->lg = std::__lg(n) + 1;
        this->N = 1 << lg;
        info.assign(N << 1, Info{});
        tag.assign(N << 1, Tag{});
    }

    void init(const std::vector<Info>&arr){
        init(arr.size());
        std::copy(arr.begin(), arr.end(), info.begin() + N);
        for(int x=N - 1; x; x--){
            info[x] = info[x * 2] + info[x * 2 + 1];
        }
    }

    void push_down(int x) {
        info[x * 2] += tag[x]; tag[x * 2] += tag[x];
        info[x * 2 + 1] += tag[x]; tag[x * 2 + 1] += tag[x];
        tag[x] = Tag{};
    }

    void apply(int l, int r, const Tag&t){
        if(l > r) return;
        l += N, r += N;
        for(auto i=lg; i; i--) push_down(l >> i), push_down(r >> i);
        info[l] += t;
        if(l < r){
            info[r] += t;
            for(; l ^ r ^ 1; ){
                if(~l & 1) info[l ^ 1] += t, tag[l ^ 1] += t;
                if( r & 1) info[r ^ 1] += t, tag[r ^ 1] += t;
                l /= 2, r /= 2;
                info[l] = info[l * 2] + info[l * 2 + 1];
                info[r] = info[r * 2] + info[r * 2 + 1];
            }
        }
        for(; l /= 2; ) info[l] = info[l * 2] + info[l * 2 + 1];
    }

    Info ask(int l, int r){
        if(l > r) return {};
        l += N, r += N;
        for(auto i=lg; i; i--) push_down(l >> i), push_down(r >> i);
        Info infol = info[l];
        if(l < r){
            Info infor = l == r ? Info{} : info[r];
            for(; l ^ r ^ 1; ){
                if(~l & 1) infol = infol + info[l ^ 1];
                if( r & 1) infor = info[r ^ 1] + infor;
                l /= 2, r /= 2;
            }
            infol = infol + infor;
        }
        return infol;
    }
};

#endif
```

# Sparse Table Data Structure

```cpp
#ifndef SEPARATE_TABLE
#define SEPARATE_TABLE

template<typename T>
struct ST{
    unsigned N, B;
    T IE;
    std::function<T(T,T)> calc;
    std::vector<std::vector<T>> a, b;
    std::vector<T> pre, suf;

    ST(){}
    template<typename ... Args>
    ST(Args ... args){
        init(args ...);
    }
    template<typename ... Args>
    T Calc(const T&x, const T&y, Args ... args){
        return Calc(calc(x, y), args...);
    }
    T Calc(const T&x, const T&y){return calc(x, y);}
    void init(const std::vector<T>&v, 
        const std::function<T(T,T)>&calc = [](T x,T y){return std::max(x, y);},
        const T&IE = T()){
        this->N = v.size();
        if(N <= 0) return;
        this->B = sqrt(N) + 1;
        this->calc = calc;
        this->IE = IE;
        pre = suf = v;
        const int M = (N + B - 1) / B;
        const int lgM = std::__lg(M);
        const int lgB = std::__lg(B);
        a.resize(lgM+1); a[0].resize(M);
        b.resize(lgB+1); b[0] = v;
        for(int i=0; i<M; i++){
            a[0][i] = v[i * B];
            for(int j=1; j<B && i*B + j < N; j++){
                a[0][i] = Calc(a[0][i], v[i*B + j]);
            }
        }
        for(int j=0; j<lgM; j++){
            a[j+1].resize(M - (2<<j) + 1);
            for(int i=0; i+(2<<j)<=M; i++){
                a[j+1][i] = Calc(a[j][i], a[j][i + (1<<j)]);
            }
        }
        for(int j=0; j<lgB && (2<<j)<=N; j++){
            b[j+1].resize(N - (2<<j) + 1);
            for(int i=0; i+(2<<j)<=N; i++){
                b[j+1][i] = b[j][i];
                if((i+(1<<j))/B == i/B){
                    b[j+1][i] = Calc(b[j+1][i], b[j][i + (1<<j)]);
                }
            }
        }
        for(int i=1; i<N; i++){
            if(i%B != 0){
                pre[i] = Calc(pre[i-1], pre[i]);
            }
        }
        for(int i=N-2; i>=0; i--){
            if(i%B != B-1){
                suf[i] = Calc(suf[i], suf[i+1]);
            }
        }
    }
    T ask_O1(int l, int r){
        if(l > r) return IE;
        assert(0<=l && r < N);
        if(l/B != r/B){
            T ans = Calc(suf[l], pre[r]);
            l = l/B + 1;
            r = r/B - 1;
            if(l <= r){
                int k = std::__lg(r-l+1);
                ans = Calc(ans, a[k][l], a[k][r + 1 - (1<<k)]);
            }
            return ans;
        }else{
            int k = std::__lg(r-l+1);
            return Calc(b[k][l], b[k][r + 1 - (1<<k)]);
        }
    }
    T ask(int l, int r){
        if(l > r) return IE;
        assert(0<=l && r < N);
        T ans = IE;
        if(l/B != r/B){
            ans = suf[l];
            T ansr = pre[r];
            l = l/B + 1; r = r/B - 1;
            while(l <= r){
                int k = std::__lg(r-l+1);
                ans = Calc(ans, a[k][l]);
                l += (1<<k);
            }
            ans = Calc(ans, ansr);
        }else{
            while(l <= r){
                int k = std::__lg(r-l+1);
                ans = Calc(ans, b[k][l]);
                l += (1<<k);
            }
        }
        return ans;
    }
};

#endif
```

# String Algorithm

```cpp
#ifndef STRING
#define STRING

namespace String {
    std::vector<int> manacher(const std::string&s) {
        std::string t = "#";
        for(auto c : s) {
            t += c; t += '#';
        }
        int n = t.size();
        std::vector<int> r(n);
        for(int i=0, j=0; i<n; i++) {
            if(2*j-i >= 0 && j+r[j] > i) {
                r[i] = std::min(r[2*j-i], j+r[j]-i);
            }
            while(i-r[i] >= 0 && i+r[i] < n && t[i-r[i]] == t[i+r[i]]) {
                r[i]++;
            }
            if(i+r[i] > j+r[j]) { j = i; }
        }
        return r;
    }
    
    std::vector<int> exkmp(const std::string&s) {
        int n = s.size();
        std::vector<int> z(n+1, 0); z[0] = n;
        for(int i=1, j=1; i<n; i++) {
            z[i] = std::max(0, std::min(j+z[j]-i, z[i-j]));
            while(i+z[i] < n && s[z[i]] == s[i+z[i]]) {
                z[i]++;
            }
            if(i+z[i] > j+z[j]) { j = i; }
        }
        return z;
    }
    std::vector<int> exkmp(const std::string&s, const std::string&t) {
        auto z = exkmp(s + t);
        return {z.begin()+s.size(), z.end()};
    }
    
    struct SA {
        int n;
        std::vector<int> sa, rk, lc;
        SA(const std::string &s) {
            n = s.length();
            sa.resize(n);
            lc.resize(n - 1);
            rk.resize(n);
            std::iota(sa.begin(), sa.end(), 0);
            std::sort(sa.begin(), sa.end(), [&](int a, int b) {return s[a] < s[b];});
            rk[sa[0]] = 0;
            for(int i=1; i<n; i++) {
                rk[sa[i]] = rk[sa[i-1]] + (s[sa[i]] != s[sa[i-1]]);
            }
            int k = 1;
            std::vector<int> tmp, cnt(n);
            tmp.reserve(n);
            while(rk[sa[n-1]] < n-1) {
                tmp.clear();
                for(int i=0; i<k; i++) { tmp.push_back(n - k + i); }
                for(auto i : sa) if(i>=k) { tmp.push_back(i - k); }
                std::fill(cnt.begin(), cnt.end(), 0);
                for(int i=0; i<n; i++) { ++cnt[rk[i]]; }
                for(int i=1; i<n; i++) { cnt[i] += cnt[i-1]; }
                for(int i=n-1; i>=0; i--) { sa[--cnt[rk[tmp[i]]]] = tmp[i]; }
                std::swap(rk, tmp);
                rk[sa[0]] = 0;
                for(int i=1; i<n; i++) {
                    rk[sa[i]] = rk[sa[i-1]] + (tmp[sa[i-1]] < tmp[sa[i]] || sa[i-1] + k == n || tmp[sa[i-1] + k] < tmp[sa[i] + k]);
                }
                k *= 2;
            }
            for(int i=0, j=0; i < n; ++i) {
                if(rk[i] == 0) {
                    j = 0;
                } else {
                    for(j -= j > 0; i+j < n && sa[rk[i]-1]+j < n && s[i+j] == s[sa[rk[i]-1] + j]; ++j);
                    lc[rk[i]-1] = j;
                }
            }
        }
    };
}

#endif
```

## including

+ Manacher's Algorithm
+ Z Algorithm(extended KMP Algorithm)
+ Suffix Array Data Structure

# String Hash Data Structure

```cpp
#ifndef STRING_HASH
#define STRING_HASH

using u64 = unsigned long long;
struct SHash{
    static constexpr std::size_t C = 2;
    static constexpr u64 M[] = {
        (-1u + 1ull),1000000007,1118872217,122420729,163227661,
        217636919,290182597,386910137,515880193,687840301,
        917120411,1222827239,1610612741,3221225473ul,4294967291ul
    };

    std::array<u64, C> val;

    SHash() { val.fill(0ull); }

    SHash(const std::initializer_list<u64>&list) {
        std::copy(list.begin(), list.begin()+C, val.begin());
    }

    bool operator==(const SHash&o)const{
        return val == o.val;
    }
    bool operator<(const SHash&o)const{
        return val < o.val;
    }
    SHash operator*(const SHash&o)const{
        SHash res = *this;
        for(int i=0; i<C; i++){
            res.val[i] *= o.val[i];
            res.val[i] %= M[i];
        }
        return res;
    }
    SHash operator+(u64 y)const{
        SHash res = *this;
        for(auto&x : res.val) x += y;
        return res.norm();
    }
    SHash operator+(const SHash&o)const{
        SHash res;
        std::transform(val.begin(), val.end(), o.val.begin(), res.val.begin(), std::plus());
        return res.norm();
    }
    SHash operator-(u64 y)const{
        SHash res = *this;
        for(int i=0; i<C; i++){
            if(res.val[i] >= y) res.val[i] -= y;
            else res.val[i] += M[i] - y;
        }
        return res;
    }
    SHash operator-(const SHash&o)const{
        SHash res;
        for(int i=0; i<C; i++){
            if(val[i] >= o.val[i]) res.val[i] = val[i] - o.val[i];
            else res.val[i] = val[i] + M[i] - o.val[i];
        }
        return res.norm();
    }
    SHash norm(){
        for(int i=0; i<C; i++){
            if(val[i] >= M[i]) val[i] -= M[i];
        }
        return *this;
    }
};
namespace std {
    template <>
    struct hash<SHash> {
        std::size_t operator()(SHash x)const{
            return std::accumulate(x.val.begin(), x.val.end(), 0, std::bit_xor());
        }
    };
}
static const SHash B = {37, 53, 71, 97, 137, 251, 353, 491, 599, 617, 773, 853, 977, 1009};
constexpr int MAXN = 1000005;
static const std::array<SHash, MAXN+1> Bp = [](){
    std::array<SHash, MAXN+1> Bp;
    Bp[0].val.fill(1ull);
    for(int i=1; i<=MAXN; i++){
        Bp[i] = Bp[i-1] * B;
    }
    return Bp;
}();

std::vector<SHash> build_prefix(std::string s){
    std::vector<SHash> prefix(s.size() + 1);
    for(int i=1; i<=s.size(); i++){
        prefix[i] = prefix[i-1]*B + s[i-1];
    }
    return prefix;
}

SHash get(const std::vector<SHash>&prefix, int l, int r){
    if(l > r) return SHash{};
    return prefix[r] - prefix[l-1];
}

#endif
```

# Simpson Algorithm

自适应辛普森算法（定积分）

```cpp
const double Pi = std::acos(-1.0);
constexpr double EPS = 1e-9;
double v, r, d;
double f(double x) {
    double s = std::sin(x);
    return 1 / v / (std::sqrt(s * s + 3) - s);
}
double simpson(double l, double r) {
    return (f(l) + 4 * f((l + r) / 2) + f(r)) * (r - l) / 6;
}
double integral(double l, double r, double eps, double st) {
    double mid = (l + r) / 2;
    double sl = simpson(l, mid);
    double sr = simpson(mid, r);
    if (std::abs(sl + sr - st) <= 15 * eps)
        return sl + sr + (sl + sr - st) / 15;
    return integral(l, mid, eps / 2, sl) + integral(mid, r, eps / 2, sr);
}
double integral(double l, double r) {
    return integral(l, r, EPS, simpson(l, r));
}
```

$$
\text{Catalan}[n]=\frac{C(2n,n)}{n+1}
$$


# Blossom Algorithm

```cpp
struct Graph {
    int n;
    std::vector<std::vector<int>> e;
    Graph(int n) : n(n), e(n) {}
    void addEdge(int u, int v) {
        e[u].push_back(v);
        e[v].push_back(u);
    }
    std::vector<int> findMatching(int m, const auto &init) {
        std::vector<int> match(n, -1), vis(n), link(n), f(n), dep(n);
        for (auto [x, y] : init) {
            match[x] = y;
            match[y] = x;
        }
        // disjoint set union
        auto find = [&](int u) {
            while (f[u] != u)
                u = f[u] = f[f[u]];
            return u;
        };
        auto lca = [&](int u, int v) {
            u = find(u);
            v = find(v);
            while (u != v) {
                if (dep[u] < dep[v])
                    std::swap(u, v);
                u = find(link[match[u]]);
            }
            return u;
        };
        std::queue<int> que;
        auto blossom = [&](int u, int v, int p) {
            while (find(u) != p) {
                link[u] = v;
                v = match[u];
                if (vis[v] == 0) {
                    vis[v] = 1;
                    que.push(v);
                }
                f[u] = f[v] = p;
                u = link[v];
            }
        };
        // find an augmenting path starting from u and augment (if exist)
        auto augment = [&](int u) {
            while (!que.empty())
                que.pop();
            std::iota(f.begin(), f.end(), 0);
            // vis = 0 corresponds to inner vertices, vis = 1 corresponds to outer vertices
            std::fill(vis.begin(), vis.end(), -1);
            que.push(u);
            vis[u] = 1;
            dep[u] = 0;
            int y = -1;
            while (!que.empty()){
                int u = que.front();
                que.pop();
                if (u >= m) {
                    y = u;
                }
                for (auto v : e[u]) {
                    if (vis[v] == -1) {
                        vis[v] = 0;
                        link[v] = u;
                        dep[v] = dep[u] + 1;
                        // found an augmenting path
                        if (match[v] == -1) {
                            for (int x = v, y = u, temp; y != -1; x = temp, y = x == -1 ? -1 : link[x]) {
                                temp = match[y];
                                match[x] = y;
                                match[y] = x;
                            }
                            return;
                        }
                        vis[match[v]] = 1;
                        dep[match[v]] = dep[u] + 2;
                        que.push(match[v]);
                    } else if (vis[v] == 1 && find(v) != find(u)) {
                        // found a blossom
                        int p = lca(u, v);
                        blossom(u, v, p);
                        blossom(v, u, p);
                    }
                }
            }
            if (y != -1) {
                for (int x = -1, temp; y != -1; x = temp, y = x == -1 ? -1 : link[x]) {
                    temp = match[y];
                    if (x != -1) {
                        match[x] = y;
                    }
                    match[y] = x;
                }
            }
        };
        // find a maximal matching greedily (decrease constant)
        // auto greedy = [&]() {
        //     for (int u = 0; u < n; ++u) {
        //         if (match[u] != -1)
        //             continue;
        //         for (auto v : e[u]) {
        //             if (match[v] == -1) {
        //                 match[u] = v;
        //                 match[v] = u;
        //                 break;
        //             }
        //         }
        //     }
        // };
        // greedy();
        for (int u = 0; u < m; ++u)
            if (match[u] == -1)
                augment(u);
        return match;
    }
};
```

# Splay Tree Data Structure

```cpp
namespace Splay{
    constexpr int D = 27;
    struct Info {
        int up[D][2] {};
        int down[D][2] {};
        int t = 0;
        i64 ans = 0;
    };

    Info operator+(const Info &a, const Info &b) {
        Info c;
        c.t = a.t ^ b.t;
        c.ans = a.ans + b.ans;
        for (int i = 0; i < D; i++) {
            for (int j = 0; j < 2; j++) {
                c.ans += (1LL << i) * a.down[i][j] * b.up[i][j ^ 1];
                c.up[i][j] += a.up[i][j] + b.up[i][j ^ (a.t >> i & 1)];
                c.down[i][j] += b.down[i][j] + a.down[i][j ^ (b.t >> i & 1)];
            }
        }
        return c;
    }
    struct Node {
        Node *ch[2], *p;
        Info val;
        Info tot;
        int cnt[D][2];
        i64 pair[D][2];
        i64 sum;
        Node() : ch{nullptr, nullptr}, p(nullptr), cnt {}, pair {}, sum {} {}
    };
    void pull(Node *t) {
        t->tot = (t->ch[0] ? t->ch[0]->tot : Info {}) + t->val + (t->ch[1] ? t->ch[1]->tot : Info {});
    }
    bool isroot(Node *t) {
        return t->p == nullptr || (t->p->ch[0] != t && t->p->ch[1] != t);
    }
    int pos(Node *t) {
        return t->p->ch[1] == t;
    }
    void rotate(Node *t) {
        Node *q = t->p;
        int x = !pos(t);
        q->ch[!x] = t->ch[x];
        if (t->ch[x]) {
            t->ch[x]->p = q;
        }
        t->p = q->p;
        if (!isroot(q)) {
            q->p->ch[pos(q)] = t;
        }
        t->ch[x] = q;
        q->p = t;
        pull(q);
    }
    void update(Node *t) {
        t->val.ans = t->val.t + t->sum;
        for (int i = 0; i < D; i++) {
            t->val.ans += (1LL << i) * t->pair[i][t->val.t >> i & 1];
            for (int j = 0; j < 2; j++) {
                t->val.up[i][j] = t->cnt[i][j ^ (t->val.t >> i & 1)];
                t->val.down[i][j] = t->cnt[i][j ^ (t->val.t >> i & 1)];
            }
            t->val.up[i][t->val.t >> i & 1]++;
            t->val.down[i][t->val.t >> i & 1]++;
        }
        pull(t);
    }
    void splay(Node *t) {
        while (!isroot(t)) {
            if (!isroot(t->p)) {
                if (pos(t) == pos(t->p)) {
                    rotate(t->p);
                } else {
                    rotate(t);
                }
            }
            rotate(t);
        }
        pull(t);
    }
    void add(Node *t, Info s) {
        for (int i = 0; i < D; i++) {
            for (int x = 0; x < 2; x++) {
                t->pair[i][x] += s.up[i][1 ^ x];
                for (int j = 0; j < 2; j++) {
                    t->pair[i][x] += t->cnt[i][j] * s.up[i][j ^ 1 ^ x];
                }
            }
            for (int j = 0; j < 2; j++) {
                t->cnt[i][j] += s.up[i][j];
            }
        }
        t->sum += s.ans;
    }
    void del(Node *t, Info s) {
        t->sum -= s.ans;
        for (int i = 0; i < D; i++) {
            for (int j = 0; j < 2; j++) {
                t->cnt[i][j] -= s.up[i][j];
            }
            for (int x = 0; x < 2; x++) {
                for (int j = 0; j < 2; j++) {
                    t->pair[i][x] -= t->cnt[i][j] * s.up[i][j ^ 1 ^ x];
                }
                t->pair[i][x] -= s.up[i][1 ^ x];
            }
        }
    }
    void access(Node *t, int v) {
        Info lst;
        for (Node *i = t, *q = nullptr; i; q = i, i = i->p) {
            splay(i);
            if (i->ch[1]) {
                add(i, i->ch[1]->tot);
            }
            i->ch[1] = q;
            if (q) {
                del(i, lst);
            } else {
                i->val.t = v;
            }
            lst = i->tot;
            update(i);
        }
        splay(t);
    }
}
```

# Gauss–Jordan Elimination Algorithm

```cpp
struct Matrix{
    const double eps = 1e-9;
    int n, m;
    double val[maxn][maxn];
    Matrix(int n, int m){
        this -> n = n;
        this -> m = m;
        memset(val, 0, sizeof(val));
    }
    double* operator [](int x){
        return val[x];
    }

    void swap(double &x, double &y){
        double temp = x;
        x = y;
        y = temp;
    }

    void swapLine(int x, int y){
        for(int i = 0; i < m; i++) swap(val[x][i], val[y][i]);
    }

    bool GE(void){
        if(n != m - 1) return false;
        for(int i = 0; i < n; i++){
            int maxArg = i;
            for(int j = i + 1; j < n; j++)
                if(std::abs(val[j][i]) > std::abs(val[maxArg][i])) maxArg = j;
            swapLine(maxArg, i);
            if(std::abs(val[i][i]) < eps) return false;
            for(int j = 0; j < n; j++){
                if(i == j) continue;
                double c = val[j][i] / val[i][i];
                for(int k = i; k < m; k++) val[j][k] -= val[i][k] * c;
            }
        }
        for(int i = 0; i < n; i++) val[i][n] /= val[i][i], val[i][i] = 1.0;
        return true;
    }
};
```

# Licao Segment Tree Data Structure

```cpp
template <typename Line>
struct LCTree{
    //一个类Line，和一个比较函数cmp，cmp(x, y, t)表示在t的时候，Line x是否优于Line y
    int root, n, tot;
    std::vector<int> ls, rs;
    std::vector<Line> line;
    std::function<bool(Line x, Line y, int t)> cmp = [](Line x, Line y, int t) -> bool{return true;};

    LCTree() = default;

    LCTree(int n, std::function<bool(Line x, Line y, int t)> cmp){
        this -> n = n;
        this -> tot = 0;
        this -> cmp = cmp;
        ls.resize(n << 1);
        rs.resize(n << 1);
        line.resize(n << 1);
        clear();
    }

    int newNode(Line v){
        line[++tot] = v;
        ls[tot] = rs[tot] = 0;
        return tot;
    }

    int insert(Line t, int now, int l, int r)
    {
        if(r < t.k) return 0; //注意，仅用在决策单调性中，如果不需要注意删除！！
        if(!now) return newNode(t);
        int mid = (l + r) >> 1;
        bool lc1 = cmp(t, line[now], l), rc1 = cmp(t, line[now], r);
        bool lc2 = cmp(line[now], t, l), rc2 = cmp(line[now], t, r);
        if(!lc1 && !rc1){return now;}
        if(!lc2 && !rc2){line[now] = t; return now;}
        if(lc1){
            if(cmp(t, line[now], mid)) rs[now] = insert(line[now], rs[now], mid + 1, r), line[now] = t;
            else ls[now]=insert(t, ls[now], l, mid);
        }
        else if(rc1){
            if(cmp(t, line[now], mid + 1)) ls[now] = insert(line[now], ls[now], l, mid), line[now] = t;
            else rs[now] = insert(t, rs[now], mid + 1, r);
        }
        return now;
    }

    void insert(Line t){
        root = insert(t, root, 1, n);
    }

    Line better(Line x, Line y, int t){
        return cmp(x, y, t) ? x : y;
    }

    Line ask(int x, int now, int l, int r){
        if(!now) return Line();
        int mid = (l + r) >> 1;
        if(x <= mid) return ls[now] ? better(ask(x, ls[now], l, mid), line[now], x) : line[now];
        else return rs[now] ? better(ask(x, rs[now], mid + 1, r), line[now], x) : line[now];
    }

    Line ask(int x){
        return ask(x, root, 1, n);
    }

    void clear(void){
        root = 0, tot = 0;
    }
};

int c[maxn], n, sum[maxn], mc, id[maxn];
unsigned long long dp[maxn];

unsigned long long pow(unsigned long long x){
    return x * x;
}

struct Line{
    int k, col;
    unsigned long long b;
    Line(void){
        k = col = b = 0;
    }
    Line(int k, unsigned long long b, int col){
        this -> k = k;
        this -> b = b;
        this -> col = col;
    }
    unsigned long long calc(int x){
        return b + pow(x + 1llu - k) * col;
    }
};

auto cmp = [](Line x, Line y, int t) -> bool{
    if(x.k > t || y.k > t) return x.calc(std::max(x.k, y.k)) > y.calc(std::max(x.k, y.k));
    return x.calc(t) > y.calc(t);
};
```

# Linear Basis Data Structure

```cpp
template<int S>
struct Hamel{
    std::bitset<S> val[S];
    ll cnt[S];
    Hamel(void) {}

    bool insert(std::bitset<S> x, ll c) {
        for(int i = S - 1; i >= 0; i--){
            if(x[i]){
                if(val[i] == 0){
                    val[i] = x;
                    cnt[i] = c;
                    return true;
                }
                else{
                    x = x ^ val[i];
                    if(x == 0) {
                        cnt[i] += c;
                        return false;
                    }
                }
            }
        }
        return false;
    }
    std::bitset<S> getMax(std::bitset<S> ans){
        for(int i = S - 1; i >= 0; i--){
            if(ans[i]) continue;
            ans = ans ^ val[i];
        }
        return ans;
    }
    std::bitset<S> getMin(std::bitset<S> ans){
        for(int i = S - 1; i >= 0; i--){
            if(ans[i]){
                ans = ans ^ val[i];
            }
        }
        return ans;
    }
    ll query(std::bitset<S> x) {
        ll ans = 1ll;
        for(int i = S - 1; i >= 0; i--){
            if(x[i]){
                if(val[i] == 0){
                    return 0ll;
                }
                else{
                    x = x ^ val[i];
                    ans = ans * (pow2[cnt[i] - 1]) % mod;
                }
            }
            else{
                if(cnt[i] != 0)
                    ans = ans * (pow2[cnt[i] - 1]) % mod;
            }
        }
        return ans;
    }
};
```

# Min Cost Max Flow Algorithm

```cpp
using ll = long long;
struct MCFGraph {
    struct Edge {
        int v, c, f;
        Edge(int v, int c, int f) : v(v), c(c), f(f) {}
    };
    const int n;
    std::vector<Edge> e;
    std::vector<std::vector<int>> g;
    std::vector<ll> h, dis;
    std::vector<int> pre;
    bool dijkstra(int s, int t) {
        dis.assign(n, std::numeric_limits<ll>::max());
        pre.assign(n, -1);
        std::priority_queue<std::pair<ll, int>, std::vector<std::pair<ll, int>>, std::greater<>> que;
        dis[s] = 0;
        que.emplace(0, s);
        while (!que.empty()) {
            ll d = que.top().first;
            int u = que.top().second;
            que.pop();
            if (dis[u] < d) continue;
            for (int i : g[u]) {
                int v = e[i].v;
                int c = e[i].c;
                int f = e[i].f;
                if (c > 0 && dis[v] > d + h[u] - h[v] + f) {
                    dis[v] = d + h[u] - h[v] + f;
                    pre[v] = i;
                    que.emplace(dis[v], v);
                }
            }
        }
        return dis[t] != std::numeric_limits<ll>::max();
    }
    MCFGraph(int n) : n(n), g(n) {}
    void addEdge(int u, int v, int c, int f) {
        g[u].push_back(e.size());
        e.emplace_back(v, c, f);
        g[v].push_back(e.size());
        e.emplace_back(u, 0, -f);
    }
    std::pair<int, ll> flow(int s, int t) {
        int flow = 0;
        ll cost = 0;
        h.assign(n, 0);
        while (dijkstra(s, t)) {
            for (int i = 0; i < n; ++i) h[i] += dis[i];
            int aug = std::numeric_limits<int>::max();
            for (int i = t; i != s; i = e[pre[i] ^ 1].v) aug = std::min(aug, e[pre[i]].c);
            for (int i = t; i != s; i = e[pre[i] ^ 1].v) {
                e[pre[i]].c -= aug;
                e[pre[i] ^ 1].c += aug;
            }
            flow += aug;
            cost += ll(aug) * h[t];
        }
        return std::make_pair(flow, cost);
    }
};

/** 费用流Double版 **/
struct MCFGraphDouble {
    const double eps = 1e-6;
    struct Edge {
        int v, c;
        double f;
        Edge(int v, int c, double f) : v(v), c(c), f(f) {}
    };
    const int n;
    std::vector<Edge> e;
    std::vector<std::vector<int>> g;
    std::vector<double> h, dis;
    std::vector<int> pre;
    bool dijkstra(int s, int t) {
        dis.assign(n, std::numeric_limits<double>::max());
        pre.assign(n, -1);
        std::priority_queue<std::pair<double, int>, std::vector<std::pair<double, int>>, std::greater<>> que;
        dis[s] = 0;
        que.emplace(0, s);
        while (!que.empty()) {
            double d = que.top().first;
            int u = que.top().second;
            que.pop();
            if (dis[u] < d) continue;
            for (int i : g[u]) {
                int v = e[i].v;
                int c = e[i].c;
                double f = e[i].f;
                if (c > 0 && dis[v] - eps > d + h[u] - h[v] + f) {
                    dis[v] = d + h[u] - h[v] + f;
                    pre[v] = i;
                    que.emplace(dis[v], v);
                }
            }
        }
        return dis[t] != std::numeric_limits<double>::max();
    }
    MCFGraph(int n) : n(n), g(n) {}
    void addEdge(int u, int v, int c, double f) {
        g[u].push_back(e.size());
        e.emplace_back(v, c, f);
        g[v].push_back(e.size());
        e.emplace_back(u, 0, -f);
    }
    std::pair<int, double> flow(int s, int t) {
        int flow = 0;
        double cost = 0;
        h.assign(n, 0);
        while (dijkstra(s, t)) {
            for (int i = 0; i < n; ++i) h[i] += dis[i];
            int aug = std::numeric_limits<int>::max();
            for (int i = t; i != s; i = e[pre[i] ^ 1].v) aug = std::min(aug, e[pre[i]].c);
            for (int i = t; i != s; i = e[pre[i] ^ 1].v) {
                e[pre[i]].c -= aug;
                e[pre[i] ^ 1].c += aug;
            }
            flow += aug;
            cost += double(aug) * h[t];
        }
        return std::make_pair(flow, cost);
    }
};
```

# Max Flow Algorithm

```cpp
constexpr int inf = 1E9;
template<class T>
struct MaxFlow {
    struct _Edge {
        int to;
        T cap;
        _Edge(int to, T cap) : to(to), cap(cap) {}
    };

    int n;
    std::vector<_Edge> e;
    std::vector<std::vector<int>> g;
    std::vector<int> cur, h;

    MaxFlow() {}
    MaxFlow(int n) {
        init(n);
    }

    void init(int n) {
        this->n = n;
        e.clear();
        g.assign(n, {});
        cur.resize(n);
        h.resize(n);
    }

    bool bfs(int s, int t) {
        h.assign(n, -1);
        std::queue<int> que;
        h[s] = 0;
        que.push(s);
        while (!que.empty()) {
            const int u = que.front();
            que.pop();
            for (int i : g[u]) {
                auto [v, c] = e[i];
                if (c > 0 && h[v] == -1) {
                    h[v] = h[u] + 1;
                    if (v == t) {
                        return true;
                    }
                    que.push(v);
                }
            }
        }
        return false;
    }

    T dfs(int u, int t, T f) {
        if (u == t) {
            return f;
        }
        auto r = f;
        for (int &i = cur[u]; i < int(g[u].size()); ++i) {
            const int j = g[u][i];
            auto [v, c] = e[j];
            if (c > 0 && h[v] == h[u] + 1) {
                auto a = dfs(v, t, std::min(r, c));
                e[j].cap -= a;
                e[j ^ 1].cap += a;
                r -= a;
                if (r == 0) {
                    return f;
                }
            }
        }
        return f - r;
    }
    void addEdge(int u, int v, T c) {
        g[u].push_back(e.size());
        e.emplace_back(v, c);
        g[v].push_back(e.size());
        e.emplace_back(u, 0);
    }
    T flow(int s, int t) {
        T ans = 0;
        while (bfs(s, t)) {
            cur.assign(n, 0);
            ans += dfs(s, t, std::numeric_limits<T>::max());
        }
        return ans;
    }

    std::vector<bool> minCut() {
        std::vector<bool> c(n);
        for (int i = 0; i < n; i++) {
            c[i] = (h[i] != -1);
        }
        return c;
    }

    struct Edge {
        int from;
        int to;
        T cap;
        T flow;
    };
    std::vector<Edge> edges() {
        std::vector<Edge> a;
        for (int i = 0; i < e.size(); i += 2) {
            Edge x;
            x.from = e[i + 1].to;
            x.to = e[i].to;
            x.cap = e[i].cap + e[i + 1].cap;
            x.flow = e[i + 1].cap;
            a.push_back(x);
        }
        return a;
    }
};
```

# Subset Dynamic Programming Algorithm

```cpp
const long long INF = 0x7ffffffffffffll;
using ll = long long;
const ll MOD = 1e9 + 7;

template<typename T>
void son_dp(std::vector<T> &dp, int limit) {
    limit = 1 << limit;
    for(int k = limit >> 1; k > 0; k >>= 1) {
        for(int i = 0; i < limit; i++) {
            if(i & k) {
                dp[i] += dp[i ^ k];
            }
        }
    }
} //输入初始数组dp和他的长度(2 ^ limit)，对它的每一个子集求和

using std::cin, std::cout, std::endl, std::vector;
int sz[1 << 21];
ll pow2[1 << 21];

//求解从n个物品中，选出全集的方案数，使用容斥+subset_dp
int main() {
    int n, N;
    cin >> n >> N;
    pow2[0] = 1ll;
    for(int i = 1; i < (1 << N); i++) {
        sz[i] = sz[i - (i & -i)] + 1;
    }
    std::vector<int> dp(1 << N, 0);
    for(int x, s, i = 1; i <= n; i++) {
        cin >> x;
        s = 0;
        for(int t, j = 0; j < x; j++) {
            cin >> t;
            s |= (1 << (t - 1));
        }
        dp[s]++;
        pow2[i] = pow2[i - 1] * 2 % MOD;
    }
    son_dp(dp, N);
    ll ans = 0;
    for(int i = 0; i < (1 << N); i++) {
        if((sz[i] & 1) ^ (sz[(1 << N) - 1] & 1)) {
            ans = (ans - pow2[dp[i]] + MOD) % MOD;
        }
        else {
            ans = (ans + pow2[dp[i]]) % MOD;
        }
    }
    cout << ans << endl;
}
```

# Fast Fourier Transform Algorithm

```cpp
using ll = long long;
using std::vector, std::swap, std::cin, std::cout, std::endl;

struct Complex {
    double real = 0, imag = 0;

    ~Complex() = default;

    Complex operator+ (const Complex &b) const {
        return {real + b.real, imag + b.imag};
    }

    Complex operator- (const Complex &b) const {
        return {real - b.real, imag - b.imag};
    }

    Complex operator* (const Complex &b) const {
        return {real * b.real - imag * b.imag, real * b.imag + imag * b.real};
    }
};

struct Poly : public std::vector<Complex> {

    using std::vector<Complex>::vector;

    bool operator<(const Poly &b) const {
        return size() > b.size();
    }
};

void FFT(Poly &A, vector<int> &R, int lim, int type) {
    for (int i = 0; i < lim; i++) {
        if (i < R[i]) {
            swap(A[i], A[R[i]]);
        }
    }
    for (int mid = 1; mid < lim; mid <<= 1) {
        Complex wn(cos(M_PI / mid), type * sin(M_PI / mid));
        for (int len = mid << 1, pos = 0; pos < lim; pos += len) {
            Complex w(1.0, 0.0);
            for (int k = 0; k < mid; ++k, w = w * wn) {
                Complex x = A[pos + k];
                Complex y = w * A[pos + mid + k];
                A[pos + k] = x + y;
                A[pos + mid + k] = x - y;
            }
        }
    }
    if (type == -1) {
        Complex inv(1.0 / lim, 0.0);
        for (int i = 0; i < lim; i++) {
            A[i] = A[i] * inv;
        }
    }
}

Poly operator*(Poly A, Poly B) {
    int lim = 1, lim_bit = 0;
    auto sz = A.size() + B.size() - 1;
    while (lim <= A.size() + B.size()) {
        lim <<= 1, lim_bit++;
    }
    A.resize(lim);
    B.resize(lim);
    vector<int> R(lim, 0);
    for (int i = 1; i < lim; i++) {
        R[i] = (R[i >> 1] >> 1) | ((i & 1) << (lim_bit - 1));
    }
    FFT(A, R, lim, 1);
    FFT(B, R, lim, 1);
    Poly res(lim);
    for (int i = 0; i < lim; i++) {
        res[i] = A[i] * B[i];
    }
    FFT(res, R, lim, -1);
    res.resize(sz);
    return res;
}
```

# Number Theoretic Transform Algorithm

```cpp
const int G = 3, Gi = 332748118, MOD = 998244353;

using ll = long long;
using std::vector, std::swap, std::cin, std::cout, std::endl;

ll power(ll a, ll b = MOD - 2, ll p = MOD) {
    ll res = 1;
    while (b) {
        if (b & 1) {
            res = res * a % p;
        }
        a = a * a % p;
        b >>= 1;
    }
    return res;
}

template<typename T>
void NTT(T &A, vector<int> &R, int lim, int type) {
    for (int i = 0; i < lim; i++) {
        if (i < R[i]) {
            swap(A[i], A[R[i]]);
        }
    }
    for (int mid = 1; mid < lim; mid <<= 1) {
        ll wn = power(type == 1 ? G : Gi, (MOD - 1) / (mid << 1));
        for (int len = mid << 1, pos = 0; pos < lim; pos += len) {
            long long w = 1;
            for (int k = 0; k < mid; ++k, w = w * wn % MOD) {
                long long x = A[pos + k];
                long long y = w * A[pos + mid + k] % MOD;
                A[pos + k] = (x + y) % MOD;
                A[pos + mid + k] = (x - y + MOD) % MOD;
            }
        }
    }
    if (type == -1) {
        long long inv = power(lim, MOD - 2);
        for (int i = 0; i < lim; i++) {
            A[i] = (A[i] * inv % MOD + MOD) % MOD;
        }
    }
}

struct Poly : public std::vector<long long> {

    using std::vector<long long>::vector;

    bool operator<(const Poly &b) const {
        return size() > b.size();
    }

    Poly deriv() const {
        if (this->empty()) {
            return Poly();
        }
        Poly res(this->size() - 1);
        for (int i = 0; i < this->size() - 1; ++i) {
            res[i] = (i + 1) * (*this)[i + 1] % MOD;
        }
        return res;
    }

    Poly integr() const {
        Poly res(this->size() + 1);
        for (int i = 0; i < this->size(); ++i) {
            res[i + 1] = (*this)[i] * power(i + 1) % MOD;
        }
        return res;
    }

    Poly shift(int k) const {
        if (k >= 0) {
            auto b = *this;
            b.insert(b.begin(), k, 0);
            return b;
        } else if (this->size() <= -k) {
            return {};
        } else {
            return Poly(this->begin() + (-k), this->end());
        }
    }

    Poly trunc(int k) const {
        Poly f = *this;
        f.resize(k);
        return f;
    }

    Poly operator+(const Poly &b) const {
        Poly ans(std::max(size(), b.size()), 0ll);
        for(int i = 0; i < size(); i++) {
            ans[i] += (*this)[i];
        }
        for(int i = 0; i < b.size(); i++) {
            ans[i] = (ans[i] + b[i]) % MOD;
        }
        return ans;
    }

    Poly operator-(const Poly &b) const {
        Poly ans(std::max(size(), b.size()), 0ll);
        for(int i = 0; i < size(); i++) {
            ans[i] += (*this)[i];
        }
        for(int i = 0; i < b.size(); i++) {
            ans[i] = (ans[i] - b[i]) % MOD;
        }
        return ans;
    }

    Poly operator*(Poly B) const {
        Poly A = *this;
        int lim = 1, lim_bit = 0;
        auto sz = A.size() + B.size() - 1;
        while (lim <= A.size() + B.size()) {
            lim <<= 1, lim_bit++;
        }
        A.resize(lim);
        B.resize(lim);
        vector<int> R(lim, 0);
        for (int i = 1; i < lim; i++) {
            R[i] = (R[i >> 1] >> 1) | ((i & 1) << (lim_bit - 1));
        }
        NTT(A, R, lim, 1);
        NTT(B, R, lim, 1);
        Poly res(lim, 0);
        for (int i = 0; i < lim; i++) {
            res[i] = A[i] * B[i] % MOD;
        }
        NTT(res, R, lim, -1);
        res.resize(sz);
        return res;
    }

    Poly operator* (ll v) const {
        auto ans = *this;
        for(auto &x : ans) {
            x = x * v % MOD;
        }
        return ans;
    }

    Poly inv(int m) const {
        Poly x{power((*this)[0], MOD - 2)}, poly_two{2ll};
        int k = 1;
        while(k < m) {
            k <<= 1;
            x = x * (poly_two - trunc(k) * x);
            x.resize(k);
        }
        x.resize(m);
        return x;
    }

    Poly log(int m) const {
        auto ans = (deriv() * inv(m)).integr();
        ans.resize(m);
        return ans;
    }

    Poly exp(int m) const {
        Poly x{1};
        int k = 1;
        while(k < m) {
            k <<= 1;
            x = x * (Poly{1} - x.log(k) + trunc(m));
            x.resize(k);
        }
        x.resize(m);
        return x;
    }

    Poly pow(int k, int m) const {
        int i = 0;
        while (i < this->size() && (*this)[i] == 0) {
            i++;
        }
        if (i == this->size() || 1LL * i * k >= m) {
            return Poly(m);
        }
        auto v = (*this)[i];
        auto f = shift(-i) * power(v);
        return (f.log(m - i * k) * k).exp(m - i * k).shift(i * k) * power(v, k);
    }
};

// 求原根
template<int P>
constexpr ll findPrimitiveRoot() {
    ll i = 2;
    int k = __builtin_ctz(P - 1);
    while (true) {
        if (power(i, (P - 1) / 2, P) != 1) {
            break;
        }
        i += 1;
    }
    return power(i, (P - 1) >> k, P);
}

int main() {

}
```

# Palindrome Automata Data Structure

```cpp
struct PAM {
    struct Node {
        int son[26], len, fail;
    };

    int node_tot, str_tot, last;
    Node node[maxn << 1];
    char str[maxn]{};

    PAM() {
        clear();
    }

    void clear() {
        node_tot = -1;
        str_tot = 0;
        last = 0;
        str[0] = '$';
        node[new_node()].len = 0;
        node[new_node()].len = -1;
        node[0].fail = 1;
    }

    int new_node() {
        node[++node_tot] = {};
        return node_tot;
    }

    int get_fail(int x, int index) {
        while (str[index - node[x].len - 1] != str[index]) {
            x = node[x].fail;
        }
        return x;
    }

    void insert(char ch) {
        str[++str_tot] = ch;
        int cur = get_fail(last, str_tot);
        int branch = ch - 'a';
        if(!node[cur].son[branch]) {
            int node_t = new_node();
            node[node_t].len = node[cur].len + 2;
            node[node_t].fail = node[get_fail(node[cur].fail, str_tot)].son[branch];
            node[cur].son[branch] = node_t;
        }
        last = node[cur].son[branch];
    }
} pam;
```

# Aho–Corasick Automaton Data Structure

```cpp
struct PAM {
    struct Node {
        int son[26], len, fail;
    };

    int node_tot, str_tot, last;
    Node node[maxn << 1];
    char str[maxn]{};

    PAM() {
        clear();
    }

    void clear() {
        node_tot = -1;
        str_tot = 0;
        last = 0;
        str[0] = '$';
        node[new_node()].len = 0;
        node[new_node()].len = -1;
        node[0].fail = 1;
    }

    int new_node() {
        node[++node_tot] = {};
        return node_tot;
    }

    int get_fail(int x, int index) {
        while (str[index - node[x].len - 1] != str[index]) {
            x = node[x].fail;
        }
        return x;
    }

    void insert(char ch) {
        str[++str_tot] = ch;
        int cur = get_fail(last, str_tot);
        int branch = ch - 'a';
        if(!node[cur].son[branch]) {
            int node_t = new_node();
            node[node_t].len = node[cur].len + 2;
            node[node_t].fail = node[get_fail(node[cur].fail, str_tot)].son[branch];
            node[cur].son[branch] = node_t;
        }
        last = node[cur].son[branch];
    }
} pam;
```

# Centroid Decomposition Algorithm

```cpp
const int N = 100005;
int n, m, ans[N];

std::vector<std::pair<int, int>> tr[N];
using std::cin, std::cout, std::endl;

namespace PointDivide {

    int vis[N], size[N], max_size[N], judge[10000005], dis[N], que[N];
    std::vector<int> dis_list;

    int get_root(int x, int p, int s) {
        int root = 0;
        std::function<void(int, int)> dfs = [&](int x, int p) {
            size[x] = 1;
            max_size[x] = 0;
            for (auto [y, v] : tr[x]) {
                if (vis[y] || y == p) {
                    continue;
                }
                dfs(y, x);
                size[x] += size[y];
                max_size[x] = std::max(max_size[x], size[y]);
            }
            max_size[x] = std::max(max_size[x], s - size[x]);
            if (!root || max_size[x] < max_size[root]) {
                root = x;
            }
        };
        dfs(x, p);
        return root;
    }

    void update(int x, int p) {
        if(dis[x] <= 1e7) {
            dis_list.push_back(dis[x]);
        }
        for (auto [y, v] : tr[x]) {
            if (vis[y] || y == p) {
                continue;
            }
            dis[y] = dis[x] + v;
            update(y, x);
        }
    }

    void solve(int x, int s) {
        vis[x] = 1;
        std::vector<int> temp_list;
        for (auto [y, v] : tr[x]) {
            if (vis[y]) continue;
            dis[y] = v;
            dis_list.clear();
            update(y, x);
            judge[0] = 1;
            for (auto dis_y : dis_list) {
                for (int k = 1; k <= m; k++) {
                    if (que[k] < dis_y) {
                        continue;
                    }
                    if (que[k] - dis_y == 0 || judge[que[k] - dis_y]) {
                        ans[k] = 1;
                    }
                }
            }
            temp_list.insert(temp_list.end(), dis_list.begin(), dis_list.end());
            for (auto dis_y : dis_list) {
                judge[dis_y] = 1;
            }
        }
        for (auto dis_t : temp_list) {
            judge[dis_t] = 0;
        }
        for (auto [y, v] : tr[x]) {
            if (vis[y]) {
                continue;
            }
            auto size_y = size[y] > size[x] ? s - size[x] : size[y];
            auto root = get_root(y, x, size_y);
            solve(root, size_y);
        }
    }
}

using namespace PointDivide;

//m个询问，每次询问树上有没有长度为k的路径

int main() {
    std::ios::sync_with_stdio(false);
    cin.tie(nullptr);

    cin >> n >> m;
    for (int x, y, z, i = 1; i < n; i++) {
        cin >> x >> y >> z;
        tr[x].emplace_back(y, z);
        tr[y].emplace_back(x, z);
    }
    for (int x, i = 1; i <= m; i++) {
        cin >> x;
        que[i] = x;
    }
    auto root = get_root(1, 0, n);
    solve(root, n);
    for (int i = 1; i <= m; i++) {
        cout << (ans[i] ? "AYE" : "NAY") << endl;
    }
}
```

# Suffix Automaton Data Structure

```cpp
const int maxn = 1e4 + 5;
struct Node {
	int link, len;
	long long size;
	int son[27];
	
	Node() {
		link = len = 0;
		memset(son, 0, sizeof(son));
	}
};

struct SAM {
	int root, tot, last;
	Node node[maxn << 1];
	
	SAM() {
		tot = 0;
		root = newNode();
		last = root;
	}
	
	int newNode() {
		node[++tot] = {};
		return tot;
	}
	
	void insert(char ch) {
		auto cur = newNode();
		int branch = ch - 'a';
		node[cur].size = 1;
		node[cur].len = node[last].len + 1;
		while(last != 0 && (node[last].son[branch] == 0)) {
			node[last].son[branch] = cur;
			last = node[last].link;
		}
		if(last == 0) {
			node[cur].link = root;
		}
		else {
			int q = node[last].son[branch];
			if(node[q].len == node[last].len + 1) {
				node[cur].link = q;
			}
			else {
				int clone = newNode();
				node[clone] = node[q];
				node[clone].len = node[last].len + 1;
				node[clone].size = 0;
				while(last != 0 && node[last].son[branch] == q) {
					node[last].son[branch] = clone;
					last = node[last].link;
				}
				node[q].link = clone;
				node[cur].link = clone;
			}
		}
		last = cur;
	}
	
	void getSize() {
		std::vector<int> cnt(tot + 1, 0);
		for(int i = root + 1; i <= tot; i++) {
			cnt[node[i].link]++;
		}
		
		std::queue<int> q;
		
		for(int i = root + 1; i <= tot; i++) {
			if(cnt[i] == 0) {
				q.push(i);
			}
		}
		while(!q.empty()) {
			int x = q.front();
			q.pop();
			int y = node[x].link;
			node[y].size += node[x].size;
			if(--cnt[y] == 0) {
				q.push(y);
			}
		}
	}
	
	void clear() {
		tot = 0;
		root = newNode();
		last = root;
	}
} sam;
```

# Border Series Algorithm

```cpp
using std::vector;
using std::string;

vector<int> getFail(const string &str) {
    vector<int> fail(str.size(), 0);
    fail[0] = -1;
    for(int i = 1; i < str.size(); i++) {
        int now = fail[i - 1];
        while(now >= -1 && str[now + 1] != str[i]) {
            now = now > -1 ? fail[now] : -2;
        }
        fail[i] = now + 1;
    }
    return fail;
}

int gcd(int a, int b, const vector<int> &fail) {
    a = fail[a], b = fail[b];
    while(a != b && (a != -1) && (b != -1)) {
        if(a < b) {
            std::swap(a, b);
        }
        if(fail[a] + 1 > (a + 1) / 2) {
            int d = a - fail[a];
            if((a + 1) % d == (b + 1) % d) {
                return b + 1;
            }
            else {
                a = (a + 1) % d + d - 1;
            }
        }
        else {
            a = fail[a];
        }
    }
    return std::min(a, b) + 1;
}

int main() {
    string s;
    std::cin >> s;
    auto fail = getFail(s);
    int n;
    std::cin >> n;
    while(n--) {
        int x, y;
        std::cin >> x >> y;
        std::cout << gcd(x - 1, y - 1, fail) << '\n';
    }
    return 0;
}
```

