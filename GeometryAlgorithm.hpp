#ifndef GEOMETRY_ALGORITHM_HPP
#define GEOMETRY_ALGORITHM_HPP

namespace GeometryAlgorithm{
    using namespace Geometry;
    using namespace ConvexHull;

    namespace CheckInside{
        template<typename ItA>
        bool inside_polygon(const Point&p, const ItA&pts){
            static std::mt19937 eng(std::random_device{}());
            static std::uniform_real_distribution<> dis(EPS,INF);

            auto begin = pts.begin(), end = pts.end();
            f80 K = dis(eng);
            f80 B = p.y - K*p.x;
            Line L{p, p + Point{1.l, K + B}};

            int cnt = 0;
            for(auto i=begin; i!=end; i++){
                auto j = std::next(i) == end? begin : std::next(i);
                Point pi{*i}, pj{*j};
                if(p == pi) return true;
                cnt += intersection_ray(L, Line{pi, pj}) && intersection_ray(L, Line{pj, pi});
            }
            return cnt % 2 == 1;
        }
    }

    namespace Diameter{
        std::tuple<int, int, f80> rotating_calipers(const Chull&H){
            int n = (int)H.size();
            if(n < 3) return {0, n - 1, (H.back() - H.front()).norm()};
            std::pair<int,int> diameter;
            f80 d = 0.l;
            for(int i=0, j=2; i<n; i++){
                while(fabsl(cross(H[(i + 1) % n] - H[i], H[j] - H[(i + 1) % n])) <=
                      fabsl(cross(H[(i + 1) % n] - H[i], H[(j + 1) % n] - H[(i + 1) % n]))) {
                    j = (j+1)%n;
                }
                auto tmp = std::max<std::pair<f80, int>>(
                        {(H[i] - H[j]).norm(), i},
                        {(H[(i + 1) % n] - H[j]).norm(), (i + 1) % n}
                );
                if(tmp.first > d) d = tmp.first, diameter = {tmp.second, j};
            }
            return {diameter.first, diameter.second, d};
        }
    }

    namespace ClosetPair{
        template<typename RA>
        std::tuple<int, int, f80> recursive(const RA&pts){
            std::pair<int, int> cls_pair;
            f80 d = INF;
            std::vector<int> sorted(pts.size());
            std::iota(sorted.begin(),sorted.end(),0);
            std::sort(sorted.begin(),sorted.end(),[&pts](auto a, auto b){
                return pts[a] < pts[b];
            });
            auto rec = [&, BF_n = std::max<int>(64, std::log2(pts.size()))](auto&&rec, int l,int r)->void{
                if(r - l <= BF_n){
                    std::sort(sorted.begin()+l, sorted.begin()+r, [&pts](auto a, auto b){
                        return pts[a].y < pts[b].y;
                    });
                    for(int i=l; i<r; i++){
                        while(pts[sorted[l]].y + d < pts[sorted[i]].y) l++;
                        for(int j=l; j<i; j++){
                            f80 tmp = (pts[sorted[i]] - pts[sorted[j]]).norm();
                            if(tmp < d) d = tmp, cls_pair = {sorted[i], sorted[j]};
                        }
                    }
                }else{
                    int mid = (l + r) / 2;
                    auto&middot = pts[sorted[mid]];
                    rec(rec, l, mid); rec(rec, mid, r);
                    std::inplace_merge(sorted.begin()+l, sorted.begin()+mid, sorted.begin()+r,
                                       [&pts](auto a, auto b){
                                           return pts[a].y < pts[b].y;
                                       });

                    std::array<std::deque<int>, 2> st;
                    for(int i=l; i<r; i++){
                        auto&dot = pts[sorted[i]];
                        if(fabsl(dot.x - middot.x) > d)continue;
                        bool lr = dot < middot;
                        while(!st[lr].empty() && pts[st[lr].front()].y + d < dot.y) st[lr].pop_front();
                        for(auto it=st[lr].rbegin(); it!=st[lr].rend(); it++){
                            f80 tmp = (dot - pts[*it]).norm();
                            if(tmp < d) d = tmp, cls_pair = {sorted[i], *it};
                        }
                        st[!lr].push_back(sorted[i]);
                    }
                }
            };
            rec(rec, 0, pts.size());
            return {cls_pair.first, cls_pair.second, d};
        }
        template<typename RA>
        std::tuple<int, int, f80> multiset(const RA&points){
            std::pair<int, int> cls_pair;
            f80 d = INF;
            std::vector<int> sortx(points.size());
            std::iota(sortx.begin(), sortx.end(), 0);
            std::sort(sortx.begin(), sortx.end(), [&points](auto a,auto b){
                return points[a] < points[b];
            });
            auto cmpy = [&points](auto a,auto b){
                return points[a].y == points[b].y ? points[a].x < points[b].x : points[a].y < points[b].y;
            };
            auto sorty = sortx;
            std::sort(sorty.begin(), sorty.end(), cmpy);
            std::multiset<int,decltype(cmpy)> s(cmpy);
            for(int i=0, l=0; i<points.size(); i++){
                auto&dot = points[sortx[i]];
                while(l<i && points[l].x + d < dot.x) s.erase(l++);
                auto it = s.lower_bound(
                        *lower_bound(sorty.begin(), sorty.end(), dot.y - d,
                                     [&points](auto a, f80 b){
                                         return points[a].y < b;
                                     })
                );
                while(it != s.end() && points[*it].y <= dot.y + d){
                    f80 tmp = (dot - points[*it]).norm();
                    if(tmp < d) d = tmp, cls_pair = {sortx[i], *it};
                    it++;
                }
                s.insert(sortx[i]);
            }
            return {cls_pair.first, cls_pair.second, d};
        }
    }

    namespace MinCoverage{
        std::pair<Points, f80> rectangle_coverage(const Chull&H){
            int n = (int)H.size();
            if(n <= 2) return {H, 0.l};
            int i = 0, j = 1;
            int k1 = 0, k2 = 0, k3 = 0;
            f80 ans = INF;
            Points ans_p(4);
            auto update = [&](){
                f80 H1 = dist(Line{H[i], H[j]}, H[k2]);
                f80 W1 = dist(Line{H[i], H[i] + (H[j] - H[i]).rotate(PI / 2)}, H[k1]);
                f80 W2 = dist(Line{H[i], H[i] + (H[j] - H[i]).rotate(PI / 2)}, H[k3]);
                f80 S = H1 * (W1 + W2);
                if(S < ans){
                    Point D = (H[j] - H[i]).unit();
                    ans_p[0] = H[i]; ans_p[0] = ans_p[0] + D * (-W1);
                    ans_p[1] = H[i]; ans_p[1] = ans_p[1] + D * (W2);
                    D = Point{-D.y, D.x};
                    ans_p[2] = ans_p[1]; ans_p[2] = ans_p[2] + D * H1;
                    ans_p[3] = ans_p[0]; ans_p[3] = ans_p[3] + D * H1;
                    ans = S;
                }
            };
            while(cross(H[j] - H[i], H[k1] - H[i]) <= cross(H[j] - H[i], H[(k1 + 1) % n] - H[i])){
                k1 = (k1+1)%n;
            }
            for(; i<n; i++, j=(j+1)%n){
                while(cross((H[j] - H[i]).rotate(PI / 2), H[k1] - H[i]) <= cross((H[j] - H[i]).rotate(PI / 2), H[(k1 + 1) % n] - H[i])){
                    k1 = (k1+1)%n;
                }
                while(cross(H[j] - H[i], H[k2] - H[i]) <= cross(H[j] - H[i], H[(k2 + 1) % n] - H[i])){
                    k2 = (k2+1)%n;
                }
                while(cross((H[j] - H[i]).rotate(PI / 2), H[k3] - H[i]) >= cross((H[j] - H[i]).rotate(PI / 2), H[(k3 + 1) % n] - H[i])){
                    k3 = (k3+1)%n;
                }
                update();
            }
            return {ans_p, ans};
        }
        template<typename ItA>
        std::pair<Circle, f80> circle_coverage(const ItA&_pts){
            Points points(_pts.begin(), _pts.end());
            std::shuffle(points.begin(), points.end(), std::random_device{});
            Circle C;
            for(auto i=points.begin(); i!=points.end(); i++){
                if(dist(C.c, *i) <= C.r) continue;
                C = {*i, 0.};
                for(auto j=points.begin(); j!=i; j++){
                    if(dist(C.c, *j) <= C.r) continue;
                    C.c = (*i + *j) * 0.5;
                    C.r = dist(C.c, *j);
                    for(auto k=points.begin(); k!=j; k++){
                        if(dist(C.c, *k) <= C.r) continue;
                        C = Circle{std::array<Point, 3>{*i, *j, *k}};
                    }
                }
            }
            return {C, C.area()};
        }
    }

    namespace HalfPlane{
        template<typename ItA>
        Points overlap(const ItA&_lines){
            Lines lines(_lines.begin(), _lines.end());
            std::sort(lines.begin(), lines.end(), [](auto a,auto b){
                return Point(a).theta() < Point(b).theta();
            });
            if(lines.empty() || Point(lines[0]).norm() < EPS) return {};

            std::deque<Point> chull;
            std::deque<Line> dq;
            for(auto L : lines){
                while(!chull.empty() && sign(cross(Point(L), chull.back() - L[0])) < 0){
                    chull.pop_back(); dq.pop_back();
                }
                while(!chull.empty() && sign(cross(Point(L), chull.front() - L[0])) < 0){
                    chull.pop_front(); dq.pop_front();
                }
                dq.push_back(L);
                if(dq.size() > 1){
                    if(sign(cross(Point(dq.back()), Point(dq[dq.size()-2]))) == 0){
                        Line bk = dq.back(); dq.pop_back();
                        if(sign(cross(Point(bk), Point(dq.back()[0] - bk[0]))) < 0){
                            if(sign(dot(Point(bk), Point(dq.back()))) < 0) return {};
                            dq.back() = L;
                        }
                        if(!chull.empty()) chull.pop_back();
                    }
                    if(dq.size() > 1) chull.push_back(intersection(dq.back(), dq[dq.size() - 2]));
                }
            }
            while(!chull.empty() && sign(cross(Point(dq.front()), chull.back() - dq.front()[0])) < 0){
                chull.pop_back(); dq.pop_back();
            }
            if(dq.size() > 1) chull.push_front(intersection(dq.back(), dq.front()));
            assert(chull.size() >= 3 && sign(cross(Point(dq.front()), Point(dq.back()))) <= 0);

            return Points{chull.begin(), chull.end()};
        }
    }

    namespace Area{
        template<typename ItA>
        f80 polygon_area(const ItA&pts){
            auto begin = pts.begin(), end = pts.end();
            if(begin == end) return 0.l;
            f80 res = 0;
            Point o = *begin;
            for(auto it=begin; std::next(it)!=end; it++){
                Point vec = (*std::next(it) - *it);
                res += cross(*it - o, vec) / 2.l;
            }
            return res;
        }
    }
}

#endif