#ifndef CONVEX_HULL_HPP
#define CONVEX_HULL_HPP

namespace ConvexHull{
    using namespace Geometry;

    template<typename ItA>
    std::array<Points, 2> Graham_chull(const ItA&_pts){
        Points pts(_pts.begin(), _pts.end());
        Point c = *std::min_element(pts.begin(), pts.end());
        std::sort(pts.begin(), pts.end(), [&c](auto&a, auto&b){return c.cmp_theta(a, b);});
        Points H, nH;
        for(auto&p : pts){
            while(H.size() > 1 && cross(H.back() - H[H.size() - 2], p - H.back()) <= 0){
                nH.push_back(H.back()); H.pop_back();
            }
            H.push_back(p);
        }
        return {H, nH};
    }
    template<typename ItA>
    std::array<Points, 2> Andrew_chull(const ItA&_pts){
        Points pts(_pts.begin(), _pts.end());
        std::sort(pts.begin(), pts.end());
        std::vector<int> I, nI;
        for(int i=0; i<pts.size(); i++){
            while(I.size() > 1 && Line{pts[I[I.size() - 2]], pts[I.back()]}.onLeft(pts[i]) <= 0){
                nI.push_back(I.back()); I.pop_back();
            }
            I.push_back(i);
        }
        std::size_t s = I.size();
        for(int i=(int)pts.size() - 2; i>=0; i--){
            while(I.size() > s && Line{pts[I[I.size() - 2]], pts[I.back()]}.onLeft(pts[i]) <= 0){
                nI.push_back(I.back()); I.pop_back();
            }
            I.push_back(i);
        }
        I.pop_back();
        Points H(I.size()), nH(nI.size());
        for(int i=0; i<I.size(); i++) H[i] = pts[I[i]];
        for(int i=0; i<nI.size(); i++) nH[i] = pts[nI[i]];
        return {H, nH};
    }
    template<typename ItA>
    Points minkowski_sum(const ItA&poly1, const ItA&poly2) {
        Points a{poly1.begin(), poly1.end()},
                b{poly2.begin(), poly2.end()};
        Point a0 = a[0], b0 = b[0];
        Points c{a0 + b0};
        for(int i = 0; i + 1 < a.size(); i++) a[i] = a[i + 1] - a[i];
        for(int i = 0; i + 1 < b.size(); i++) b[i] = b[i + 1] - b[i];
        a.back() = a0 - a.back(), b.back() = b0 - b.back();
        c.resize(a.size() + b.size() + 1);
        merge(a.begin(), a.end(), b.begin(), b.end(), c.begin() + 1,
              [](auto&L1, auto&L2){ return sign(cross(L1, L2)) > 0; });
        std::partial_sum(c.begin(), c.end(), c.begin());
        return c;
    }

    struct Chull : public Points{
        using Points::vector;
        Chull() = default;
        Chull(const Chull&H) = default;
        explicit Chull(const Points&pts) : Points(pts) {}

        Chull& init(){
            return *this = Chull{Andrew_chull(*this)[0]};
        }
        friend Chull operator+(const Chull&A, const Chull&B){
            return Chull{minkowski_sum(A, B)}.init();
        }
        bool inside(const Point&p) const{
            auto&H = *this;
            if(H[0].cmp_angular(p, H[1]) || H[0].cmp_angular(H.back(), p)) return false;
            auto itr = std::lower_bound(H.begin() + 2, std::prev(H.end()), 0, [&](auto&q, auto t){
                return H[0].cmp_angular(q, p);
            });
            auto p1 = *std::prev(itr), p2 = *itr;
            return Line{p1, p2}.onLeft(p) >= 0;
        }
        f80 area() const{
            return std::accumulate(begin(), end(), 0.l, [last_p=back()](auto s, const auto&p)mutable{
                s += cross(p - last_p, Point{} - last_p) / 2.l;
                last_p = p;
                return s;
            });
        }
        Point center() const{
            return std::accumulate(begin(), end(), Point{}, [last_p=back()](auto c, auto&p)mutable{
                c = c + (p + last_p) / 3.l * cross(p - last_p, Point{} - last_p) / 2.l;
                last_p = p;
                return c;
            }) / area();
        }
        std::array<std::size_t, 2> tangency_line(const Point&p) const{
            auto&H = *this;
            auto up = [&](int i, int j){
                if(i >= size()) i -= (int)size();
                if(j >= size()) j -= (int)size();
                return Line{H[i], H[j]}.onLeft(p);
            };
            auto binary_search = [&](const auto&&func){
                std::size_t l = 0, r = size() - 1;
                while(l < r){
                    std::tie(l, r) = func(l, r, (l + r) / 2);
                }
                return l;
            };
            auto p1 = binary_search([&](int l, int r, int m){
                bool f1 = up(l, l + 1) <= 0, f2 = up(m, m + 1) < 0, f3 = up(l, m) <= 0;
                if(!f1 && (up(r, r + 1) < 0 || up(l, r) > 0)) return std::make_pair(l, l);
                if((f1 && f2 && f3) || (!f1 && (f2 || !f3))) return std::make_pair(m + 1, r);
                return std::make_pair(l, m);
            });
            auto p2 = binary_search([&](int l, int r, int m){
                bool f1 = up(l, l + 1) >= 0, f2 = up(m, m + 1) >= 0, f3 = up(l, m) >= 0;
                if(!f1 && (up(r, r + 1) >= 0 || up(l, r) < 0)) return std::make_pair(l, l);
                if((f1 && f2 && f3) || (!f1 && (f2 || !f3))) return std::make_pair(m + 1, r);
                return std::make_pair(l, m);
            });
            return {p1, p2};
        }
    };
    std::array<Line, 4> tangency_line(const Chull&A, const Chull&B){
        if(A.size() > B.size()) return tangency_line(B, A);

        for(auto&p : A){
            auto[p1, p2] = B.tangency_line(p);
        }
        // todo  tangency_lines of two Convex Hull
        return {noLine, noLine, noLine, noLine};
    }
}

#endif