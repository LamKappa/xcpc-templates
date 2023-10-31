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
                        if(abs(dot.x-middot.x)>resv)continue;
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
    namespace HalfPlane{
        template<typename LineContainer>
        Points getIntersection(const LineContainer&_lines){
            Lines lines(_lines.begin(), _lines.end());
            std::sort(lines.begin(), lines.end(), [](auto a,auto b){
                return theta(lineVec(a)) < theta(lineVec(b));
            });
            if(lines.empty() || lineVec(lines[0]).norm() < EPS) return {};

            std::deque<Point> res;
            std::deque<Line> dq;
            for(auto l : lines){
                while(!res.empty() && cross(lineVec(l), res.back() - l[0]) < -EPS){
                    res.pop_back(); dq.pop_back();
                }
                while(!res.empty() && cross(lineVec(l), res.front() - l[0]) < -EPS){
                    res.pop_front(); dq.pop_front();
                }
                dq.push_back(l);
                if(dq.size()>1){
                    if(std::abs(cross(lineVec(dq.back()), lineVec(dq[dq.size()-2]))) < EPS){
                        Line bk = dq.back(); dq.pop_back();
                        if(cross(lineVec(bk), dq.back()[0] - bk[0]) < -EPS){
                            if(lineVec(bk) * lineVec(dq.back()) < -EPS) return {};
                            dq.back() = l;
                        }
                        if(!res.empty()) res.pop_back();
                    }
                    if(dq.size()>1) res.push_back(intersection(dq.back(), dq[dq.size()-2]));
                }
            }
            while(!res.empty() && cross(lineVec(dq.front()), res.back() - dq.front()[0]) < -EPS){
                res.pop_back(); dq.pop_back();
            }
            assert(cross(lineVec(dq.front()), lineVec(dq.back())) <= EPS);
            if(dq.size()>1) res.push_front(intersection(dq.back(), dq.front()));

            return Points(res.begin(), res.end());
        }
    }
}

#endif