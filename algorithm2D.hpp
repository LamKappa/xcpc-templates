#ifndef ALGORITHM2D
#define ALGORITHM2D
// 平面点集算法

#ifndef BASIS2D
#include "basis2D.hpp"
#endif

namespace Algorithm2D{
    using namespace Basis2D;
    namespace Check_inside{
        bool __check_inside_ray(Point&p, const Points&dots){
            static std::mt19937 eng(std::random_device{}());
            static std::uniform_real_distribution<> dis(0,PI/2.);
            double K = dis(eng);
            double B = p.y - K*p.x;

            int cnt = 0;
            for(int i=0;i<dots.size();i++){
                int j = i+1==dots.size()?0:i+1;
                Point pi{dots[i]}, pj{dots[j]}, pk;

                double k = (pi.y-pj.y) / (pi.x-pj.x);
                double b = pi.y - k*pi.x;

                pk.x = (b-B) / (K-k);
                if(!finite(k)) pk.x = pi.x;
                pk.y = K*pk.x + B;

                if(pk.x>=p.x && pk.between(pi,pj))cnt++;
            }
            return cnt&1;
        }
        bool solve(Point&p, const Points&dots){
            return __check_inside_ray(p, dots);
        }
    }
    namespace Closet_pair{
        typedef std::pair<int,int> resultType;
        std::pair<resultType,double> __closet_pair_rec(const Points&dots){
            resultType res;
            double resv = INF;
            std::vector<int> sorted(dots.size());
            std::iota(sorted.begin(),sorted.end(),0);
            sort(sorted.begin(),sorted.end(),[&dots](auto a,auto b){
                return dots[a] < dots[b];
            });
            int BF_n = std::max<int>(64, std::log2(dots.size()));
            std::function<void(int,int)> rec = [&](int l,int r){
                if(r-l<=BF_n){
                    std::sort(sorted.begin()+l,sorted.begin()+r,[&dots](auto a,auto b){
                        return dots[a].y < dots[b].y;
                    });
                    for(int i=l;i<r;i++){
                        while(dots[sorted[l]].y+resv<dots[sorted[i]].y) l++;
                        for(int j=l;j<i;j++){
                            double tmp = (dots[sorted[i]]-dots[sorted[j]]).norm();
                            if(tmp<resv)resv=tmp,res={sorted[i],sorted[j]};
                        }
                    }
                }else{
                    int mid = (l+r) / 2;
                    auto&middot = dots[sorted[mid]];
                    rec(l, mid); rec(mid, r);
                    std::inplace_merge(sorted.begin()+l,sorted.begin()+mid,sorted.begin()+r,
                        [&dots](auto a,auto b){
                            return dots[a].y < dots[b].y;
                    });

                    std::array<std::deque<int>, 2> st;
                    for(int i=l;i<r;i++){
                        auto&dot = dots[sorted[i]];
                        if(abs(dot.x-middot.x)>resv)continue;
                        bool lr = dot < middot;
                        while(!st[lr].empty()&&dots[st[lr].front()].y+resv<dot.y) st[lr].pop_front();
                        for(auto it=st[lr].rbegin();it!=st[lr].rend();it++){
                            double tmp = (dot-dots[*it]).norm();
                            if(tmp<resv)resv=tmp,res={sorted[i],*it};
                        }
                        st[!lr].push_back(sorted[i]);
                    }
                }
            };
            rec(0, dots.size());
            return {res,resv};
        }
        std::pair<resultType,double> __closet_pair_multiset(const Points&dots){
            resultType res;
            double resv = INF;
            std::vector<int> sortx(dots.size());
            std::iota(sortx.begin(),sortx.end(),0);
            std::sort(sortx.begin(),sortx.end(),[&dots](auto a,auto b){
                return dots[a] < dots[b];
            });
            auto cmpy = [&dots](auto a,auto b){
                return dots[a].y==dots[b].y?dots[a].x<dots[b].x:dots[a].y<dots[b].y;
            };
            auto sorty = sortx;
            std::sort(sorty.begin(),sorty.end(),cmpy);
            std::multiset<int,decltype(cmpy)> s(cmpy);
            for(int i=0,l=0;i<dots.size();i++){
                auto&dot = dots[sortx[i]];
                while(l<i&&dots[l].x+resv<dot.x) s.erase(l++);
                auto it = s.lower_bound(*lower_bound(sorty.begin(),sorty.end(),dot.y-resv,
                    [&dots](auto a,double b){
                        return dots[a].y<b;
                }));
                while(it!=s.end()&&dots[*it].y<=dot.y+resv){
                    double tmp = (dot-dots[*it]).norm();
                    if(tmp<resv)resv=tmp,res={sortx[i],*it};
                    it++;
                }
                s.insert(sortx[i]);
            }
            return {res,resv};
        }
        std::pair<resultType,double> __solve(const Points&dots){
            return __closet_pair_rec(dots);
        }
        resultType solve_pair(const Points&dots){
            return __solve(dots).first;
        }
        double solve_dist(const Points&dots){
            return __solve(dots).second;
        }
    }
}

#endif