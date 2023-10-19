#ifndef ALGORITHM2D
#define ALGORITHM2D
// 平面点集算法

// depends basis2D.hpp

namespace Algorithm2D{
    using namespace Basis2D;
    namespace Check_inside{
        template<typename ForwardIt>
        bool __check_inside_ray(const Point&p,
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
                if((pi-p).cross(pi-pj)) return false;
            }
            return true;
        }
        template<typename Iterable>
        bool solve(const Point&p, const Iterable&dots, bool isChull = false){
            if(dots.size()<3)return false;
            auto check_inside = __check_inside_ray<typename Iterable::iterator>;
            if(isChull) check_inside = __check_inside_chull<typename Iterable::iterator>;
            return check_inside(p, dots.begin(), dots.end());
        }
    }
    namespace Closet_pair{
        typedef std::pair<int,int> resultType;
        template<typename RandomAccessable>
        std::pair<resultType,double> __closet_pair_rec(const RandomAccessable&dots){
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
        template<typename RandomAccessable>
        std::pair<resultType,double> __closet_pair_multiset(const RandomAccessable&dots){
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
        template<typename RandomAccessable>
        std::pair<resultType,double> __solve(const RandomAccessable&dots){
            return __closet_pair_rec(dots);
        }
        template<typename RandomAccessable>
        resultType solve_pair(const RandomAccessable&dots){
            return __solve(dots).first;
        }
        template<typename RandomAccessable>
        double solve_dist(const RandomAccessable&dots){
            return __solve(dots).second;
        }
    }
}

#endif