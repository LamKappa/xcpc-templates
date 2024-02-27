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