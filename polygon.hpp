#ifndef POLYGON
#define POLYGON
// 平面凸包

#ifndef BASIS2D
#include "basis2D.hpp"
#endif

namespace Polygon{
    using namespace Basis2D;
    namespace Chull{
        Points __Graham_chull(const Points&_dots){
            Points dots{_dots};
            Point c = *std::min_element(dots.begin(), dots.end());
            std::sort(dots.begin(), dots.end(), [&c](auto a,auto b){
                auto ac = a - c, bc = b - c;
                double delta = ac.theta() - bc.theta();
                return std::abs(delta)>EPS?delta<0:ac.norm()<bc.norm();
            });
            Points H;
            for(int i=0;i<dots.size();i++){
                while(H.size()>1&&(H.back()-H[H.size()-2]).cross(dots[i]-H.back())<=0){
                    H.pop_back();
                }
                H.push_back(dots[i]);
            }
            return H;
        }
        Points __Andrew_chull(const Points&_dots){
            Points dots{_dots};
            std::sort(dots.begin(),dots.end());
            std::vector<int> I;
            std::vector<bool> used(dots.size(), true);
            for(int i=0;i<dots.size();i++){
                while(I.size()>1&&(dots[I.back()]-dots[I[I.size()-2]]).cross(dots[i]-dots[I.back()])<=0){
                    used[I.back()] = false; I.pop_back();
                }
                I.push_back(i);
            }
            used[0] = false;
            std::size_t s = I.size();
            for(int i=dots.size()-2;i>=0;i--){
                if(used[i]) continue;
                while(I.size()>s&&(dots[I.back()]-dots[I[I.size()-2]]).cross(dots[i]-dots[I.back()])<=0){
                    I.pop_back();
                }
                I.push_back(i);
            }
            Points H(I.size()-1);
            for(int i=0;i<I.size()-1;i++) H[i] = dots[I[i]];
            return H;
        }
        Points solve(const Points&dots){
            if(dots.size()<3)return dots;
            return __Andrew_chull(dots);
        }
    }
    namespace Get_diameter{
        typedef std::pair<int,int> resultType;
        std::pair<resultType,double> __get_diameter_rotating_calipers(const Points&dots){
            std::size_t n = dots.size();
            if(n<3)return {{0,n-1},(dots.back()-dots.front()).norm()};
            resultType res;
            double resv = 0.;
            int j = 2;
            for(int i=0;i<n;i++){
                while(std::abs((dots[(i+1)%n]-dots[i]).cross(dots[j]-dots[(i+1)%n])) <=
                    std::abs((dots[(i+1)%n]-dots[i]).cross(dots[(j+1)%n]-dots[(i+1)%n]))) j = (j+1)%n;
                auto tmp = std::max<std::pair<double,int>>({(dots[i]-dots[j]).norm(),i},
                    {(dots[(i+1)%n]-dots[j]).norm(),(i+1)%n});
                if(tmp.first>resv)resv=tmp.first,res={tmp.second,j};
            }
            return {res,resv};
        }
        std::pair<resultType,double> __solve(const Points&dots, bool isPolygon){
            auto get_diameter = __get_diameter_rotating_calipers;
            if(!isPolygon)return get_diameter(Chull::solve(dots));
            return get_diameter(dots);
        }
        resultType solve_pair(const Points&dots, bool isPolygon = false){
            return __solve(dots, isPolygon).first;
        }
        double solve_dist(const Points&dots, bool isPolygon = false){
            return __solve(dots, isPolygon).second;
        }
    }
}

#endif