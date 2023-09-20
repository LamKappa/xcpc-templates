#ifndef POLYGON
#define POLYGON
// 平面凸包

#include "basis2D.hpp"

namespace Polygon{
    using namespace Basis2D;
    namespace Graham{
        Points chull(const Points&_dots){
            Points dots{_dots};
            Point c = *std::min_element(dots.begin(), dots.end());
            std::sort(dots.begin(), dots.end(), [&c](auto a,auto b){
                auto ac = a - c, bc = b - c;
                double delta = ac.theta() - bc.theta();
                return delta?delta<0:ac.norm()<bc.norm();
            });
            Points H{dots[0],dots[1]};
            for(int i=2;i<dots.size();i++){
                while(H.size()>1&&(H.back()-H[H.size()-2]).cross(dots[i]-H.back())<=0){
                    H.pop_back();
                }
                H.push_back(dots[i]);
            }
            return H;
        }
    }
    namespace Andrew{
        Points chull(const Points&_dots){
            Points dots{_dots};
            std::sort(dots.begin(),dots.end());
            std::vector<int> I{0};
            std::vector<bool> used(dots.size(), true);
            for(int i=1;i<dots.size();i++){
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
    }
}

#endif