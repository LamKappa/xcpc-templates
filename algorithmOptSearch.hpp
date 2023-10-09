#ifndef ALGORITHM_OPT_SEARCH
#define ALGORITHM_OPT_SEARCH
// 优化查找通用算法

#ifndef RANDOM
#include "random.hpp"
#endif

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