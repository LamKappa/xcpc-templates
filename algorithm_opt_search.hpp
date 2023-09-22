#ifndef ALGORITHM_OPT_SEARCH
#define ALGORITHM_OPT_SEARCH
// 优化查找通用算法

namespace Search{
    constexpr double EPS = 1e-8;
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
}

#endif