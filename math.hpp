#ifndef MATH_CPP
#define MATH_CPP
// 数学

namespace Math{
    using i64 = long long;
    const i64 INF = -1ULL >> 3;
    const long double EPS = 1e-16;

    namespace Basic{
        i64 gcd(i64 a, i64 b){
            if((a&-a) < (b&-b)) std::swap(a, b);
            int bz = __builtin_ctzll(b);
            b >>= bz;
            while(a){
                a >>= __builtin_ctzll(a);
                i64 tmp = b - a;
                if(a < b) b = a;
                a = std::abs(tmp);
            }
            return b << bz;
        }
        std::array<std::pair<i64,i64>, 2> fraction(long double x, i64 M = INF, i64 N = INF){
            i64 m = 1, n = 1, lm=0, ln=1, rm=1, rn=0;
            while(m<=M && n<=N){
                int k = 0, fl = x*n > m;
                if(fl) lm = m, ln = n;
                else rm = m, rn = n;
                while(k>=0){
                    if(fl) m = lm + (rm<<k), n = ln + (rn<<k);
                    else m = (lm<<k) + rm, n = (ln<<k) + rn;
                    if(m>M || n>N) break;
                    if((x*n > m) == fl){
                        if(fl) lm = m, ln = n;
                        else rm = m, rn = n;
                        k++;
                    }else k--;
                }
                m=lm+rm; n=ln+rn;
            }
            return {std::make_pair(lm, ln), std::make_pair(rm, rn)};
        }
    }

    namespace Sieve{
        // for: linear_gcd
        // std::vector<std::array<int,3>> FAC_3;
        // std::vector<std::vector<int>> linear_gcd_lst;
        std::vector<int> minp, phi, primes, mu;


        void init(int N){
            // FAC_3.resize(N+1);
            minp.assign(N+1, 0);
            phi.assign(N+1, 0);
            mu.assign(N+1, 0); mu[1] = 1;
            primes.clear();
            
            for(int n=2;n<=N;n++){
                if(0==minp[n]){
                    // FAC_3[n]= {1, 1, n};
                    minp[n] = n;
                    phi[n]  = n-1;
                    mu[n]   = -1;
                    primes.push_back(n);
                }
                for(int p : primes){
                    if(n * p > N) break;
                    // FAC_3[n * p] = FAC_3[n]; FAC_3[n * p][0] *= p;
                    // std::sort(FAC_3[n * p].begin(), FAC_3[n * p].end());
                    minp[n * p] = p;
                    if(p==minp[n]){
                        phi[n * p]  = phi[n] * p;
                        mu[n * p]   = 0;
                        break;
                    }
                    phi[n * p]  = phi[n] * (p-1);
                    mu[n * p]   = mu[n] * mu[p];
                }
            }
            // int sqN = std::sqrt(N);
            // linear_gcd_lst.resize(sqN+1);
            // for(int i=1; i<=sqN; i++){
            //     linear_gcd_lst[i].resize(i+1);
            //     linear_gcd_lst[i][0] = i;
            //     for(int j=1; j<=i; j++){
            //         linear_gcd_lst[i][j] = linear_gcd_lst[j][i % j];
            //     }
            // }
        }

        // int linear_gcd(int a, int b){
        //     if(a > b) std::swap(a, b);
        //     if(a <= 1) return a ? 1 : b;
        //     int res = 1;
        //     for(auto k : FAC_3[a]){
        //         int tmp = 1;
        //         if(k > linear_gcd_lst.size()){
        //             if(b % k == 0) tmp = k;
        //         }else tmp = linear_gcd_lst[k][b % k];
        //         b /= tmp;
        //         res *= tmp;
        //     }
        //     return res;
        // }
    }

}

#endif