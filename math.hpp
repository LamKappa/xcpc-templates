#ifndef MATH
#define MATH
// 数学

namespace Math{
    using i64 = long long;
    using i128 = __int128;
    const i64 INF = -1ULL >> 3;
    const long double EPS = 1e-16;

    namespace Basic{
        i64 gcd(i64 a, i64 b){
            i64 neg = ((a < 0) ^ (b < 0)) ? -1 : 1;
            a = std::abs(a); b = std::abs(b);
            if((a&-a) < (b&-b)) std::swap(a, b);
            if(b == 0) return a;
            int bz = __builtin_ctzll(b);
            b >>= bz;
            while(a){
                a >>= __builtin_ctzll(a);
                i64 tmp = b - a;
                if(a < b) b = a;
                a = std::abs(tmp);
            }
            return neg * (b << bz);
        }
        i64 qmul(i64 a, i64 b, i64 mod){
            i64 res = a*b - mod*(i64)(1.L/mod*a*b);
            return res - mod*(res>=mod) + mod*(res<0);
        }
        i64 qpow(i64 a, i64 b, i64 mod){
            i64 res = 1ll;
            while(b){
                if(b & 1) res = qmul(res, a, mod);
                a = qmul(a, a, mod);
                b >>= 1;
            }
            return res;
        }
        bool miller_rabin(i64 x){
            constexpr std::array<i64, 7> test_p = {2, 325, 9375, 28178, 450775, 9780504, 1795265022};
            if(x<=3 || x%2 == 0) return x==2 || x==3;
            if(x%6 != 1 && x%6 != 5) return false;
            i64 k = x-1;
            int r = __builtin_ctzll(k), j;
            k >>= r;
            for(auto p : test_p){
                if(p % x == 0) continue;
                i64 v = qpow(p,k,x);
                if(v == 1) continue;
                for(j=1; j<=r; j++, v=qmul(v,v,x)) if(v==x-1) break;
                if(j > r) return false;
            }
            return true;
        }
        i64 pollard_rho(i64 num){
            if(num % 2 == 0) return 2;
            i64 val = 1ll, s = 0, t = 0;
            static std::mt19937 eng(std::random_device{}());
            i64 c = std::uniform_int_distribution<i64>(1ll, num-1)(eng);
            auto func = [](i64 x,i64 c,i64 mod){
                return (qmul(x,x,mod) + c) % mod;
            };
            for(int goal=1; ; goal<<=1, s=t, val=1ll){
                for(int steps=1; steps<=goal; ++steps){
                    t = func(t, c, num);
                    val = qmul(val, std::abs(t-s), num);
                    if(steps % 127 == 0){
                        i64 gcd_ = gcd(val, num);
                        if(gcd_ > 1) return gcd_;
                    }
                }
                i64 gcd_ = gcd(val, num);
                if(gcd_ > 1) return gcd_;
            }
            assert(false);
        }
        i64 max_prime_factor(i64 num){
            i64 max_factor = -1;
            auto dfs = [&max_factor](auto dfs, i64 num){
                if(num <= max_factor || num < 2) return;
                if(miller_rabin(num)){
                    return (void)(max_factor = std::max(max_factor, num));
                }
                i64 factor_ = pollard_rho(num);
                while(factor_ >= num) factor_ = pollard_rho(num);
                while(num % factor_ == 0) num /= factor_;
                dfs(dfs, num); dfs(dfs, factor_);
            };
            dfs(dfs, num);
            return max_factor;
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
        std::vector<int> minp, primes, phi, mu;

        void init(int N){
            minp.assign(N+1, 0);
            phi.assign(N+1, 0);
            mu.assign(N+1, 0); mu[1] = 1;
            primes.clear();
            
            for(int n=2;n<=N;n++){
                if(0==minp[n]){
                    minp[n] = n;
                    phi[n]  = n-1;
                    mu[n]   = -1;
                    primes.push_back(n);
                }
                for(int p : primes){
                    if(n * p > N) break;
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
        }

        std::vector<std::array<int,3>> FAC_3;
        std::vector<std::vector<int>> linear_gcd_lst;
        void linear_gcd_init(int N){
            FAC_3.resize(N+1);
            minp.assign(N+1, 0);
            primes.clear();
            
            for(int n=2;n<=N;n++){
                if(0==minp[n]){
                    FAC_3[n]= {1, 1, n};
                    minp[n] = n;
                    primes.push_back(n);
                }
                for(int p : primes){
                    if(n * p > N) break;
                    FAC_3[n * p] = FAC_3[n]; FAC_3[n * p][0] *= p;
                    std::sort(FAC_3[n * p].begin(), FAC_3[n * p].end());
                    minp[n * p] = p;
                    if(p==minp[n]) break;
                }
            }
            int sqN = std::sqrt(N);
            linear_gcd_lst.resize(sqN+1);
            for(int i=1; i<=sqN; i++){
                linear_gcd_lst[i].resize(i+1);
                linear_gcd_lst[i][0] = i;
                for(int j=1; j<=i; j++){
                    linear_gcd_lst[i][j] = linear_gcd_lst[j][i % j];
                }
            }
        }
        int linear_gcd(int a, int b){
            if(a > b) std::swap(a, b);
            if(a <= 1) return a ? 1 : b;
            int res = 1;
            for(auto k : FAC_3[a]){
                int tmp = 1;
                if(k > linear_gcd_lst.size()){
                    if(b % k == 0) tmp = k;
                }else tmp = linear_gcd_lst[k][b % k];
                b /= tmp;
                res *= tmp;
            }
            return res;
        }
    }

}

#endif