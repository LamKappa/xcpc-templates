#ifndef MATH
#define MATH
// 数学

namespace Math {
    using i64 = long long;
    using i128 = __int128;
    const i64 INF = -1ULL >> 3;
    const long double EPS = 1e-16;

    namespace Sieve {
        std::vector<i64> minp, primes, phi, mu;
        
        void init(int N) {
            int last_n = std::max<int>(2, minp.size());
            if(primes.empty()){
                minp.assign(N+1, 0);
                phi.assign(N+1, 0);
                mu.assign(N+1, 0); mu[1] = 1;
                primes.clear();
            }
            
            for(int n=last_n; n<=N; n++) {
                if(0==minp[n]) {
                    minp[n] = n;
                    phi[n]  = n-1;
                    mu[n]   = -1;
                    primes.push_back(n);
                }
                for(int p : primes) {
                    if(n * p > N) { break; }
                    minp[n * p] = p;
                    if(p==minp[n]) {
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
        void linear_gcd_init(int N) {
            FAC_3.resize(N+1);
            minp.assign(N+1, 0);
            primes.clear();
            
            for(int n=2; n<=N; n++) {
                if(0==minp[n]) {
                    FAC_3[n]= {1, 1, n};
                    minp[n] = n;
                    primes.push_back(n);
                }
                for(int p : primes) {
                    if(n * p > N) { break; }
                    FAC_3[n * p] = FAC_3[n]; FAC_3[n * p][0] *= p;
                    std::sort(FAC_3[n * p].begin(), FAC_3[n * p].end());
                    minp[n * p] = p;
                    if(p==minp[n]) { break; }
                }
            }
            int sqN = std::sqrt(N);
            linear_gcd_lst.resize(sqN+1);
            for(int i=1; i<=sqN; i++) {
                linear_gcd_lst[i].resize(i+1);
                linear_gcd_lst[i][0] = i;
                for(int j=1; j<=i; j++) {
                    linear_gcd_lst[i][j] = linear_gcd_lst[j][i % j];
                }
            }
        }
        int linear_gcd(int a, int b) {
            if(a > b) { std::swap(a, b); }
            if(a <= 1) { return a ? 1 : b; }
            int res = 1;
            for(auto k : FAC_3[a]) {
                int tmp = 1;
                if(k > linear_gcd_lst.size()) {
                    if(b % k == 0) { tmp = k; }
                } else { tmp = linear_gcd_lst[k][b % k]; }
                b /= tmp;
                res *= tmp;
            }
            return res;
        }
    }
    
    namespace Basic {
        i64 gcd(i64 a, i64 b) {
            i64 neg = ((a < 0) ^ (b < 0)) ? -1 : 1;
            a = std::abs(a); b = std::abs(b);
            if((a&-a) < (b&-b)) { std::swap(a, b); }
            if(b == 0) { return a; }
            int bz = __builtin_ctzll(b);
            b >>= bz;
            while(a) {
                a >>= __builtin_ctzll(a);
                i64 tmp = b - a;
                if(a < b) { b = a; }
                a = std::abs(tmp);
            }
            return neg * (b << bz);
        }
        // a*x + b*y = gcd(a, b) returns min positive x
        std::array<i64, 3> exgcd(i64 a, i64 b){
            auto __exgcd = [](auto&&__exgcd, auto a, auto b, auto&x, auto&y)->i64{
                if(b==0){
                    x = 1; y = 0; return a;
                }
                auto d = __exgcd(__exgcd, b, a%b, x, y);
                std::tie(x, y) = std::make_pair(y, x - a / b * y);
                return d;
            };
            i64 x, y;
            auto d = __exgcd(__exgcd, a, b, x, y);
            auto tx = b / d;
            auto ty = a / d;
            x = x + tx * (i64)ceil((1. - x) / tx);
            y = (d - a * x) / b;
            return {x, y, d};
        }
        i64 qmul(i64 a, i64 b, i64 mod) {
            i64 res = a*b - mod*(i64)(1.L/mod*a*b);
            return res - mod*(res>=mod) + mod*(res<0);
        }
        i64 qpow(i64 a, int b){
            i64 res = 1;
            while(b) {
                if(b & 1) res *= a;
                a *= a;
                b >>= 1;
            }
            return res;
        }
        i64 qpow(i64 a, i64 b, i64 mod, i64 phi_mod = 0) {
            a %= mod;
            if(phi_mod > 0) b = std::min(b, b % phi_mod + phi_mod);
            i64 res = 1ll;
            while(b) {
                if(b & 1) { res = qmul(res, a, mod); }
                a = qmul(a, a, mod);
                b >>= 1;
            }
            return res;
        }
        bool miller_rabin(i64 x) {
            constexpr std::array<i64, 7> test_p = {2, 325, 9375, 28178, 450775, 9780504, 1795265022ll};
            if(x<=3 || x%2 == 0) { return x==2 || x==3; }
            if(x%6 != 1 && x%6 != 5) { return false; }
            i64 k = x-1;
            int r = __builtin_ctzll(k), j;
            k >>= r;
            for(auto p : test_p) {
                if(p % x == 0) { continue; }
                i64 v = qpow(p,k,x);
                if(v == 1) { continue; }
                for(j=1; j<=r; j++, v=qmul(v,v,x)) if(v==x-1) { break; }
                if(j > r) { return false; }
            }
            return true;
        }
        i64 pollard_rho(i64 num) {
            if(num % 2 == 0) { return 2; }
            i64 val = 1ll, s = 0, t = 0;
            static std::mt19937 eng(std::random_device{}());
            i64 c = std::uniform_int_distribution<i64>(1ll, num-1)(eng);
            auto func = [](i64 x,i64 c,i64 mod) {
                return (qmul(x,x,mod) + c) % mod;
            };
            for(int goal=1; ; goal<<=1, s=t, val=1ll) {
                for(int steps=1; steps<=goal; ++steps) {
                    t = func(t, c, num);
                    val = qmul(val, std::abs(t-s), num);
                    if(steps % 127 == 0) {
                        i64 gcd_ = gcd(val, num);
                        if(gcd_ > 1) { return gcd_; }
                    }
                }
                i64 gcd_ = gcd(val, num);
                if(gcd_ > 1) { return gcd_; }
            }
            assert(false);
        }
        i64 max_prime_factor(i64 num) {
            i64 max_factor = 1;
            auto dfs = [&max_factor](auto dfs, i64 num) {
                if(num <= max_factor || num < 2) { return; }
                if(miller_rabin(num)) {
                    return (void)(max_factor = std::max(max_factor, num));
                }
                i64 factor_ = pollard_rho(num);
                while(factor_ >= num) { factor_ = pollard_rho(num); }
                while(num % factor_ == 0) { num /= factor_; }
                dfs(dfs, num); dfs(dfs, factor_);
            };
            dfs(dfs, num);
            return max_factor;
        }
        std::vector<std::pair<i64, int>> factorize(i64 x){
            std::vector<std::pair<i64, int>> fac;
            while(x > 1){
                i64 p = (Sieve::minp.size() > x) ? Sieve::minp[x] : max_prime_factor(x), cnt = 0;
                for(auto pp=p, k=1ll; x%pp==0 || ((pp=p,k=1ll) && x%p==0); pp*=pp, k<<=1) x/=pp, cnt+=k;
                fac.emplace_back(p, cnt);
            }
            return fac;
        }
        i64 phi(i64 x, const std::vector<std::pair<i64, int>>&fac={}) {
            if(Sieve::phi.size() > x) return Sieve::phi[x];
            i64 res = x;
            for(auto [p, _] : fac.empty() ? factorize(x) : fac){
                res = res / p * (p-1);
            }
            return res;
        }
        i64 inv(i64 a, i64 mod) {
            auto[x, y, d] = exgcd(a, mod);
            if(d > 1) return -1;
            return x;
            // return qpow(a, phi(mod) - 1, mod);
        }
        std::array<i64, 2> primitive_root(i64 x) {
            // x is 2 || 4 || prime^k || 2*prime^k
            auto fac_x = factorize(x);
            if(fac_x.size() > 1){
                if(fac_x[0].first > fac_x[1].first) std::swap(fac_x[0], fac_x[1]);
                if(fac_x.size() != 2 || fac_x[0].second != 1){
                    return {-1};
                }
            }else if(x % 2 == 0){
                if(x <= 4) return {x-1, x/2};
                return {-1};
            }
            auto phi_x = phi(x, fac_x);
            auto fac_phi_x = factorize(phi_x);
            for(i64 g=2; g<x; g++){
                bool fl = true;
                if(gcd(g, x) > 1) continue;
                for(auto [fac, _] : fac_phi_x){
                    if(qpow(g, phi_x/fac, x) == 1){
                        fl = false;
                        break;
                    }
                }
                if(fl) return {g, phi_x};
            }
            return {-1};
        }
        // x == ai (mod mi)
        std::array<i64, 2> exCRT(const std::vector<std::pair<i64, i64>>&eq){
            if(eq.empty()) return {-1};
            auto[a0 ,m0] = eq.front(); a0 %= m0;
            for(int i=1; i<eq.size(); i++){
                auto[ai, mi] = eq[i]; ai %= mi;
                // (ai - a0) == x * m0 - y * mi
                auto[x, y, d] = exgcd(m0, mi);
                if((ai - a0) % d != 0) return {-1};
                mi = m0 / d * mi;
                x = qmul((ai - a0) / d, x, mi);
                a0 = (a0 + qmul(x, m0, mi)) % mi;
                m0 = mi;
            }
            return {a0, m0};
        }
        // ax == b (mod m)
        i64 liEu(i64 a, i64 b, i64 m){
            a %= m; b %= m;
            if(a == 0) return b == 0 ? 0 : -1;
            auto d = gcd(a, m);
            if(b % d != 0) return -1;
            // x = x0 + i * n / d
            return (b / d) * inv(a / d, m / d) % (m / d);
        }
        // c * a^x == b (mod m) exBSGS
        i64 exBSGS(i64 a, i64 b, i64 m, i64 c = 1) {
            a %= m; b %= m; c %= m;
            if(b == c || m == 1) return 0;
            // a^k/D * a^(x-k) == b/D (mod m/D)
            i64 k = 0, D = 1, akD = 1;
            for(i64 d; (d=gcd(a, m/D)) > 1; ){
                if(b/D % d != 0) return -1;
                k++; D *= d; akD = qmul(akD, (a / d), m);
                if(qmul(c, qmul(akD, D, m), m) == b){
                    return k;
                }
            }
            m = m / D; b = b / D;
            a %= m; b %= m; c = qmul(c, akD, m);
            
            i64 sqm = std::ceil(std::sqrt(m));

            std::unordered_map<i64, i64> baby_steps;
            for(i64 i=0, ax=1; i<sqm; i++){
                i64 baby_step = qmul(ax, b, m);
                baby_steps[baby_step] = i;
                ax = qmul(ax, a, m);
            }
            i64 x = -1;
            for(i64 j=1, asqm=qpow(a, sqm, m), ax=asqm; j<=sqm; j++){
                i64 giant_step = qmul(ax, c, m);
                if(baby_steps.count(giant_step)){
                    x = j * sqm - baby_steps[giant_step];
                    break;
                }
                ax = qmul(ax, asqm, m);
            }
            if(x == -1) return -1;
            return x + k;
        }
        // g^x == {b0, b1, b2, ...} (mod m)
        std::pair<i64, std::vector<i64>> logs(const std::vector<i64>&b_vec, i64 m, i64 g=-1, i64 phi_m=-1){
            // assert(gcd(bi, m) == 1)
            if(g == -1){
                std::tie(g, phi_m) = std::tuple_cat(primitive_root(m));
            }else if(phi_m == -1){
                phi_m = phi(m);
            }

            auto n = b_vec.size();
            i64 B = std::sqrt(phi_m * n);
            i64 BB = (phi_m + B - 1) / B;

            std::unordered_map<i64, i64> baby_steps;
            for(i64 i=0, gx=1; i<B; i++){
                baby_steps[gx] = i;
                gx = qmul(gx, g, m);
            }
            std::vector<i64> x(n, -1);
            for(int i=0; i<n; i++){
                auto b = b_vec[i];
                if(b == 1 || m == 1) { x[i] = 0; continue; }
                for(i64 j=1, gB=qpow(g, B, m), gx=qmul(gB, inv(b, m), m); j<=BB; j++){
                    if(baby_steps.count(gx)){
                        x[i] = j * B - baby_steps[gx];
                        break;
                    }
                    gx = qmul(gx, gB, m);
                }
            }
            return {g, x};
        }
        // x^a == b (mod m)
        // TODO: b==0 || m not prime
        std::vector<i64> LOG_BSGS(i64 a, i64 b, i64 m) {
            // g^(a*c) == b (mod m)
            // ac == t (mod phi(m))
            auto[g, phi_m] = primitive_root(m);
            if(g == -1) return {};

            auto t = exBSGS(g, b, m);
            if(t == -1) return {};
            auto c = liEu(a, t, phi_m);
            if(c == -1) return {};
            auto delta = phi_m / gcd(a, phi_m);
            std::vector<i64> ans;
            for(auto cur=c%delta; cur<phi_m; cur+=delta){
                ans.push_back(qpow(g, cur, m));
            }
            return ans;
        }
        std::array<std::pair<i64,i64>, 2> fraction(long double x, i64 M = INF, i64 N = INF) {
            i64 m = 1, n = 1, lm=0, ln=1, rm=1, rn=0;
            while(m<=M && n<=N) {
                int k = 0, fl = x*n > m;
                if(fl) { lm = m, ln = n; }
                else { rm = m, rn = n; }
                while(k>=0) {
                    if(fl) { m = lm + (rm<<k), n = ln + (rn<<k); }
                    else { m = (lm<<k) + rm, n = (ln<<k) + rn; }
                    if(m>M || n>N) { break; }
                    if((x*n > m) == fl) {
                        if(fl) { lm = m, ln = n; }
                        else { rm = m, rn = n; }
                        k++;
                    } else { k--; }
                }
                m=lm+rm; n=ln+rn;
            }
            return {std::make_pair(lm, ln), std::make_pair(rm, rn)};
        }
        
        struct Lucas{
            i64 mod;
            std::vector<i64> fac, inv, inv_fac;
            Lucas(i64 _mod) : mod(_mod){
                fac.resize(mod); 
                inv.resize(mod);
                inv_fac.resize(mod);
                fac[1] = inv[1] = inv_fac[1] = 1;
                for(int i=2; i<mod; i++){
                    fac[i] = qmul(fac[i-1], i, mod);
                    inv[i] = qmul(mod - mod / i, inv[mod % i], mod);
                    inv_fac[i] = qmul(inv_fac[i-1], inv[i], mod);
                }
            }
            i64 C(i64 n, i64 m)const{
                if(n <= m || m <= 0){
                    return n == m || m == 0;
                }
                if(n < mod){
                    return qmul(qmul(fac[n], inv_fac[m], mod), inv_fac[n - m], mod);
                }else{
                    return qmul(C(n % mod, m % mod), C(n / mod, m / mod), mod);
                }
            }
            i64 operator()(i64 n, i64 m)const{
                return C(n, m);
            }
        };

        struct exLucas{
            i64 mod;
            std::vector<std::pair<i64, int>> fac;
            std::vector<i64> pk;
            std::vector<std::vector<i64>> prod;
            exLucas(i64 _mod) : mod(_mod), fac(factorize(mod)){
                pk.resize(fac.size());
                prod.resize(fac.size());
                for(int i=0; i<fac.size(); i++){
                    auto[p, k] = fac[i];
                    pk[i] = qpow(p, k);
                    prod[i].resize(pk[i]); prod[i][0] = 1;
                    for(int j=0; j<pk[i]; j+=p){
                        if(j) prod[i][j] = prod[i][j-1];
                        for(int J=1; J<p; J++){
                            prod[i][j + J] = qmul(prod[i][j + J - 1], j + J, pk[i]);
                        }
                    }
                }
            }
            i64 operator()(i64 n, i64 m)const{
                if(n <= m || m <= 0){
                    return n == m || m == 0;
                }
                std::vector<std::pair<i64, i64>> eq(fac.size());
                for(int i=0; i<fac.size(); i++){
                    auto[p, k] = fac[i];
                    auto exponent = [&](i64 n)->i64{
                        i64 res = 0;
                        while(n /= p) res += n;
                        return res;
                    };
                    auto product = [&](i64 n)->i64{
                        i64 res = 1;
                        do{
                            if((n / pk[i]) % 2) res = pk[i] - res;
                            res = qmul(res, prod[i][n % pk[i]], pk[i]);
                        }while(n /= p);
                        return res;
                    };
                    auto e = exponent(n) - exponent(m) - exponent(n - m);
                    if(e >= k){
                        eq[i] = {0, pk[i]};
                    }else{
                        i64 c = qpow(p, e, pk[i], pk[i]-pk[i]/k);
                        c = qmul(c, product(n), pk[i]);
                        c = qmul(c, inv(product(m), pk[i]), pk[i]);
                        c = qmul(c, inv(product(n - m), pk[i]), pk[i]);
                        eq[i] = {c, pk[i]};
                    }
                }
                return exCRT(eq)[0];
            }
        };

        template<std::size_t B>
        struct BSGS{
            i64 m, g, cyc_g;
            std::unordered_map<i64, i64> baby_steps;
            BSGS(i64 _m, i64 _g=-1, i64 _cyc_g=-1) : m(_m), g(_g), cyc_g(_cyc_g){
                // g is generator of m with cyclic of cyc_g
                if(g == -1){
                    std::tie(g, cyc_g) = std::tuple_cat(primitive_root(m));
                }else if(cyc_g == -1){
                    cyc_g = phi(m);
                }

                for(i64 i=0, gx=1; i<B; i++){
                    baby_steps[gx] = i;
                    gx = qmul(gx, g, m);
                }
            }
            // g^x == b (mod m)
            i64 operator()(i64 b)const{
                if(b == 1 || m == 1) return 0;
                i64 BB = (cyc_g + B - 1) / B;

                for(i64 j=1, gB=qpow(g, B, m), gx=qmul(gB, inv(b, m), m); j<=BB; j++){
                    if(baby_steps.count(gx)){
                        return j * B - (baby_steps.find(gx)->second);
                    }
                    gx = qmul(gx, gB, m);
                }
                return -1;
            }
            // a^x == b (mod m)
            // xa*x == xb (mod cyc_g)
            i64 operator()(i64 a, i64 b)const{
                auto&bsgs = *this;
                return liEu(bsgs(a), bsgs(b), cyc_g);
            }
        };

        // based on Index Calculus
        template<std::size_t B>
        struct IndexCalculus{
            i64 mod, g, phi_m;
            std::vector<i64> primes, x_primes;
            IndexCalculus(i64 _mod, i64 _g=-1, i64 _phi_m=-1) : mod(_mod), g(_g), phi_m(_phi_m){
                if(g == -1){
                    std::tie(g, phi_m) = std::tuple_cat(primitive_root(mod));
                }else if(phi_m == -1){
                    phi_m = phi(mod);
                }

                Sieve::init(B);
                auto pn = std::upper_bound(Sieve::primes.begin(), Sieve::primes.end(), B) - Sieve::primes.begin() - 1;
                primes = {Sieve::primes.begin(), Sieve::primes.begin() + pn + 1};
                x_primes = logs(primes, mod, g, phi_m).second;
            }
            // g^x == b (mod m);
            i64 operator()(i64 b)const{
                static std::mt19937 eng{std::random_device{}()};
                for(i64 y = 0; ; y = eng() % phi_m){
                    auto z = qmul(b, qpow(g, y, mod, phi_m), mod);
                    auto max_p = max_prime_factor(z);
                    if(max_p > B) continue;
                    auto fac = factorize(z / max_p);
                    sort(fac.begin(), fac.end());
                    if(fac.empty() || fac.back().first != max_p){
                        if(max_p > 1) fac.emplace_back(max_p, 1);
                    }else fac.back().second++;
                    i64 ans = 0;
                    for(int i=0, j=0; i<fac.size(); i++){
                        auto[p, k] = fac[i];
                        while(primes[j] < p) j++;
                        ans = (ans + qmul(k, x_primes[j], phi_m)) % phi_m;
                    }
                    return (ans + phi_m - y) % phi_m;
                }
            }
        };

        // based on Pohlig–Hellman
        template<std::size_t B>
        struct PohligHellman{
            i64 m, g, phi_m;
            std::vector<std::pair<i64, int>> fac;
            std::vector<i64> pk, gi;
            std::vector<BSGS<B>> bsgs;
            PohligHellman(i64 _m, i64 _g=-1, i64 _phi_m=-1) : m(_m), g(_g), phi_m(_phi_m){
                if(g == -1){
                    std::tie(g, phi_m) = std::tuple_cat(primitive_root(m));
                }else if(phi_m == -1){
                    phi_m = phi(m);
                }

                fac = factorize(phi_m);
                pk.resize(fac.size());
                gi.resize(fac.size());
                bsgs.reserve(fac.size());
                for(int i=0; i<fac.size(); i++){
                    auto[p, k] = fac[i];
                    pk[i] = qpow(p, k);
                    gi[i] = qpow(g, phi_m / pk[i], m);
                    bsgs.emplace_back(m, gi[i], pk[i]);
                }
            }
            // g^x == h (mod m);
            i64 operator()(i64 h)const{
                std::vector<std::pair<i64, i64>> eq;
                for(int i=0; i<fac.size(); i++){
                    auto[p, k] = fac[i];
                    auto hi = qpow(h, phi_m / pk[i], m);
                    auto calc = [&](i64 gi, i64 hi)->i64{
                        i64 x = 0, pe = 1, ga = qpow(gi, pk[i]/p, m);
                        for(int e=0; e<k; e++, pe*=p){
                            auto he = qpow(qmul(hi, inv(qpow(gi, x, m), m), m), pk[i]/pe/p, m);
                            x = (x + qmul(pe, bsgs[i](ga, he), pk[i])) % pk[i];
                        }
                        return x;
                    };
                    auto xi = calc(gi[i], hi);
                    eq.emplace_back(xi, pk[i]);
                }
                return exCRT(eq)[0];
            }
            i64 operator()(i64 a, i64 b)const{
                auto&bsgs = *this;
                return liEu(bsgs(a), bsgs(b), phi_m);
            }
        };
    }
}

#endif