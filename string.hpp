#ifndef STRING
#define STRING

namespace String{
    std::vector<int> manacher(const std::string&s){
        std::string t = "#";
        for(auto c : s){
            t += c; t += '#';
        }
        int n = t.size();
        std::vector<int> r(n);
        for(int i=0, j=0; i<n; i++){
            if(2*j-i >= 0 && j+r[j] > i) {
                r[i] = std::min(r[2*j-i], j+r[j]-i);
            }
            while(i-r[i] >= 0 && i+r[i] < n && t[i-r[i]] == t[i+r[i]]){
                r[i]++;
            }
            if(i+r[i] > j+r[j]) j = i;
        }
        return r;
    }

    std::vector<int> exkmp(const std::string&s) {
        int n = s.size();
        std::vector<int> z(n+1, 0); z[0] = n;
        for(int i=1, j=1; i<n; i++){
            z[i] = std::max(0, std::min(j+z[j]-i, z[i-j]));
            while(i+z[i] < n && s[z[i]] == s[i+z[i]]){
                z[i]++;
            }
            if(i+z[i] > j+z[j]) j = i;
        }
        return z;
    }
    std::vector<int> exkmp(const std::string&s, const std::string&t){
        auto z = exkmp(s + t);
        return {z.begin()+s.size(), z.end()};
    }

    template<int C=2,int N=int(2e4)>
    struct StringHash{
        typedef unsigned long long hash;
        typedef std::array<hash,C> Hash;
        static constexpr std::array<hash,14> Ms = {1000000007,1000000009,
            122420729,163227661,217636919,290182597,
            386910137,515880193,687840301,917120411,
            1222827239,1610612741,3221225473ul,4294967291ul};
        static constexpr hash B = 133;
        static constexpr std::array<Hash,N+1> norms = []()constexpr{
            assert(C<=Ms.size());
            std::array<Hash,N+1> norms;
            for(int mi=0;mi<C;mi++){
                auto M = Ms[mi];
                norms[0][mi] = 0ull;
                for(int i=1;i<=N;i++){
                    norms[i][mi] = norms[i-1][mi]*B%M;
                }
            }
            return norms;
        }();

        int n;
        std::vector<Hash> h;
        StringHash(){}
        StringHash(const std::string&s){
            init(s);
        }
        void init(const std::string&s){
            this->n = s.size();
            h.assign(n+1,{});
            for(int mi=0; mi<C; mi++){
                auto M = Ms[mi];
                for(int i=1; i<=n; i++){
                    h[i][mi] = (h[i-1][mi]*B%M + s[i-1])%M;
                }
            }
        }
        Hash ask(int l,int r){
            auto res = h[r+1];
            for(int mi=0; mi<C; mi++){
                auto M = Ms[mi];
                res[mi] = (res[mi] + M - h[l][mi]*norms[r-l+1][mi]%M)%M;
            }
            return res;
        }
        Hash ask(){return ask(0,n-1);}
    };

    struct SA{
        int n;
        std::vector<int> sa, rk, lc;
        SA(const std::string &s){
            n = s.length();
            sa.resize(n);
            lc.resize(n - 1);
            rk.resize(n);
            std::iota(sa.begin(), sa.end(), 0);
            std::sort(sa.begin(), sa.end(), [&](int a, int b){return s[a] < s[b];});
            rk[sa[0]] = 0;
            for(int i=1; i<n; i++){
                rk[sa[i]] = rk[sa[i-1]] + (s[sa[i]] != s[sa[i-1]]);
            }
            int k = 1;
            std::vector<int> tmp, cnt(n);
            tmp.reserve(n);
            while(rk[sa[n-1]] < n-1){
                tmp.clear();
                for(int i=0; i<k; i++) tmp.push_back(n - k + i);
                for(auto i : sa) if(i>=k) tmp.push_back(i - k);
                std::fill(cnt.begin(), cnt.end(), 0);
                for(int i=0; i<n; i++) ++cnt[rk[i]];
                for(int i=1; i<n; i++) cnt[i] += cnt[i-1];
                for(int i=n-1; i>=0; i--) sa[--cnt[rk[tmp[i]]]] = tmp[i];
                std::swap(rk, tmp);
                rk[sa[0]] = 0;
                for(int i=1; i<n; i++){
                    rk[sa[i]] = rk[sa[i-1]] + (tmp[sa[i-1]] < tmp[sa[i]] || sa[i-1] + k == n || tmp[sa[i-1] + k] < tmp[sa[i] + k]);
                }
                k *= 2;
            }
            for(int i=0, j=0; i < n; ++i){
                if(rk[i] == 0){
                    j = 0;
                }else{
                    for(j -= j > 0; i+j < n && sa[rk[i]-1]+j < n && s[i+j] == s[sa[rk[i]-1] + j]; ++j);
                    lc[rk[i]-1] = j;
                }
            }
        }
    };
}

#endif