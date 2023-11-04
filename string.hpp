#ifndef STRING_HPP
#define STRING_HPP

namespace String{
    std::vector<int> manacher(const std::string&s){
        std::string t = "#";
        for(auto c : s){
            t += c;
            t += '#';
        }
        int n = t.size();
        std::vector<int> r(n);
        for(int i=0, j=0; i < n; i++){
            if(2*j - i >= 0 && j + r[j] > i) {
                r[i] = std::min(r[2 * j - i], j + r[j] - i);
            }
            while(i - r[i] >= 0 && i + r[i] < n && t[i - r[i]] == t[i + r[i]]){
                r[i] += 1;
            }
            if(i + r[i] > j + r[j]) j = i;
        }
        return r;
    }

    std::vector<int> exkmp(const std::string&s) {
        int n = s.size();
        std::vector<int> z(n + 1); z[0] = n;
        for(int i=1, j=1; i < n; i++){
            z[i] = std::max(0, std::min(j + z[j] - i, z[i - j]));
            while(i + z[i] < n && s[z[i]] == s[i + z[i]]){
                z[i]++;
            }
            if(i + z[i] > j + z[j])j = i;
        }
        return z;
    }

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
            for(int i=1; i < n; ++i){
                rk[sa[i]] = rk[sa[i - 1]] + (s[sa[i]] != s[sa[i - 1]]);
            }
            int k = 1;
            std::vector<int> tmp, cnt(n);
            tmp.reserve(n);
            while(rk[sa[n - 1]] < n - 1){
                tmp.clear();
                for(int i=0; i < k; ++i) tmp.push_back(n - k + i);
                for(auto i : sa) if (i >= k) tmp.push_back(i - k);
                std::fill(cnt.begin(), cnt.end(), 0);
                for(int i=0; i < n; ++i) ++cnt[rk[i]];
                for(int i=1; i < n; ++i) cnt[i] += cnt[i - 1];
                for(int i=n-1; i >= 0; --i) sa[--cnt[rk[tmp[i]]]] = tmp[i];
                std::swap(rk, tmp);
                rk[sa[0]] = 0;
                for(int i=1; i < n; ++i){
                    rk[sa[i]] = rk[sa[i - 1]] + (tmp[sa[i - 1]] < tmp[sa[i]] || sa[i - 1] + k == n || tmp[sa[i - 1] + k] < tmp[sa[i] + k]);
                }
                k *= 2;
            }
            for(int i=0, j=0; i < n; ++i){
                if(rk[i] == 0){
                    j = 0;
                }else{
                    for(j -= j > 0; i + j < n && sa[rk[i] - 1] + j < n && s[i + j] == s[sa[rk[i] - 1] + j]; ++j);
                    lc[rk[i] - 1] = j;
                }
            }
        }
    };
}

#endif