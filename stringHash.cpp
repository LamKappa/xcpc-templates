#ifndef STRINGHASH
#define STRINGHASH
// 字符串哈希

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
        for(int mi=0;mi<C;mi++){
            auto M = Ms[mi];
            for(int i=1;i<=n;i++){
                h[i][mi] = (h[i-1][mi]*B%M + s[i-1])%M;
            }
        }
    }
    Hash ask(int l,int r){
        auto res = h[r+1];
        for(int mi=0;mi<C;mi++){
            auto M = Ms[mi];
            res[mi] = (res[mi] + M - h[l][mi]*norms[r-l+1][mi]%M)%M;
        }
        return res;
    }
    Hash ask(){return ask(0,n-1);}
};

#endif