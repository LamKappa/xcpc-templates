#ifndef STRING_HASH
#define STRING_HASH

using u64 = unsigned long long;
struct SHash{
    static constexpr std::size_t C = 2;
    static constexpr u64 M[] = {
        1000000007,1118872217,122420729,163227661,
        217636919,290182597,386910137,515880193,687840301,
        917120411,1222827239,1610612741,3221225473ul,4294967291ul
    };

    std::array<u64, C> val;

    SHash() { val.fill(0ull); }

    SHash(const std::initializer_list<u64>&list) {
        std::copy(list.begin(), list.begin()+C, val.begin());
    }

    bool operator==(const SHash&o)const{
        return val == o.val;
    }
    bool operator<(const SHash&o)const{
        return val < o.val;
    }
    SHash operator*(const SHash&o)const{
        SHash res = *this;
        for(int i=0; i<C; i++){
            res.val[i] *= o.val[i];
            res.val[i] %= M[i];
        }
        return res;
    }
    SHash operator+(u64 y)const{
        SHash res = *this;
        for(auto&x : res.val) x += y;
        return res.norm();
    }
    SHash operator+(const SHash&o)const{
        SHash res;
        std::transform(val.begin(), val.end(), o.val.begin(), res.val.begin(), std::plus());
        return res.norm();
    }
    SHash operator-(u64 y)const{
        SHash res = *this;
        for(int i=0; i<C; i++){
            if(res.val[i] >= y) res.val[i] -= y;
            else res.val[i] += M[i] - y;
        }
        return res;
    }
    SHash operator-(const SHash&o)const{
        SHash res;
        for(int i=0; i<C; i++){
            if(val[i] >= o.val[i]) res.val[i] = val[i] - o.val[i];
            else res.val[i] = val[i] + (M[i] - o.val[i]);
        }
        return res.norm();
    }
    SHash norm(){
        for(int i=0; i<C; i++){
            if(val[i] >= M[i]) val[i] -= M[i];
        }
        return *this;
    }
};
namespace std {
    template <>
    struct hash<SHash> {
        std::size_t operator()(SHash x)const{
            return std::accumulate(x.val.begin(), x.val.end(), 0, std::bit_xor());
        }
    };
}
static const SHash B = {37, 53, 71, 97, 137, 251, 353, 491, 599, 617, 773, 853, 977, 1009};
constexpr int MAXN = 50005;
static const std::array<SHash, MAXN+1> Bp = [](){
    std::array<SHash, MAXN+1> Bp;
    Bp[0].val.fill(1ull);
    for(int i=1; i<=MAXN; i++){
        Bp[i] = Bp[i-1] * B;
    }
    return Bp;
}();

std::vector<SHash> build_prefix(std::string s){
    std::vector<SHash> prefix(s.size() + 1);
    for(int i=1; i<=s.size(); i++){
        prefix[i] = prefix[i-1] * B + s[i-1];
    }
    return prefix;
}

SHash get(const std::vector<SHash>&prefix, int l, int r){
    if(l > r) return SHash{};
    return prefix[r + 1] - prefix[l] * Bp[r - l + 1];
}

#endif