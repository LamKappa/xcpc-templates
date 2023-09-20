#include <bits/stdc++.h>
#pragma GCC optimize(3,"Ofast","inline")
#include <ext/pb_ds/assoc_container.hpp>
#include <ext/pb_ds/tree_policy.hpp>
using namespace std;
template<typename T>
using order_statistic_tree = __gnu_pbds::tree<
T,
__gnu_pbds::null_type,
std::less<T>,
__gnu_pbds::rb_tree_tag,
__gnu_pbds::tree_order_statistics_node_update>;
#define int long long
using i64 = long long;
using f128 = long double;

namespace tools{
    template<typename T=i64>
    T read(){
        char c;
        bool neg = false;
        while(!isdigit(c=getchar()))neg=c=='-';
        T ret = c^48;
        while(isdigit(c=getchar()))ret=ret*10+c-'0';
        return neg?-ret:ret;
    }
    template<typename T>
    void write(T x){
        static char st[80];
        int top = 0;
        if(x<(T)0){putchar('-');x=-x;}
        do{st[top++]=x%10;}while(x/=10);
        while(top)putchar(st[--top]^48);
    }
    void write(char c){putchar(c);}
    void write(const char*s){printf("%s",s);}
    void write(string s){write(s.c_str());}
    void write(f128 x){printf("%.0Lf",x);}
    void write(double x){printf("%lf",x);}
    void write(float x){printf("%f",x);}
    void write(bool x){write(x?"True":"False");}
    template<typename T, typename ... Args>
    void write(const T&x, Args ... args){
        write(x);putchar(' ');write(args...);
    }
    #define writeln(...) {write(__VA_ARGS__);putchar('\n');}
    struct pint{
        i64 x,y;
        pint(i64 x=0,i64 y=0):x(x),y(y){}
        bool operator<(const pint&o)const{
            return (x==o.x)?(y<o.y):(x<o.x);
        }
        struct hash{
            static size_t hash_val(const i64&val1,const i64&val2,size_t seed=0){
                seed ^= ::hash<i64>()(val1)+0x9e3779b9+(seed << 6)+(seed >> 2);
                seed ^= ::hash<i64>()(val2)+0x3a0f43eb+(seed << 6)+(seed >> 2);
                return seed;
            }
            static size_t splitmix64(size_t x){
                x += 0x9e3779b97f4a7c15;
                x = (x^(x>>30)) * 0xbf58476d1ce4e5b9;
                x = (x^(x>>27)) * 0x94d049bb133111eb;
                return x^(x>>31);
            }
            size_t operator()(const pint&p)const{
                return splitmix64(hash_val(p.x,p.y,p.x^p.y));
            }
        };
    };
    int lowbit(int x){return x&-x;}
    int log2(i64 n){
            int ret = 0;
            if(n&0xffffffff00000000){ret += 32; n >>= 32;}
            if(n&0x00000000ffff0000){ret += 16; n >>= 16;}
            if(n&0x000000000000ff00){ret += 8; n >>= 8;}
            if(n&0x00000000000000f0){ret += 4; n >>= 4;}
            if(n&0x000000000000000c){ret += 2; n >>= 2;}
            if(n&0x0000000000000002){ret += 1; n >>= 1;}
            return ret;
    }
    template<typename T=i64>
    T qpow(T x, T y,T MOD=1e9+7) {
        if(!(x%=MOD)) return 0;
        T ret = 1;
        while(y){
            if(y & 1) ret = ret*x%MOD;
            x = x*x%MOD;
            y >>= 1;
        }
        return ret;
    }
    template<typename T=i64>
    T gcd(T a, T b){
        if(a<b) swap(a,b);
        while(b){
            T t = b;
            b = a%b;
            a = t;
        }
        return a;
    }
    template<typename T=i64>
    void exgcd(T a,T b,T&x,T&y){
        if(!b){x=1;y=0;return;}
        exgcd(b,a%b,y,x);
        y -= a/b*x;
    }
    template<typename T=i64>
    T lcm(T x,T y){
        return x/gcd(x,y)*y;
    }
}
using pii = pair<int,int>;
using namespace tools;

// #define assert(con) {if(!(con)){*(int*)nullptr=0;}}

signed main(){
    // ios::sync_with_stdio(false);
    // cin.tie(nullptr);
    
    mt19937 eng(time(NULL));
    int T = 1;
    // T = read();
    while(T--){

    }
    return 0;
}