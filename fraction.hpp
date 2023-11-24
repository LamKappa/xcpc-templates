#ifndef FRACTION
#define FRACTION

using i64 = long long;
using i128 = __int128;
using f80 = long double;
const i64 INF = -1ULL >> 15;

ostream& operator<<(ostream&out, i128 x){
    static int st[33];
    int top = 0;
    do{st[top++] = x%10;}while(x/=10);
    while(top>0) out<<st[--top];
    return out;
}

std::array<std::pair<i64,i64>, 2> fraction(f80 x, i64 M = INF, i64 N = INF){
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
struct Fraction{
    i64 a = 0, b = 1;
    i128 c = 0;

    Fraction(f80 x){
        auto[p1, p2] = fraction(std::abs(x));
        if(p1.first > p2.first) std::swap(p1, p2);
        if(p1.second == 0) std::swap(p1, p2);
        *this = Fraction((x<0?-1:1) * p1.first, p1.second);
    }
    Fraction(i128 a, i128 b, i128 c = 0){
        this->c = c + a / b;
        a %= b;
        constexpr i64 MAX_I64 = -1ULL>>1;
        auto mx = std::max(std::abs(a), std::abs(b));
        int l = 0, r = 64;
        if(mx <= MAX_I64) r = 0;
        while(l < r){
            int mid = (l+r) / 2;
            if((mx >> mid) <= MAX_I64){
                r = mid;
            }else l = mid + 1;
        }
        pack(a>>l, b>>l);
    }
    void pack(i64 a, i64 b){
        bool neg = (a<0) ^ (b<0);
        a = std::abs(a), b = std::abs(b);
        i64 g = std::__gcd(a, b);
        if(neg) a = -a;
        this->a = a / g;
        this->b = b / g;
    }

    operator f80()const{
        return c + (f80)a / b;
    }
    friend std::ostream& operator<<(std::ostream&out, const Fraction&f){
        return out<<"("<<f.c<<","<<f.a<<"/"<<f.b<<")";
    }
    Fraction operator+(const Fraction&o)const{
        i64 g = std::__gcd(b, o.b);
        return Fraction{a*(o.b/g) + o.a*(b/g), b*(o.b/g), c + o.c};
    }
    Fraction operator-(const Fraction&o)const{
        return *this + (-o);
    }
    Fraction operator-()const{
        return Fraction{-a, b, -c};
    }
    Fraction operator*(const Fraction&o)const{
        i64 g1 = std::__gcd(a, o.b);
        i64 g2 = std::__gcd(b, o.a);
        return Fraction{(a/g1) * (o.a/g2), (b/g2) * (o.b/g1), c * o.c} +
            Fraction{c*o.a, o.b} + Fraction{o.c*a, b};
    }
    Fraction operator/(const Fraction&o)const{
        return *this * o.inv();
    }
    Fraction inv()const{
        return Fraction{b, c*b + a};
    }

    Fraction operator+(const f80&o)const{
        return *this + Fraction(o);
    }
    Fraction operator-(const f80&o)const{
        return *this - Fraction(o);
    }
    Fraction operator*(const f80&o)const{
        return *this * Fraction(o);
    }
    Fraction operator/(const f80&o)const{
        return *this / Fraction(o);
    }
};

#endif