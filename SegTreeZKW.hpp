#ifndef SEGTREE_ZKW_HPP
#define SEGTREE_ZKW_HPP

struct Info{
    int g = 0;

    Info operator+(const Info&o){
        return {__gcd(g, o.g)};
    }
};

template<class Info>
struct SegTreeZKW{
    int n, N;
    std::vector<Info> info;

    SegTreeZKW() = default;
    template<typename ...Args>
    explicit SegTreeZKW(Args&&... args){ init(std::forward<Args>(args)...); }

    void init(int n){
        this->n = n;
        this->N = 1 << (__lg(n) + 1);
        info.assign(N << 1, Info{});
    }

    void init(const std::vector<Info>&arr){
        init(arr.size());
        std::copy(arr.begin(), arr.end(), info.begin() + N);
        for(int x=N - 1; x; x--){
            info[x] = info[x * 2] + info[x * 2 + 1];
        }
    }

    void apply(int x, const Info&v){
        info[x += N] = v;
        while(x /= 2){
            info[x] = info[x * 2] + info[x * 2 + 1];
        }
    }

    Info ask(int l, int r){
        if(l > r) return {};
        l += N, r += N;
        Info infol = info[l], infor = l == r ? Info{} : info[r];
        for(; l + 1 < r; l /= 2, r /= 2){
            if(~l & 1) infol = infol + info[l ^ 1];
            if(r & 1) infor = info[r ^ 1] + infor;
        }
        return infol + infor;
    }
};

#endif