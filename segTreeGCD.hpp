#ifndef SEGTREE_GCD
#define SEGTREE_GCD
// GCD线段树

template<typename T = long long>
struct SegTreeGCD{
    struct Info{
        T val = 0;
        Info operator+(const Info&o){
            return{
                std::abs(std::__gcd(val, o.val)),
            };
        }
    };
    int l_bound, r_bound;
    struct Node{
        Info info;
        int ch[2] ={-1, -1};
    };
    std::vector<Node> node;

    SegTreeGCD(){}
    SegTreeGCD(const std::vector<int>&initarr):SegTreeGCD(){init(initarr);}
    void init(int l_bound, int r_bound){
        this->l_bound = l_bound;
        this->r_bound = r_bound;
        node.emplace_back();
    }
    void init(const std::vector<int>&initarr){
        init(0,initarr.size());
        node.reserve(4<<std::__lg(r_bound - l_bound));
        std::function<void(int,int,int)> build = [&](int rt,int l,int r){
            if(l+1==r) return(void)(node[rt].info = Info{initarr[l] - (l?initarr[l-1]:0)});
            __push(rt);
            build(node[rt].ch[0], l,(l+r)/2); build(node[rt].ch[1],(l+r)/2, r);
            __pull(rt);
        };
        build(0, l_bound, r_bound);
    }

    void __push(int rt){
        for(auto&ch : node[rt].ch){
            if(ch<0) ch = node.size(), node.emplace_back();
        }
    }
    void __pull(int rt){
        node[rt].info = node[node[rt].ch[0]].info + node[node[rt].ch[1]].info;
    }
    template<typename Modifier>
    void __modify(int rt, int l, int r, int L, int R, const Modifier&func){
        if(R<=l || r<=L) return;
        if(L<=l && r<=R) return(void)(func(rt, l, r));
        __push(rt);
        if(L<(l+r)/2) __modify(node[rt].ch[0], l,(l+r)/2, L, R, func);
        if(R>(l+r)/2) __modify(node[rt].ch[1],(l+r)/2, r, L, R, func);
        __pull(rt);
    }
    void add(int X, T v){
        if(X < l_bound || r_bound <= X) return;
        __modify(0, l_bound, r_bound, X, X+1, [&](int rt,int l,int r){
            node[rt].info.val += v;
        });
    }
    void add(int L, int R, T v){
        if(L>R) return;
        add(L, v); add(R+1, -v);
    }
    template<typename Filter>
    Info __ask(int rt, int l, int r, int L, int R, const Filter&func){
        if(R<=l || r<=L) return Info();
        if(L<=l && r<=R) return func(rt, l, r);
        __push(rt);
        return __ask(node[rt].ch[0], l,(l+r)/2, L, R, func) +
            __ask(node[rt].ch[1],(l+r)/2, r, L, R, func);
    }
    T ask(int L, int R){
        if(L > R) std::swap(L, R);
        return __ask(0, l_bound, r_bound, L, R+1, [&](int rt,int l,int r){
            return node[rt].info;
        }).val;
    }
};

#endif