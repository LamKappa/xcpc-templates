#ifndef SEGTREE_H
#define SEGTREE_H
// 矩阵线段树 动态开点

#ifndef MATRIX
#include "matrix.hpp"
#endif

struct Tag{
    typedef Matrix<2,2> M;
    M mx;

    Tag(const M&_mx = {{
            {1, 0},
            {0, 1}
        }}) : mx(_mx){}
    
    void operator+=(const Tag&t){
        mx = mx * t.mx;
    }
};

struct Info{
    typedef Matrix<1,2> M;
    M mx;

    Info(const M&_mx = {{
            {0, 0}
        }}) : mx(_mx){}
    
    Info operator+(const Info&o){
        return mx + o.mx;
    }
    void operator+=(const Tag&t){
        mx = mx * t.mx;
    }
};

template<class Info, class Tag>
struct SegTree{
    int l_bound, r_bound;
    struct Node{
        Info info;
        Tag tag;
        int chl=-1, chr=-1;
    };
    std::vector<Node> node;

    SegTree(){}
    SegTree(int n):SegTree(){init(0,n);}
    SegTree(int l_bound, int r_bound):SegTree(){init(l_bound,r_bound);}
    SegTree(int n, Info v):SegTree(){init(std::vector(n,v));}
    SegTree(const std::vector<Info>&initarr):SegTree(){init(initarr);}
    void init(int l_bound, int r_bound){
        this->l_bound = l_bound;
        this->r_bound = r_bound;
        node.emplace_back();
    }
    void init(const std::vector<Info>&initarr){
        init(0,initarr.size());
        node.reserve(4<<std::__lg(r_bound - l_bound));
        std::function<void(int,int,int)> build = [&](int rt,int l,int r){
            if(l+1==r) return (void)(node[rt].info = initarr[l]);
            __push(rt, l, r);
            build(node[rt].chl, l, (l+r)/2); build(node[rt].chr, (l+r)/2, r);
            __pull(rt);
        };
        build(0, l_bound, r_bound);
    }

    void __apply(int rt, const Tag&v){
        node[rt].info += v;
        node[rt].tag += v;
    }
    void __push(int rt, int l, int r){
        if(node[rt].chl<0) node[rt].chl = node.size(), node.emplace_back();
        if(node[rt].chr<0) node[rt].chr = node.size(), node.emplace_back();
        __apply(node[rt].chl, node[rt].tag);
        __apply(node[rt].chr, node[rt].tag);
        node[rt].tag = Tag();
    }
    void __pull(int rt){
        node[rt].info = node[node[rt].chl].info + node[node[rt].chr].info;
    }
    template<typename Modifier>
    void __modify(int rt, int l, int r, int L, int R, const Modifier&func){
        if(R<=l || r<=L) return;
        if(L<=l && r<=R) return (void)(func(rt, l, r));
        __push(rt, l, r);
        if(L<(l+r)/2) __modify(node[rt].chl, l, (l+r)/2, L, R, func);
        if(R>(l+r)/2) __modify(node[rt].chr, (l+r)/2, r, L, R, func);
        __pull(rt);
    }
    void assign(int L, int R, const Info&v){
        __modify(0, l_bound, r_bound, L, R+1, [&](int rt,int l,int r){
            node[rt].info = v;
        });
    }
    void apply(int L, int R, const Tag&t){
        __modify(0, l_bound, r_bound, L, R+1, [&](int rt,int l,int r){
            __apply(rt, t);
        });
    }

    template<typename Filter>
    Info __ask(int rt, int l, int r, int L, int R, const Filter&func){
        if(R<=l || r<=L) return Info();
        if(L<=l && r<=R) return func(rt, l, r);
        __push(rt, l, r);
        return __ask(node[rt].chl, l, (l+r)/2, L, R, func) +
            __ask(node[rt].chr, (l+r)/2, r, L, R, func);
    }
    Info ask(int L, int R){
        return __ask(0, l_bound, r_bound, L, R+1, [&](int rt,int l,int r){
            return node[rt].info;
        });
    }

    static const int npos = -1;
    template<class Predicator>
    int __find(int rt, int l, int r, int L, int R, const Predicator&pred, bool forward){
        if(R<=l || r<=L || !pred(node[rt].info)) return npos;
        if(l+1==r) return l;
        __push(rt, l, r);
        int chl = node[rt].chl, chr = node[rt].chr;
        if(forward) std::swap(chl, chr);
        int pos = __find(chl, l, (l+r)/2, L, R, pred, forward);
        if(pos == npos){
            pos = __find(chr, (l+r)/2, r, L, R, pred, forward);
        }
        return pos;
    }
    template<class Predicator>
    int findFront(int L, int R, const Predicator&pred){
        return __find(0, l_bound, r_bound, L, R+1, pred, true);
    }
    template<class Predicator>
    int findBack(int L, int R, const Predicator&pred){
        return __find(0, l_bound, r_bound, L, R+1, pred, false);
    }
};

#endif