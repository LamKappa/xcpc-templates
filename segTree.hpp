#ifndef SEGTREE
#define SEGTREE
// 矩阵线段树 动态开点

// depends matrix.hpp

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
        int ch[2] = {-1, -1};
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
            __push(rt);
            build(node[rt].ch[0], l, (l+r)/2); build(node[rt].ch[1], (l+r)/2, r);
            __pull(rt);
        };
        build(0, l_bound, r_bound);
    }

    void __apply(int rt, const Tag&v){
        node[rt].info   += v;
        node[rt].tag    += v;
    }
    void __push(int rt){
        if(node[rt].ch[0]<0) node[rt].ch[0] = node.size(), node.emplace_back();
        if(node[rt].ch[1]<0) node[rt].ch[1] = node.size(), node.emplace_back();
        for(auto&ch : node[rt].ch){
            __apply(ch, node[rt].tag);
        }
        node[rt].tag = Tag();
    }
    void __pull(int rt){
        node[rt].info = node[node[rt].ch[0]].info + node[node[rt].ch[1]].info;
    }
    template<typename Modifier>
    void __modify(int rt, int l, int r, int L, int R, const Modifier&func){
        if(R<=l || r<=L) return;
        if(L<=l && r<=R) return (void)(func(rt, l, r));
        __push(rt);
        if(L<(l+r)/2) __modify(node[rt].ch[0], l, (l+r)/2, L, R, func);
        if(R>(l+r)/2) __modify(node[rt].ch[1], (l+r)/2, r, L, R, func);
        __pull(rt);
    }
    void assign(int L, int R, const Info&v){
        if(L > R) std::swap(L, R);
        __modify(0, l_bound, r_bound, L, R+1, [&](int rt,int l,int r){
            node[rt].info = v;
        });
    }
    void apply(int L, int R, const Tag&t){
        if(L > R) std::swap(L, R);
        __modify(0, l_bound, r_bound, L, R+1, [&](int rt,int l,int r){
            __apply(rt, t);
        });
    }

    template<typename Filter>
    Info __ask(int rt, int l, int r, int L, int R, const Filter&func){
        if(R<=l || r<=L) return Info();
        if(L<=l && r<=R) return func(rt, l, r);
        __push(rt);
        return __ask(node[rt].ch[0], l, (l+r)/2, L, R, func) +
            __ask(node[rt].ch[1], (l+r)/2, r, L, R, func);
    }
    Info ask(int L, int R){
        if(L > R) std::swap(L, R);
        return __ask(0, l_bound, r_bound, L, R+1, [&](int rt,int l,int r){
            return node[rt].info;
        });
    }

    static const int npos = -1;
    template<class Predicator>
    int __find(int rt, int l, int r, int L, int R, const Predicator&pred, bool forward){
        if(R<=l || r<=L || !pred(node[rt].info)) return npos;
        if(l+1==r) return l;
        __push(rt);
        int chl = node[rt].ch[0], chr = node[rt].ch[1], pos = npos;
        if(forward){
            pos = __find(chl, l, (l+r)/2, L, R, pred, forward);
            if(pos == npos){
                pos = __find(chr, (l+r)/2, r, L, R, pred, forward);
            }
        }else{
            pos = __find(chr, (l+r)/2, r, L, R, pred, forward);
            if(pos == npos){
                pos = __find(chl, l, (l+r)/2, L, R, pred, forward);
            }
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