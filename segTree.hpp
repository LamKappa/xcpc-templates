#ifndef SEGTREE_H
#define SEGTREE_H
// 线段树

struct Tag{
    int mul = 1;
    int add = 0;
    
    void operator+=(const Tag&t){
        mul *= t.mul;
        add *= t.mul;
        add += t.add;
    }
};

struct Info{
    int sum = 0;
    int cnt = 0;
    
    Info operator+(const Info&o){
        return {
            sum+o.sum,
            cnt+o.cnt
        };
    }
    void operator+=(const Tag&t){
        sum *= t.mul;
        sum += t.add * cnt;
    }
};

template<class Info, class Tag>
struct SegTree{
    int n, l_bound, r_bound;
    struct Node{
        Info info;
        Tag tag;
        int chl=-1, chr=-1;
    };
    std::vector<Node> node;

    SegTree(){}
    SegTree(int n):SegTree(){init(n,0,n);}
    SegTree(int l_bound, int r_bound):SegTree(){init(r_bound-l_bound,l_bound,r_bound);}
    SegTree(int n, Info v):SegTree(){init(std::vector(n,v));}
    SegTree(const std::vector<Info>&initarr):SegTree(){init(initarr);}
    void init(int n, int l_bound, int r_bound){
        this->n = n;
        this->l_bound = l_bound;
        this->r_bound = r_bound;
        node.emplace_back();
    }
    void init(const std::vector<Info>&initarr){
        init(initarr.size(),0,initarr.size());
        node.reserve(4<<std::__lg(n));
        std::function<void(int,int,int)> build = [&](int rt,int l,int r){
            if(l+1==r) return (void)(node[rt].info = initarr[l]);
            int mid = (l+r)/2;
            __push(rt);
            build(node[rt].chl, l, mid); build(node[rt].chr, mid, r);
            __pull(rt);
        };
        build(0, l_bound, r_bound);
    }

    void __apply(int rt, const Tag&v){
        node[rt].info += v;
        node[rt].tag += v;
    }
    void __push(int rt){
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
        if(L<=l && r<=R) return (void)(func(rt,l,r));
        int mid = (l+r)/2;
        __push(rt);
        if(L<mid) __modify(node[rt].chl, l, mid, L, R, func);
        if(R>mid) __modify(node[rt].chr, mid, r, L, R, func);
        __pull(rt);
    }
    void assign(int L, int R, const Info&v){
        __modify(0, l_bound, r_bound, L, R+1, [&](int x,int l,int r){
            node[x].info = v;
        });
    }
    void mul(int L, int R, int v){
        __modify(0, l_bound, r_bound, L, R+1, [&](int x,int l,int r){
            __apply(x,Tag{v,0});
        });
    }
    void add(int L, int R, int v){
        __modify(0, l_bound, r_bound, L, R+1, [&](int x,int l,int r){
            __apply(x,Tag{1,v});
        });
    }

    Info __ask(int rt, int l, int r, int L, int R){
        if(R<=l || r<=L) return Info();
        if(L<=l && r<=R) return node[rt].info;
        int mid = (l+r)/2;
        __push(rt);
        return __ask(node[rt].chl, l, mid, L, R) + __ask(node[rt].chr, mid, r, L, R);
    }
    Info ask(int L, int R){
        return __ask(0, l_bound, r_bound, L, R+1);
    }

    static const int npos = -1;
    template<class Predicator>
    int __find(int rt, int l, int r, int L, int R, const Predicator&pred, bool forward){
        if(R<=l || r<=L || !pred(node[rt])) return npos;
        if(l+1==r) return l;
        int mid = (l+r)/2;
        __push(rt);
        int pos = __find(forward?node[rt].chl:node[rt].chr, l, mid, L, R, pred, forward);
        if(pos == npos){
            pos = __find(forward?node[rt].chr:node[rt].chl, mid, r, L, R, pred, forward);
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