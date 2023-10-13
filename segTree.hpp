#ifndef SEGTREE_H
#define SEGTREE_H
// 线段树
// 动态开点
// 值域+权值
// 区间delayTag

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
    int l_bound, r_bound;
    struct Node{
        Info info;
        Tag tag;
        int l, r;
        int chl=-1, chr=-1;
        Node(int l, int r){
            this->l = l;
            this->r = r;
        }
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
        node.emplace_back(l_bound, r_bound);
    }
    void init(const std::vector<Info>&initarr){
        init(0,initarr.size());
        node.reserve(4<<std::__lg(r_bound - l_bound));
        std::function<void(int,int,int)> build = [&](int rt,int l,int r){
            if(l+1==r) return (void)(node[rt].info = initarr[l]);
            __push(rt);
            build(node[rt].chl, l, (l+r)/2); build(node[rt].chr, (l+r)/2, r);
            __pull(rt);
        };
        build(0, l_bound, r_bound);
    }

    void __apply(int rt, const Tag&v){
        node[rt].info += v;
        node[rt].tag += v;
    }
    void __push(int rt){
        if(node[rt].chl<0) node[rt].chl = node.size(), node.emplace_back(node[rt].l, (node[rt].l+node[rt].r)/2);
        if(node[rt].chr<0) node[rt].chr = node.size(), node.emplace_back((node[rt].l+node[rt].r)/2, node[rt].r);
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
        if(L<=l && r<=R) return (void)(func(rt));
        __push(rt);
        if(L<(l+r)/2) __modify(node[rt].chl, l, (l+r)/2, L, R, func);
        if(R>(l+r)/2) __modify(node[rt].chr, (l+r)/2, r, L, R, func);
        __pull(rt);
    }
    void assign(int L, int R, const Info&v){
        __modify(0, l_bound, r_bound, L, R+1, [&](int rt){
            node[rt].info = v;
        });
    }
    void mul(int L, int R, int v){
        __modify(0, l_bound, r_bound, L, R+1, [&](int rt){
            __apply(rt,Tag{v,0});
        });
    }
    void add(int L, int R, int v){
        __modify(0, l_bound, r_bound, L, R+1, [&](int rt){
            __apply(rt,Tag{1,v});
        });
    }

    Info __ask(int rt, int l, int r, int L, int R){
        if(R<=l || r<=L) return Info();
        if(L<=l && r<=R) return node[rt].info;
        __push(rt);
        return __ask(node[rt].chl, l, (l+r)/2, L, R) + __ask(node[rt].chr, (l+r)/2, r, L, R);
    }
    Info ask(int L, int R){
        return __ask(0, l_bound, r_bound, L, R+1);
    }

    static const int npos = -1;
    template<class Predicator>
    int __find(int rt, int l, int r, int L, int R, const Predicator&pred, bool forward){
        if(R<=l || r<=L || !pred(node[rt].info)) return npos;
        if(l+1==r) return l;
        __push(rt);
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