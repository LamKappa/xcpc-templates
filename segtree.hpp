#ifndef SEGTREE_H
#define SEGTREE_H
// 线段树

struct Tag {
    double mul = 1;
    double add = 0;
    
    void operator+=(const Tag&t){
        mul *= t.mul;
        add *= t.mul;
        add += t.add;
    }
};

struct Info {
    double sum = 0;
    int cnt = 0;
    
    Info operator+(const Info&o){
        return {sum+o.sum, cnt+o.cnt};
    }
    void operator+=(const Tag&t){
        sum *= t.mul;
        sum += t.add * cnt;
    }
};

template<class Info, class Tag>
struct SegTree{
    int n;
    std::vector<Info> node;
    std::vector<Tag> delay;

    SegTree(){}
    SegTree(int n):SegTree(){init(n);}
    SegTree(int n, Info v):SegTree(){init(std::vector(n,v));}
    SegTree(const std::vector<Info>&initarr):SegTree(){init(initarr);}
    void init(int n){
        this->n = n;
        node.assign(4<<std::__lg(n), Info());
        delay.assign(4<<std::__lg(n), Tag());
    }
    void init(const std::vector<Info>&initarr){
        init(initarr.size());
        std::function<void(int,int,int)> build = [&](int rt,int l,int r){
            if(l+1==r) return (void)(node[rt] = initarr[l]);
            int chl = 2*rt, chr = chl+1, mid = (l+r)/2;
            build(chl, l, mid); build(chr, mid, r);
            __pull(rt);
        };
        build(1, 0, n);
    }

    void __apply(int rt, const Tag&v){
        node[rt] += v;
        delay[rt] += v;
    }
    void __push(int rt){
        int chl = 2*rt, chr = chl+1;
        __apply(chl, delay[rt]);
        __apply(chr, delay[rt]);
        delay[rt] = Tag();
    }
    void __pull(int rt){
        int chl = 2*rt, chr = chl+1;
        node[rt] = node[chl] + node[chr];
    }
    template<typename Modifier>
    void __modify(int rt, int l, int r, int L, int R, Modifier func){
        if(R<=l || r<=L) return;
        if(L<=l && r<=R) return (void)(func(rt,l,r));
        int chl = 2*rt, chr = chl+1, mid = (l+r)/2;
        __push(rt);
        if(L<mid) __modify(chl, l, mid, L, R, func);
        if(R>mid) __modify(chr, mid, r, L, R, func);
        __pull(rt);
    }
    void mul(int L, int R, int v){
        __modify(1, 0, n, L, R+1, [&](int x,int l,int r){
            __apply(x,Tag{v,0});
        });
    }
    void add(int L, int R, int v){
        __modify(1, 0, n, L, R+1, [&](int x,int l,int r){
            __apply(x,Tag{1,v});
        });
    }

    Info __ask(int rt, int l, int r, int L, int R){
        if(R<=l || r<=L) return Info();
        if(L<=l && r<=R) return node[rt];
        int chl = 2*rt, chr = chl+1, mid = (l+r)/2;
        __push(rt);
        return __ask(chl, l, mid, L, R) + __ask(chr, mid, r, L, R);
    }
    Info ask(int L, int R){
        return __ask(1, 0, n, L, R+1);
    }

    static const int npos = -1;
    template<class Func>
    int __find(int rt, int l, int r, int L, int R, Func pred, bool forward){
        if(R<=l || r<=L || !pred(node[rt])) return npos;
        if(l+1==r) return l;
        int chl = 2*rt, chr = chl+1, mid = (l+r)/2;
        int pos = __find(chl^forward, l, mid, L, R, pred, forward);
        if(pos == npos){
            pos = __find(chr^forward, mid, r, L, R, pred, forward);
        }
        return pos;
    }
    template<class Func>
    int findFront(int L, int R, Func pred){
        return __find(1, 0, n, L, R+1, pred, false);
    }
    template<class Func>
    int findBack(int L, int R, Func pred){
        return __find(1, 0, n, L, R+1, pred, true);
    }
};

#endif