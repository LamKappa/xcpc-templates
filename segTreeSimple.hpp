#ifndef SEGTREE_SIMPLE
#define SEGTREE_SIMPLE
// 线段树

int p;

struct Tag {
    int mul = 1;
    int add = 0;
    
    void operator+=(const Tag&t){
        mul = (mul * t.mul)%p;
        add = (add * t.mul)%p;
        add = (add + t.add)%p;
    }
};

struct Info {
    int sum = 0;
    int cnt = 0;
    
    Info operator+(const Info&o){
        return {(sum+o.sum)%p, cnt+o.cnt};
    }
    void operator+=(const Tag&t){
        sum = (sum * t.mul)%p;
        sum = (sum + t.add*cnt)%p;
    }
};

struct SegTree{
    int n;
    std::vector<Info> node;
    std::vector<Tag> delay;
#define mid ((l+r)/2)
    SegTree(const std::vector<Info>&initarr){
        n = initarr.size();
        node.assign(4<<std::__lg(n), Info());
        delay.assign(4<<std::__lg(n), Tag());
        auto build = [&](auto&&build,int rt,int l,int r)->void{
            if(l+1==r) return (void)(node[rt] = initarr[l]);
            build(build, rt<<1, l, mid); build(build, rt<<1|1, mid, r);
            pull(rt);
        };
        build(build, 1, 0, n);
    }

    void push(int rt){
        node[rt<<1] += delay[rt];
        delay[rt<<1] += delay[rt];
        node[rt<<1|1] += delay[rt];
        delay[rt<<1|1] += delay[rt];
        delay[rt] = Tag();
    }
    void pull(int rt){
        node[rt] = node[rt<<1] + node[rt<<1|1];
    }
    void apply(int rt, int l, int r, int L, int R, const Tag&t){
        if(R<=l || r<=L) return;
        if(L<=l && r<=R){
            node[rt] += t;
            delay[rt] += t;
            return;
        }
        push(rt);
        apply(rt<<1, l, mid, L, R, t);
        apply(rt<<1|1, mid, r, L, R, t);
        pull(rt);
    }
    void apply(int L, int R, const Tag&t){
        if(L > R) return;
        apply(1, 0, n, L, R+1, t);
    }

    Info ask(int rt, int l, int r, int L, int R){
        if(R<=l || r<=L) return Info();
        if(L<=l && r<=R) return node[rt];
        push(rt);
        return ask(rt<<1, l, mid, L, R) +
            ask(rt<<1|1, mid, r, L, R);
    }
    Info ask(int L, int R){
        if(L > R) return Info();
        return ask(1, 0, n, L, R+1);
    }
#undef mid
};

#endif