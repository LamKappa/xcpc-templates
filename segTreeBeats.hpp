#ifndef SEGTREE_BEATS
#define SEGTREE_BEATS
// 吉司机线段树 动态开点

const long long INF = 4e18;

struct Tag{
    long long addmax = 0, addsec = 0,
              hisaddmax = 0, hisaddsec = 0;

    void operator+=(const Tag&t){
        hisaddmax = std::max(hisaddmax, addmax + t.hisaddmax);
        hisaddsec = std::max(hisaddsec, addsec + t.hisaddsec);

        addmax += t.addmax;
        addsec += t.addsec;
    }
};

struct Info{
    long long sum = 0, max = -INF, sec = -INF, maxc = 0, secc = 0, maxhis = -INF;

    Info operator+(const Info&o){
        return {
            sum + o.sum,
            std::max(max, o.max),
            std::max(max == o.max ? -INF : std::min(max, o.max), std::max(sec, o.sec)),
            max == o.max ? maxc + o.maxc : (max > o.max ? maxc : o.maxc),
            (max == o.max ? 0 : (max > o.max ? o.maxc : maxc)) + secc + o.secc,
            std::max(maxhis, o.maxhis)
        };
    }
    void operator+=(const Tag&t){
        sum += t.addmax * maxc;
        sum += t.addsec * secc;

        maxhis = std::max(maxhis, max + t.hisaddmax);
        maxhis = std::max(maxhis, sec + t.hisaddsec);
        
        max += t.addmax;
        sec += t.addsec;
    }
};

template<class Info, class Tag>
struct SegTreeBeats{
    int left_margin, right_margin;
    struct Node{
        Info info;
        Tag tag;
        int ch[2] = {-1, -1};
    };
    std::vector<Node> node;

    SegTreeBeats(){}
    SegTreeBeats(int n):SegTreeBeats(){init(0,n);}
    SegTreeBeats(int left_margin, int right_margin):SegTreeBeats(){init(left_margin,right_margin);}
    SegTreeBeats(int n, Info v):SegTreeBeats(){init(std::vector(n,v));}
    SegTreeBeats(const std::vector<Info>&initarr):SegTreeBeats(){init(initarr);}
    void init(int left_margin, int right_margin){
        this->left_margin = left_margin;
        this->right_margin = right_margin;
        node.emplace_back();
    }
    void init(const std::vector<Info>&initarr){
        init(0,initarr.size());
        node.reserve(4<<std::__lg(right_margin - left_margin));
        auto build = [&](auto&&build,int rt,int l,int r)->void{
            if(l+1==r) return (void)(node[rt].info = initarr[l]);
            __push(rt);
            build(build, __reserve(node[rt].ch[0]), l, (l+r)/2);
            build(build, __reserve(node[rt].ch[1]), (l+r)/2, r);
            __pull(rt);
        };
        build(build, 0, left_margin, right_margin);
    }

    void __apply(int rt, const Tag&v){
        node[rt].info   += v;
        node[rt].tag    += v;
    }
    void __push(int rt){
        if(node[rt].ch[0]<0) node[rt].ch[0] = node.size(), node.emplace_back();
        if(node[rt].ch[1]<0) node[rt].ch[1] = node.size(), node.emplace_back();
        long long maxx = std::max(node[node[rt].ch[0]].info.max,
                node[node[rt].ch[1]].info.max);
        for(auto&ch : node[rt].ch){
            Tag tag = node[rt].tag;
            if(maxx > node[ch].info.max){
                tag.addmax = tag.addsec;
                tag.hisaddmax = tag.hisaddsec;
            }
            __apply(ch, tag);
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
    void add(int L, int R, long long v){
        if(L > R) return;
        __modify(0, left_margin, right_margin, L, R+1, [&](int rt,int l,int r){
            __apply(rt, Tag{
                v, v, v, v
            });
        });
    }
    void min(int L, int R, long long v){
        if(L > R) return;
        std::function<void(int,int,int)> func = [&](int rt,int l,int r){
            if(node[rt].info.sec >= v){
                __push(rt);
                func(node[rt].ch[0], l, (l+r)/2);
                func(node[rt].ch[1], (l+r)/2, r);
                __pull(rt);
            }else if(v < node[rt].info.max){
                int nv = v - node[rt].info.max;
                __apply(rt, Tag{
                    nv, 0, 0, 0
                });
            }
        };
        __modify(0, left_margin, right_margin, L, R+1, func);
    }
    void assign(int L, int R, long long v){
        if(L > R) return;
        add(L, R, INF);
        min(L, R, v);
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
        if(L > R) return Info();
        return __ask(0, left_margin, right_margin, L, R+1, [&](int rt,int l,int r){
            return node[rt].info;
        });
    }
};

#endif