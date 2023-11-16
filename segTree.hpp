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
    struct Node{
        Info info;
        Tag tag;
        int ch[2] = {-1, -1};
    };
    int root = -1;
    static int left_margin, right_margin;
    static std::vector<Node> node;
    static int drop;

    SegTree(){init(left_margin,right_margin);}
    SegTree(int n):SegTree(){init(0,n);}
    SegTree(int left_margin, int right_margin):SegTree(){init(left_margin,right_margin);}
    SegTree(int n, Info v):SegTree(){init(std::vector(n,v));}
    SegTree(const std::vector<Info>&initarr):SegTree(){init(initarr);}
    static void resize(int left_margin, int right_margin){
        SegTree::left_margin = std::min(SegTree::left_margin, left_margin);
        SegTree::right_margin = std::max(SegTree::right_margin, right_margin);
    }
    void init(int left_margin, int right_margin){
        resize(left_margin, right_margin);
        __reserve(root);
    }
    void init(const std::vector<Info>&initarr){
        init(0,initarr.size());
        node.reserve(4<<std::__lg(right_margin - left_margin));
        std::function<void(int,int,int)> build = [&](int rt,int l,int r){
            if(l+1==r) return (void)(node[rt].info = initarr[l]);
            __push(rt);
            build(__reserve(node[rt].ch[0]), l, (l+r)/2);
            build(__reserve(node[rt].ch[1]), (l+r)/2, r);
            __pull(rt);
        };
        build(root, left_margin, right_margin);
    }

    static void __drop(int x){
        node[x].ch[0] = drop;
        drop = x;
    }
    static int __reserve(int&x){
        if(x != -1) return x;
        int res = node.size();
        if(drop == -1){
            x = res;
            node.emplace_back();
        }else{
            x = res = drop;
            drop = node[x].ch[0];
            node[x] = Node();
        }
        return res;
    }
    static void __apply(int rt, const Tag&v){
        if(rt == -1) return;
        node[rt].info   += v;
        node[rt].tag    += v;
    }
    static void __push(int rt){
        for(auto&ch : node[rt].ch){
            __apply(ch, node[rt].tag);
        }
        node[rt].tag = Tag();
    }
    static void __pull(int rt){
        node[rt].info = (node[rt].ch[0]==-1 ? Info() : node[node[rt].ch[0]].info) +
            (node[rt].ch[1]==-1 ? Info() : node[node[rt].ch[1]].info);
    }

    template<typename Modifier>
    static void __modify(int rt, int l, int r, int L, int R, const Modifier&func){
        if(R<=l || r<=L) return;
        if(L<=l && r<=R) return (void)(func(rt, l, r));
        if(L<(l+r)/2) __reserve(node[rt].ch[0]);
        if(R>(l+r)/2) __reserve(node[rt].ch[1]);
        __push(rt);
        __modify(node[rt].ch[0], l, (l+r)/2, L, R, func);
        __modify(node[rt].ch[1], (l+r)/2, r, L, R, func);
        __pull(rt);
    }
    void assign(int L, int R, const Info&v){
        if(L > R) return;
        __modify(root, left_margin, right_margin, L, R+1, [&](int rt,int l,int r){
            node[rt].info = v;
        });
    }
    void append(int L, int R, const Info&v){
        if(L > R) return;
        __modify(root, left_margin, right_margin, L, R+1, [&](int rt,int l,int r){
            node[rt].info += v;
        });
    }
    void apply(int L, int R, const Tag&t){
        if(L > R) return;
        __modify(root, left_margin, right_margin, L, R+1, [&](int rt,int l,int r){
            __apply(rt, t);
        });
    }

    template<typename Filter>
    static Info __ask(int rt, int l, int r, int L, int R, const Filter&func){
        if(R<=l || r<=L || rt==-1) return Info();
        if(L<=l && r<=R) return func(rt, l, r);
        __push(rt);
        return __ask(node[rt].ch[0], l, (l+r)/2, L, R, func) +
            __ask(node[rt].ch[1], (l+r)/2, r, L, R, func);
    }
    Info ask(int L, int R){
        if(L > R) return Info();
        return __ask(root, left_margin, right_margin, L, R+1, [&](int rt,int l,int r){
            return node[rt].info;
        });
    }

    static const int npos = -1;
    template<class Predicator>
    static int __find(int rt, int l, int r, int L, int R, const Predicator&pred, bool forward){
        if(R<=l || r<=L || rt==-1 || !pred(node[rt].info)) return npos;
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
        return __find(root, left_margin, right_margin, L, R+1, pred, true);
    }
    template<class Predicator>
    int findBack(int L, int R, const Predicator&pred){
        return __find(root, left_margin, right_margin, L, R+1, pred, false);
    }
    
    static int __merge(int rt1, int rt2, int l, int r){
        if(rt1 == -1) return rt2;
        if(rt2 == -1) return rt1;
        if(l+1==r){
            node[rt1].info += node[rt2].info;
            return rt1;
        }
        __push(rt1); __push(rt2);
        node[rt1].ch[0] = __merge(node[rt1].ch[0], node[rt2].ch[0], l, (l+r)/2);
        node[rt1].ch[1] = __merge(node[rt1].ch[1], node[rt2].ch[1], (l+r)/2, r);
        __pull(rt1); __drop(rt2);
        return rt1;
    }
    void operator+=(const SegTree&o){
        __merge(root, o.root, left_margin, right_margin);
    }
};

template<class Info, class Tag>
int SegTree<Info, Tag>::left_margin;

template<class Info, class Tag>
int SegTree<Info, Tag>::right_margin;

template<class Info, class Tag>
std::vector<typename SegTree<Info, Tag>::Node> SegTree<Info, Tag>::node;

template<class Info, class Tag>
int SegTree<Info, Tag>::drop = -1;

#endif