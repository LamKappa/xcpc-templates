#ifndef LCSEGTREE_CPP
#define LCSEGTREE_CPP
// 李超线段树

#define maxn 200005

template<typename Line>
struct LCTree{
    int root, n, tot;
    int ls[maxn << 3], rs[maxn << 3];
    Line line[maxn << 3];
    std::function<bool(Line x, Line y, int t)> cmp = [](Line x, Line y, int t) -> bool{return true;};
    
    LCTree(int n, std::function<bool(Line x, Line y, int t)> cmp){
        this -> n = n;
        this -> tot = 0;
        this -> cmp = cmp;
        clear();
    }
    
    int newNode(Line v){
        line[++tot] = v;
        ls[tot] = rs[tot] = 0;
        return tot;
    }
    
    int insert(Line t, int now, int l, int r){
        if(!now) return newNode(t);
        int mid = (l + r) >> 1;
        bool lc2 = cmp(line[now], t, l), rc2 = cmp(line[now], t, r);
        bool lc1 = cmp(t, line[now], l), rc1 = cmp(t, line[now], r);
        if(!lc1 && !rc1){return now;}
        if(!lc2 && !rc2){line[now] = t; return now;}
        if(lc1){
            if(cmp(t, line[now], mid)) rs[now] = insert(line[now], rs[now], mid + 1, r), line[now] = t;
            else ls[now]=insert(t, ls[now], l, mid);
        }
        else if(rc1){
            if(cmp(t, line[now], mid + 1)) ls[now] = insert(line[now], ls[now], l, mid), line[now] = t;
            else rs[now] = insert(t, rs[now], mid + 1, r);
        }
        return now;
    }
    
    void insert(Line t){
        root = insert(t, root, 1, n);
    }
    
    Line better(Line x, Line y, int t){
        return cmp(x, y, t) ? x : y;
    }
    
    Line ask(int x, int now, int l, int r){
        if(!now) return Line();
        int mid = (l + r) >> 1;
        if(x <= mid) return ls[now] ? better(ask(x, ls[now], l, mid), line[now], x) : line[now];
        else return rs[now] ? better(ask(x, rs[now], mid + 1, r), line[now], x) : line[now];
    }
    
    Line ask(int x){
        return ask(x, root, 1, n);
    }
    
    void clear(void){
        root = 0, tot = 0;
    }
};

#endif