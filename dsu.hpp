#ifndef DSU_CPP
#define DSU_CPP
// 并查集

struct DSU{
    std::vector<int> fa, sz;
    DSU(){}
    DSU(int n){init(n);}
    void init(int n){
        fa.resize(n);
        std::iota(fa.begin(), fa.end(), 0);
        sz.assign(n, 1);
    }
    int find(int x){
        while(x!=fa[x]){
            x = fa[x] = fa[fa[x]];
        }
        return x;
    }
    bool same(int x, int y){
        return find(x) == find(y);
    }
    bool merge(int x, int y){
        x = find(x);y = find(y);
        if(x==y)return false;
        if(sz[x]<sz[y])swap(x,y);
        sz[x] += sz[y];
        fa[y] = x;
        return true;
    }
    int get_size(int x){
        return sz[find(x)];
    }
};

#endif