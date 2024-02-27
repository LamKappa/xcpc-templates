#ifndef DISJOINT_SET_UNION
#define DISJOINT_SET_UNION
// 并查集

struct DSU {
    std::vector<int> fa, sz, st;
    DSU() {}
    DSU(int n) {init(n);}
    void init(int n) {
        fa.resize(n);
        std::iota(fa.begin(), fa.end(), 0);
        sz.assign(n, 1);
        st.clear();
    }
    int find(int x) {
        while(x!=fa[x]) {
            x = fa[x] = fa[fa[x]];
        }
        return x;
    }
    bool same(int x, int y) {
        return find(x) == find(y);
    }
    bool merge(int x, int y) {
        x = find(x); y = find(y);
        if(x==y) { return false; }
        if(sz[x]<sz[y]) { std::swap(x,y); }
        sz[x] += sz[y];
        st.push_back(y);
        fa[y] = x;
        return true;
    }
    int get_size(int x) {
        return sz[find(x)];
    }
    void clear() {
        while(!st.empty()) {
            int x = st.back(), p = find(x);
            fa[p] = p; sz[p] = 1;
            fa[x] = x; sz[x] = 1;
            st.pop_back();
        }
    }
};

#endif