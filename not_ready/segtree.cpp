#include <bits/stdc++.h>

using i64 = long long;
std::vector<int> minp, primes, phi;

void sieve(int n) {
    minp.assign(n + 1, 0);
    phi.assign(n + 1, 0);
    primes.clear();
    phi[1] = 1;
    
    for (int i = 2; i <= n; i++) {
        if (minp[i] == 0) {
            minp[i] = i;
            phi[i] = i - 1;
            primes.push_back(i);
        }
        
        for (auto p : primes) {
            if (i * p > n) {
                break;
            }
            minp[i * p] = p;
            if (p == minp[i]) {
                phi[i * p] = phi[i] * p;
                break;
            }
            phi[i * p] = phi[i] * (p - 1);
        }
    }
}

template<class Info,
    class Merge = std::plus<Info>>
struct SegmentTree {
    const int n;
    const Merge merge;
    std::vector<Info> info;
    SegmentTree(int n) : n(n), merge(Merge()), info(4 << std::__lg(n)) {}
    SegmentTree(std::vector<Info> init) : SegmentTree(init.size()) {
        std::function<void(int, int, int)> build = [&](int p, int l, int r) {
            if (r - l == 1) {
                info[p] = init[l];
                return;
            }
            int m = (l + r) / 2;
            build(2 * p, l, m);
            build(2 * p + 1, m, r);
            pull(p);
        };
        build(1, 0, n);
    }
    void pull(int p) {
        info[p] = merge(info[2 * p], info[2 * p + 1]);
    }
    void modify(int p, int l, int r, int x, const Info &v) {
        if (r - l == 1) {
            info[p] = v;
            return;
        }
        int m = (l + r) / 2;
        if (x < m) {
            modify(2 * p, l, m, x, v);
        } else {
            modify(2 * p + 1, m, r, x, v);
        }
        pull(p);
    }
    void modify(int p, const Info &v) {
        modify(1, 0, n, p, v);
    }
    void change(int p, int l, int r, int x, int v) {
        // std::cerr << "v : " << v << " " << info[p].a.back()[0] << "\n";
        if (info[p].a == v) {
            info[p].a = phi[info[p].a];
            info[p].sum += info[p].cnt;
        }
        info[p].sum -= 1;
        if (r - l == 1) {
            return;
        }
        int m = (l + r) / 2;
        if (x < m) {
            change(2 * p, l, m, x, v);
        } else {
            change(2 * p + 1, m, r, x, v);
        }
        // pull(p);
    }
    void change(int p, int v) {
        change(1, 0, n, p, v);
    }
    Info rangeQuery(int p, int l, int r, int x, int y) {
        if (l >= y || r <= x) {
            return Info();
        }
        if (l >= x && r <= y) {
            return info[p];
        }
        int m = (l + r) / 2;
        return merge(rangeQuery(2 * p, l, m, x, y), rangeQuery(2 * p + 1, m, r, x, y));
    }
    Info rangeQuery(int l, int r) {
        return rangeQuery(1, 0, n, l, r);
    }
    template<class F>
    int findFirst(int p, int l, int r, int x, int y, F pred) {
        if (l >= y || r <= x || !pred(info[p])) {
            return -1;
        }
        if (r - l == 1) {
            return l;
        }
        int m = (l + r) / 2;
        int res = findFirst(2 * p, l, m, x, y, pred);
        if (res == -1) {
            res = findFirst(2 * p + 1, m, r, x, y, pred);
        }
        return res;
    }
    template<class F>
    int findFirst(int l, int r, F pred) {
        return findFirst(1, 0, n, l, r, pred);
    }
};

struct Info {
    int a;
    int cnt;
    int sum;
    
    Info() : a{} {}
    Info(int x) : a{x}, cnt{1}, sum{0} {}
};

Info operator+(Info a, Info b) {
    if (!a.a) {
        return b;
    }
    if (!b.a) {
        return a;
    }
    Info c;
    while (a.a != b.a) {
        if (a.a > b.a) {
            a.sum += a.cnt;
            a.a = phi[a.a];
        } else {
            b.sum += b.cnt;
            b.a = phi[b.a];
        }
    }
    c.a = a.a;
    c.cnt = a.cnt + b.cnt;
    c.sum = a.sum + b.sum;
    return c;
}
struct DSU {
    std::vector<int> f, siz;
    
    DSU() {}
    DSU(int n) {
        init(n);
    }
    
    void init(int n) {
        f.resize(n);
        std::iota(f.begin(), f.end(), 0);
        siz.assign(n, 1);
    }
    
    int leader(int x) {
        while (x != f[x]) {
            x = f[x] = f[f[x]];
        }
        return x;
    }
    
    bool same(int x, int y) {
        return leader(x) == leader(y);
    }
    
    bool merge(int x, int y) {
        x = leader(x);
        y = leader(y);
        if (x == y) {
            return false;
        }
        siz[x] += siz[y];
        f[y] = x;
        return true;
    }
    
    int size(int x) {
        return siz[leader(x)];
    }
};
int main() {
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);
    
    sieve(5E6);
    
    int n, m;
    std::cin >> n >> m;
    
    std::vector<int> a(n);
    for (int i = 0; i < n; i++) {
        std::cin >> a[i];
    }
    
    SegmentTree seg(std::vector<Info>(a.begin(), a.end()));
    
    DSU dsu(n + 1);
    
    while (m--) {
        int t, l, r;
        std::cin >> t >> l >> r;
        
        l--;
        if (t == 1) {
            for (int i = dsu.leader(l); i < r; i = dsu.leader(i + 1)) {
                if (a[i] == 1) {
                    dsu.merge(i + 1, i);
                } else {
                    seg.change(i, a[i]);
                    a[i] = phi[a[i]];
                }
            }
        } else {
            auto info = seg.rangeQuery(l, r);
            std::cout << info.sum << "\n";
        }
    }
    
    return 0;
}