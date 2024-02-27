#ifndef SEGTREE_ZKW
#define SEGTREE_ZKW
// ZKW (+, max)

using i64 = long long;
constexpr i64 INF = 1e18;

struct SegTree {
    int N;
    std::vector<i64> info, tag;
    SegTree(int n) {
        N = 1 << __lg(n+5) + 1;
        info.assign(N<<1, 0);
        tag.assign(N<<1, 0);
    }
    SegTree(const std::vector<int>&initarr):SegTree(initarr.size()) {
        for(int i=1; i<=initarr.size(); i++) { info[i + N] = initarr[i-1]; }
        for(int i=N; --i;) { info[i] = info[i << 1] + info[i << 1 | 1]; }
    }
    
    void update(int l, int r, i64 d) {
        for(l += N - 1, r += N + 1; l ^ r ^ 1; l >>= 1, r >>= 1) {
            if(l < N) {
                info[l] = std::max(info[l << 1], info[l << 1 | 1]) + tag[l];
                info[r] = std::max(info[r << 1], info[r << 1 | 1]) + tag[r];
            }
            if(~l & 1) { info[l ^ 1] += d, tag[l ^ 1] += d; }
            if(r & 1) { info[r ^ 1] += d, tag[r ^ 1] += d; }
        }
        for(; l; l >>= 1, r >>= 1) {
            if(l < N) {
                info[l] = std::max(info[l << 1], info[l << 1 | 1]) + tag[l];
                info[r] = std::max(info[r << 1], info[r << 1 | 1]) + tag[r];
            }
        }
    };
    i64 query(int l, int r) {
        i64 maxl = -INF, maxr = -INF;
        for(l += N - 1, r += N + 1; l ^ r ^ 1; l >>= 1, r >>= 1) {
            maxl += tag[l], maxr += tag[r];
            if(~l & 1) { maxl = std::max(maxl, info[l ^ 1]); }
            if(r & 1) { maxr = std::max(maxr, info[r ^ 1]); }
        }
        for(; l; l >>= 1, r >>= 1) {
            maxl += tag[l], maxr += tag[r];
        }
        return std::max(maxl, maxr);
    };
};

#endif