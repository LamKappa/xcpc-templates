#ifndef HLD
#define HLD
// 重链剖分

// depends SegTree.hpp

template<class Info, class Tag>
struct HLD{
    int n;
    std::vector<int> sz, top, dpt, fa, in, out, seq;
    std::vector<std::vector<int>> adj;
    SegTree<Info, Tag> tree;
    
    HLD(int n=0){init(n);}
    void init(int n){
        this->n = n;
        sz.resize(n); top.resize(n); dpt.resize(n);
        fa.resize(n); in.resize(n); out.resize(n);
        seq.resize(n); adj.assign(n,{});
    }
    void init_val(const std::vector<Info>&initarr){
        std::vector<Info> mapping(n);
        for(int i=0;i<n;i++){
            mapping[in[i]] = initarr[i];
        }
        tree.init(mapping);
    }
    void add_edge(int u,int v){
        adj[u].push_back(v);
        adj[v].push_back(u);
    }
    void work(int rt=0){
        top[rt] = rt; fa[rt] = -1;
        int cur = dpt[rt] = 0;
        std::function<void(int)> dfs1 = [&](int u){
            sz[u] = 1;
            for(auto v : adj[u]){
                if(fa[u] == v) continue;
                fa[v] = u;
                dpt[v] = dpt[u] + 1;
                dfs1(v);
                sz[u] += sz[v];
            }
            if(adj[u].empty()) return;
            std::swap(*std::max_element(adj[u].begin(), adj[u].end(), [&](int x,int y){
                return fa[u]==x || sz[x]<sz[y];
            }), adj[u].front());
        };
        std::function<void(int)> dfs2 = [&](int u){
            in[u] = cur++;
            seq[in[u]] = u;
            for(auto v : adj[u]){
                if(fa[u] == v) continue;
                top[v] = (v==adj[u][0] ? top[u] : v);
                dfs2(v);
            }
            out[u] = cur;
        };
        dfs1(rt); dfs2(rt);
    }

    void apply(int u, int v, const Tag&t){
        while(top[u]!=top[v]){
            if(dpt[top[u]] < dpt[top[v]]) std::swap(u, v);
            tree.apply(in[top[u]], in[u], t);
            u = fa[top[u]];
        }
        tree.apply(in[u], in[v], t);
    }

    void apply_sub(int rt, const Tag&t){
        tree.apply(in[rt], out[rt]-1, t);
    }
    
    Info ask(int u, int v){
        Info res = Info();
        while(top[u]!=top[v]){
            if(dpt[top[u]] < dpt[top[v]]) std::swap(u, v);
            res = res + tree.ask(in[top[u]], in[u]);
            u = fa[top[u]];
        }
        return res + tree.ask(in[u], in[v]);
    }

    Info ask_sub(int rt){
        return tree.ask(in[rt], out[rt]-1);
    }

    int lca(int u, int v){
        while(top[u] != top[v]){
            if(dpt[top[u]] < dpt[top[v]]) std::swap(u, v);
            u = fa[top[u]];
        }
        return dpt[u]<dpt[v] ? u : v;
    }
    
    int dist(int u, int v){
        return dpt[u] + dpt[v] - dpt[lca(u, v)]*2;
    }

    int jump(int u, int k){
        if(dpt[u] < k) return -1;
        int d = dpt[u]-k;
        while(dpt[top[u]] > d){
            u = fa[top[u]];
        }
        return seq[in[u] - dpt[u]+d];
    }
    
    bool check_anc(int u, int v){
        return in[u]<=in[v] && in[v]<out[u];
    }
    
    int rooted_fa(int u, int rt){
        if(u == rt) return -1;
        if(!check_anc(u, rt)) return fa[u];
        auto it = std::upper_bound(adj[u].begin(), adj[u].end(), rt, [&](int x, int y){
            return in[x]<in[y];
        }) - 1;
        return *it;
    }
    
    int rooted_size(int u, int rt){
        if(u == rt) return n;
        if(!check_anc(u,rt)) return sz[u];
        return n - sz[rooted_fa(u, rt)];
    }
    
    int rooted_lca(int a, int b, int rt){
        return lca(a, b) ^ lca(b, rt) ^ lca(rt, a);
    }
};

#endif