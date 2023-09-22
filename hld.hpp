#ifndef HLD_CPP
#define HLD_CPP
// 重链剖分

struct HLD{
    int n;
    std::vector<int> sz, top, dpt, fa, in, out, seq;
    std::vector<std::vector<int>> adj;
    int cur;
    
    HLD(){}
    HLD(int n){
        init(n);
    }
    void init(int n){
        this->n = n;
        sz.resize(n); top.resize(n); dpt.resize(n);
        fa.resize(n); in.resize(n); out.resize(n);
        seq.resize(n); adj.assign(n,{});
        cur = 0;
    }
    void add_edge(int u,int v){
        adj[u].push_back(v);
        adj[v].push_back(u);
    }
    void work(int rt=0){
        top[rt] = rt;
        dpt[rt] = 0;
        fa[rt] = -1;
        dfs1(rt);
        dfs2(rt);
    }
    void dfs1(int u){
        if(fa[u]!=-1){
            adj[u].erase(std::find(adj[u].begin(), adj[u].end(), fa[u]));
        }
        
        sz[u] = 1;
        for(auto&v:adj[u]){
            fa[v] = u;
            dpt[v] = dpt[u] + 1;
            dfs1(v);
            sz[u] += sz[v];
        }
        std::swap(*std::max_element(adj[u].begin(),adj[u].end(),[&](int x,int y){
            return sz[x]<sz[y];
        }), adj[u].front());
    }
    void dfs2(int u){
        in[u] = cur++;
        seq[in[u]] = u;
        for(auto v:adj[u]){
            top[v] = v==adj[u][0]?top[u]:v;
            dfs2(v);
        }
        out[u] = cur;
    }
    int lca(int u,int v){
        while(top[u]!=top[v]){
            if(dpt[top[u]]>dpt[top[v]]){
                u = fa[top[u]];
            }else{
                v = fa[top[v]];
            }
        }
        return dpt[u]<dpt[v]?u:v;
    }
    
    int dist(int u,int v){
        return dpt[u] + dpt[v] - dpt[lca(u, v)]*2;
    }
    
    int jump(int u,int k){
        if(dpt[u]<k) return -1;
        int d = dpt[u]-k;
        while(dpt[top[u]]>d){
            u = fa[top[u]];
        }
        return seq[in[u]-dpt[u]+d];
    }
    
    bool check_anc(int u,int v){
        return in[u]<=in[v] && in[v]<out[u];
    }
    
    int rooted_fa(int u,int rt){
        if(u==rt) return -1;
        if(!check_anc(u,rt)) return fa[u];
        auto it = std::upper_bound(adj[u].begin(),adj[u].end(),rt,[&](int x, int y){
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