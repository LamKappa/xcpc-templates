#ifndef EDGE_BICONNECTED_COMPONENT
#define EDGE_BICONNECTED_COMPONENT

std::set<std::pair<int, int>> E;

struct EBCC{
    int n, cur, cnt;
    std::vector<std::vector<int>> adj;
    std::vector<int> stk;
    std::vector<int> dfn, low, bel;
    
    EBCC() = default;
    EBCC(int n){init(n);}
    
    void init(int n){
        this->n = n;
        adj.assign(n, {});
        dfn.assign(n, -1);
        low.resize(n);
        bel.assign(n, -1);
        stk.clear();
        cur = cnt = 0;
    }
    
    void add_edge(int u, int v){
        adj[u].push_back(v);
        adj[v].push_back(u);
    }
    
    void dfs(int u, int fa){
        dfn[u] = low[u] = cur++;
        stk.push_back(u);
        
        for (auto v : adj[u]){
            if(v == fa){
                fa = -1;
                continue;
            }
            if(dfn[v] == -1){
                E.emplace(u, v);
                dfs(v, u);
                low[u] = std::min(low[u], low[v]);
            }else if(bel[v] == -1 && dfn[v] < dfn[u]){
                E.emplace(u, v);
                low[u] = std::min(low[u], dfn[v]);
            }
        }
        
        if(dfn[u] == low[u]){
            for(int v=-1; v!=u; stk.pop_back()){
                v = stk.back();
                bel[v] = cnt;
            }
            cnt++;
        }
    }
    
    std::vector<int> work(){
        for(int i=0; i<n; i++){
            if(dfn[i] != -1) continue;
            dfs(i, -1);
        }
        return bel;
    }
};

#endif