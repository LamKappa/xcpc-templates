#ifndef STRONGLY_CONNECTED_COMPONENT
#define STRONGLY_CONNECTED_COMPONENT

struct SCC{
    int n, cur, cnt;
    std::vector<std::vector<int>> adj;
    std::vector<int> stk;
    std::vector<int> dfn, low, bel;
    
    SCC() = default;
    SCC(int n){init(n);}
    
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
    }
    
    void dfs(int u){
        dfn[u] = low[u] = cur++;
        stk.push_back(u);
        
        for(auto v : adj[u]){
            if(dfn[v] == -1){
                dfs(v);
                low[u] = std::min(low[u], low[v]);
            }else if(bel[v] == -1){
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
            dfs(i);
        }
        return bel;
    }
};

#endif