#ifndef EULER_LOOP
#define EULER_LOOP

struct EulerLoop{
    std::vector<std::unordered_map<int, int>> adj;

    EulerLoop(int n){
        adj.assign(n, {});
    }

    void add_edge(int u, int v){
        adj[u][v]++;
    }

    std::vector<int> find(int rt){
        std::vector<int> ans;
        // auto adj = this->adj;
        auto dfs = [&](auto&&dfs, int u)->void{
            while(!adj[u].empty()){
                auto v = adj[u].begin()->first;
                if(0 == --adj[u][v]){
                    adj[u].erase(v);
                }
                dfs(dfs, v);
            }
            ans.push_back(u);
        };
        dfs(dfs, rt);
        reverse(ans.begin(), ans.end());
        return ans;
    }
    
    bool empty()const{
        bool fl = true;
        for(const auto&x : adj){
            fl &= x.empty();
        }
        return fl;
    }
};

#endif