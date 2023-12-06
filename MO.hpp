#ifndef MO_ALG
#define MO_ALG

struct Info{
    int ans = 0;
    std::vector<int> info, cnt;
    Info(int n) : info(n), cnt(1e6+1, 0){}
};

struct MO : public Info{
    int N = 1, M = 1, B;
    int l, r, t, id;

    struct Query{
        int l, r, t, id;
    };
    std::vector<Query> query;

    struct Change{
        int x, v;
    };
    std::vector<Change> change;


    MO(int n) : Info(n){};

    void add_query(Query&&q){
        q.t = change.size()-1;
        q.id = query.size();
        query.emplace_back(q);
    }

    void add_change(Change&&c){
        change.emplace_back(c);
    }

    void init(){
        B = std::pow(std::pow(info.size(), 2) * (change.size()+1) / (query.size()+1), 1./3.);
        std::sort(query.begin(), query.end(), [this](auto&a,auto&b){
            return std::make_tuple(a.l/B, a.r/B, a.t/B) < 
                std::make_tuple(b.l/B, b.r/B, b.t/B);
        });
        r = t = -1, l = id = 0;
    }

    void redo(int i){
        if(0==cnt[info[i]]++){
            ans++;
        }
    }

    void undo(int i){
        if(0==--cnt[info[i]]){
            ans--;
        }
    }

    void update(int i){
        auto&[x, v] = change[i];
        if(l<=x && x<=r) undo(x);
        std::swap(v, info[x]);
        if(l<=x && x<=r) redo(x);
    }

    Query next(){
        auto x = query[id++];
        while(l < x.l) undo(l++);
        while(l > x.l) redo(--l);
        while(r < x.r) redo(++r);
        while(r > x.r) undo(r--);
        while(t < x.t) update(++t);
        while(t > x.t) update(t--);
        return x;
    }

    bool finish(){
        return id == query.size();
    }
};

#endif