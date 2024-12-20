#ifndef DELAUNAY_HPP
#define DELAUNAY_HPP

namespace Delaunay{
    using namespace Geometry;
    bool inside_Circular(const std::array<Point, 3>&cir, const Point&d){
        auto&[a, b, c] = cir;
        if(sign(cross(b - a, c - a)) < 0) return inside_Circular({a, c, b}, d);
        using P3 = std::array<f80, 3>;
        P3 ba = {b.x - a.x, b.y - a.y, b.square() - a.square()};
        P3 ca = {c.x - a.x, c.y - a.y, c.square() - a.square()};
        P3 da = {d.x - a.x, d.y - a.y, d.square() - a.square()};
        P3 abc = {
                ba[1] * ca[2] - ba[2] * ca[1],
                ba[2] * ca[0] - ba[0] * ca[2],
                ba[0] * ca[1] - ba[1] * ca[0]
        };
        return sign(abc[0] * da[0] + abc[1] * da[1] + abc[2] * da[2]) < 0;
    }

    struct Graph{
        template<typename T>
        using List = std::list<T>;
        struct Edge{
            int tar;
            List<Edge>::iterator rev;
        };
        std::vector<int> id;
        std::vector<List<Edge>> adj;

        auto add_edge(int u, List<Edge>::iterator posu, int v, List<Edge>::iterator posv){
            posu = adj[u].insert(posu, Edge{v});
            posv = adj[v].insert(posv, Edge{u});
            posu->rev = posv; posv->rev = posu;
            return std::make_pair(posu, posv);
        }

        void work(const Points&pts){
            id.resize(pts.size());
            std::iota(id.begin(), id.end(), 0);
            std::sort(id.begin(), id.end(), [&](auto i, auto j){
                return pts[i] < pts[j];
            });
            adj.resize(pts.size());

            auto rec = [&](auto&&rec, int l, int r)->List<int>{
                if(l + 1 == r) return {l};
                if(l + 2 == r){
                    add_edge(id[l], adj[id[l]].begin(), id[l+1], adj[id[l+1]].end());
                    return {l, l + 1};
                }
                int m = (l + r) / 2;
                List<int> I;
                I.splice(I.end(), rec(rec, l, m));
                auto it = std::prev(I.end());
                I.splice(I.end(), rec(rec, m, r));

                int pl, pr;
                for(it=std::next(it); it!=I.end(); it++){
                    auto it1 = std::prev(it);
                    auto it2 = std::prev(it1);
                    while(it1 != I.end() && it2 != I.end() &&
                        Line{pts[id[*it2]], pts[id[*it1]]}.onLeft(pts[id[*it]]) < 0){
                        I.erase(it1);
                        std::tie(it1, it2) = std::make_pair(it2, std::prev(it2));
                    }
                    if(*it1 < m && m <= *it) std::tie(pl, pr) = std::make_pair(id[*it1], id[*it]);
                    if(*it1 >= m) break;
                }
                auto[itl, itr] = add_edge(pl, adj[pl].begin(), pr, adj[pr].end());

                for(;;){
                    itl = std::next(itl); itr = std::prev(itr);
                    auto L = Line{pts[pl], pts[pr]};
                    int npl = -1, npr = -1;

                    while(itl != adj[pl].end()){
                        if(L.onLeft(pts[itl->tar]) <= 0) break;
                        auto nxt = std::next(itl);
                        if(nxt != adj[pl].end() &&
                            inside_Circular({pts[pl], pts[pr], pts[itl->tar]}, pts[nxt->tar])){
                            adj[itl->tar].erase(itl->rev);
                            itl = adj[pl].erase(itl);
                        }else{
                            npl = itl->tar;
                            break;
                        }
                    }

                    while(itr != adj[pr].end()){
                        if(L.onLeft(pts[itr->tar]) <= 0) break;
                        auto pre = std::prev(itr);
                        if(pre != adj[pr].end() &&
                            L.onLeft(pts[pre->tar]) > 0 &&
                            inside_Circular({pts[pl], pts[pr], pts[itr->tar]}, pts[pre->tar])){
                            adj[itr->tar].erase(itr->rev);
                            itr = std::prev(adj[pl].erase(itr));
                        }else{
                            npr = itr->tar;
                            break;
                        }
                    }

                    if(npl == -1 && npr == -1) break;
                    if(npr == -1 || (npl != -1 && !inside_Circular({pts[pl], pts[pr], pts[npl]}, pts[npr]))){
                        itl = std::next(itl->rev);
                        if(itl == adj[npl].end() && pts[pl].x <= pts[npl].x) itl = adj[npl].begin();
                        pl = npl;
                        std::tie(itl, itr) = add_edge(pl, itl, pr, std::next(itr));
                    }else{
                        itr = itr->rev;
                        if(itr == adj[npr].begin() && pts[pr].x > pts[npr].x) itr = adj[npr].end();
                        pr = npr;
                        std::tie(itl, itr) = add_edge(pl, itl, pr, itr);
                    }
                }

                return I;
            };

            rec(rec, 0, (int)pts.size());
        }
    };
}

#endif