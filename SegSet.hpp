#ifndef SEGSET_HPP
#define SEGSEG_HPP

using Segment = std::array<int, 2>;
constexpr int INF = 0x7fffffff;
struct SegSet : std::set<Segment>{

    SegSet(){
        std::set<Segment>::insert({-INF, -INF});
        std::set<Segment>::insert({INF, INF});
    }
    void insert(Segment sg){
        if(sg[1] <= -INF || sg[0] >= INF) return;
        auto itr = std::prev(upper_bound({sg[0], -INF}));
        if((*itr)[1] >= sg[0]){
            sg = {(*itr)[0], std::max(sg[1], (*itr)[1])};
            itr = erase(itr);
        }else itr = std::next(itr);
        while((*itr)[0] <= sg[1]){
            sg[1] = std::max(sg[1], (*itr)[1]);
            itr = erase(itr);
        }
        std::set<Segment>::insert(sg);
    }
    void insert(int l, int r){
        insert(Segment{l, r});
    }
    void insert(int x){
        insert(x, x + 1);
    }

    void merge(SegSet&t){
        auto&s = *this;
        if(s.size() < t.size()) std::swap(s, t);
        for(auto&sg : t){
            s.insert(sg);
        }
        t.clear();
    }
};

#endif