#ifndef SEGTREE_ZKW_TAG_HPP
#define SEGTREE_ZKW_TAG_HPP

struct Tag{
    char ch{-1};

    void operator+=(const Tag&t){
        if(t.ch < 0) return;
        ch = t.ch;
    }
};

struct Info{
    std::array<int,26> cnt {};
    int size = 0;

    Info operator+(const Info&o){
        Info res;
        for(int i=0;i<26;i++)res.cnt[i] = cnt[i] + o.cnt[i];
        res.size = size + o.size;
        return res;
    }
    void operator+=(const Tag&t){
        if(t.ch < 0) return;
        cnt.fill(0);
        cnt[t.ch - 'a'] = size;
    }
};

template<class Info, class Tag>
struct SegTreeZKW{
    int n, lg, N;
    std::vector<Info> info;
    std::vector<Tag> tag;

    SegTreeZKW() = default;
    template<typename ...Args>
    explicit SegTreeZKW(Args&&... args){ init(std::forward<Args>(args)...); }

    void init(int n){
        this->n = n;
        this->lg = std::__lg(n) + 1;
        this->N = 1 << lg;
        info.assign(N << 1, Info{});
        tag.assign(N << 1, Tag{});
    }

    void init(const std::vector<Info>&arr){
        init(arr.size());
        std::copy(arr.begin(), arr.end(), info.begin() + N);
        for(int x=N - 1; x; x--){
            info[x] = info[x * 2] + info[x * 2 + 1];
        }
    }

    void push_down(int x) {
        info[x * 2] += tag[x]; tag[x * 2] += tag[x];
        info[x * 2 + 1] += tag[x]; tag[x * 2 + 1] += tag[x];
        tag[x] = Tag{};
    }

    void apply(int l, int r, const Tag&t){
        if(l > r) return;
        l += N, r += N;
        for(auto i=lg; i; i--) push_down(l >> i), push_down(r >> i);
        info[l] += t;
        if(l < r){
            info[r] += t;
            for(; l ^ r ^ 1; ){
                if(~l & 1) info[l ^ 1] += t, tag[l ^ 1] += t;
                if( r & 1) info[r ^ 1] += t, tag[r ^ 1] += t;
                l /= 2, r /= 2;
                info[l] = info[l * 2] + info[l * 2 + 1];
                info[r] = info[r * 2] + info[r * 2 + 1];
            }
        }
        for(; l /= 2; ) info[l] = info[l * 2] + info[l * 2 + 1];
    }

    Info ask(int l, int r){
        if(l > r) return {};
        l += N, r += N;
        for(auto i=lg; i; i--) push_down(l >> i), push_down(r >> i);
        Info infol = info[l];
        if(l < r){
            Info infor = l == r ? Info{} : info[r];
            for(; l ^ r ^ 1; ){
                if(~l & 1) infol = infol + info[l ^ 1];
                if( r & 1) infor = info[r ^ 1] + infor;
                l /= 2, r /= 2;
            }
            infol = infol + infor;
        }
        return infol;
    }
};

#endif