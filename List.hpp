#ifndef LIST_HPP
#define LIST_HPP

template<typename T>
struct List{
    struct Node{
        T val;
        int pre, nxt;
    };
    struct iterator{
        using difference_type = std::ptrdiff_t;
        using value_type = T;
        int idx;
        T&operator*()const{
            return data[idx].val;
        }
        T*operator->()const{
            return &data[idx].val;
        }
        iterator&operator--(){
            idx = data[idx].pre;
            return *this;
        }
        iterator&operator++(){
            idx = data[idx].nxt;
            return *this;
        }
        iterator operator--(int){
            auto tmp = *this;
            --*this;
            return tmp;
        }
        iterator operator++(int){
            auto tmp = *this;
            ++*this;
            return tmp;
        }
        bool operator==(const iterator&o)const{ return idx == o.idx; }
        bool operator!=(const iterator&o)const{ return idx != o.idx; }
    };
    static std::vector<Node> data;
    static const int drop_list = 0;
    static void unlink(int pos){
        auto&cur = data[pos];
        data[cur.pre].nxt = cur.nxt;
        data[cur.nxt].pre = cur.pre;
    }
    static void link(int pre, int pos, int nxt){
        data[pos] = {data[pos].val, pre, nxt};
        data[pre].nxt = data[nxt].pre = pos;
    }
    static int new_node(){
        auto pos = data[drop_list].pre;
        if(pos == drop_list){
            pos = (int) data.size();
            data.emplace_back();
        }else unlink(pos);
        data[pos].pre = data[pos].nxt = pos;
        return pos;
    }
    static void del_node(int pos){
        if(pos < 0) return;
        unlink(pos);
        link(data[drop_list].pre, pos, drop_list);
    }

    iterator head;

    List() : head(new_node()) {}
    List(const List&o) : List() {
        for(auto&x : o) insert(end(), x);
    }
    List(const std::initializer_list<T>&ini) : List() {
        for(auto&&x : ini) insert(end(), x);
    }
    ~List() { del_node(head.idx); }
    iterator insert(iterator pos_itr, T val){
        auto new_pos = new_node();
        data[new_pos].val = val;
        link(std::prev(pos_itr).idx, new_pos, pos_itr.idx);
        return {new_pos};
    }
    iterator erase(iterator pos_itr){
        auto nxt = std::next(pos_itr).idx;
        del_node(pos_itr.idx);
        return {nxt};
    }
    iterator begin()const{
        return {data[head.idx].nxt};
    }
    iterator end()const{
        return head;
    }
    bool empty(){
        return data[head.idx].pre == data[head.idx].nxt;
    }
    void splice(iterator pos, List<T>&&o){
        int cur_out = std::prev(pos).idx, cur_in = pos.idx;
        int o_in = std::next(o.head).idx;
        int o_out = std::prev(o.head).idx;
        del_node(o.head.idx);
        data[cur_out].nxt = o_in;
        data[o_in].pre = cur_out;
        data[o_out].nxt = cur_in;
        data[cur_in].pre = o_out;
    }
};

template<typename T>
std::vector<typename List<T>::Node> List<T>::data(1, {0, 0, 0});

#endif