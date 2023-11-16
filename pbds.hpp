#ifndef STL_PBDS
#define STL_PBDS

#include <bits/extc++.h>
// using namespace __gnu_cxx;
// using namespace __gnu_pbds;

// 平衡树
template<typename K_T,
    typename Cmp = std::less<K_T>,
    typename V_T = __gnu_pbds::null_type
>
using order_set = __gnu_pbds::tree<
    K_T,                                                // key_type
    V_T,                                                // value_type
    Cmp,                                                // comparator
    __gnu_pbds::rb_tree_tag,                            // tag
    __gnu_pbds::tree_order_statistics_node_update       // policy
>;
template<typename K_T,
    typename V_T = __gnu_pbds::null_type,
    typename Cmp = std::less<K_T>
>
using order_map = order_set<K_T, Cmp, V_T>;
template<typename K_T, typename V_T, typename Cmp>
V_T& operator%(order_map<K_T,V_T,Cmp>&mp, const K_T&x){
    return mp.lower_bound(x)->second;
}
/*
tag:
rb_tree_tag                                 红黑树
splay_tree_tag                              Slpay
ov_tree_tag                                 向量树

Itr ::point_iterator

std::pair<point_iterator, bool> insert(T)   插入
bool erase(T/Itr)                           删除元素/迭代器
int order_of_key(T)                         返回排名
Itr find_by_order(T)                        排名对应元素
Itr lower_bound(T)/upper_bound(T)           二分查找
void join(order_set)                        合并
void split(T,order_set)                     保留<=,其余分离覆盖到order_set中
bool empty()                                判空
size_t size()                               大小
Itr begin()/end()                           首尾迭代器
*/

/***************/

// 堆
template<typename T,
    typename Cmp = std::greater<T>
>
using heap = __gnu_pbds::priority_queue<
    T,                                                  // type
    Cmp,                                                // comparator
    __gnu_pbds::pairing_heap_tag                        // tag
>;
/*
tag:
pairing_heap_tag        配对堆
thin_heap_tag           斐波那契堆
binomial_heap_tag       二项堆
binary_heap_tag         二叉堆

Itr ::point_iterator    可以指定为nullptr

usage:
Itr push(T)             入堆
void pop()              出堆
T top()                 堆顶
void modify(Itr, T)     修改
void join(heap)         合并,清空heap
bool empty()            判空
size_t size()           大小
void clear()            清空
*/

/***************/

// 哈希表
const int RANDOM = std::chrono::high_resolution_clock::now().time_since_epoch().count();
template<typename K_T>
struct Chash{
    int operator()(K_T x)const{return std::hash(x)^RANDOM;}
};
template<typename K_T, typename V_T, typename Hash = Chash<K_T>>
using hash_table = __gnu_pbds::cc_hash_table<K_T, V_T, Hash>;

/*
tag:
cc_hash_table           拉链法
gp_hash_table           二次探测法

V_T& operaotr[](K_T)    映射
*/


// 字典树
using trie = __gnu_pbds::trie<
    std::string,                                    //
    __gnu_pbds::null_type,                          //
    __gnu_pbds::trie_string_access_traits<>,        //
    __gnu_pbds::pat_trie_tag,                       // tag
    __gnu_pbds::trie_prefix_search_node_update      // policy
>;
/*
Itr insert(string)                          插入
void erase(string)                          删除
void join(trie)                             合并trie
std::pair<Itr, Itr> prefix_range(string)    前缀遍历[beign,end)
*/

#endif