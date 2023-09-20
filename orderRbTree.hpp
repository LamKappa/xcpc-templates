#ifndef ORDERRBTREE_CPP
#define ORDERRBTREE_CPP
// 带序红黑树

#include<ext/pb_ds/assoc_container.hpp>
#include<ext/pb_ds/tree_policy.hpp>

template<typename T>
using order_rbtree = __gnu_pbds::tree<
T,
__gnu_pbds::null_type,
std::less<T>,
__gnu_pbds::rb_tree_tag,
__gnu_pbds::tree_order_statistics_node_update>;

#endif