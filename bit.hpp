#ifndef BIT_CPP
#define BIT_CPP
// 树状数组

template<typename T=int>
struct BIT{
    int _log;
    std::vector<T> node;
    BIT(int N=16):_log(__lg(N-1)+1),node((1ll<<_log)+1,0){}
    static int lowbit(int x){return x&-x;}
    void add(int k,T v){
        if(!k)return;
        for(;k<node.size();k+=lowbit(k))node[k]+=v;
    }
    T ask(int k){
        T ret = 0;
        for(;k>0;k-=lowbit(k))ret+=node[k];
        return ret;
    }
    T ask(int l,int r){
        if(l>r)std::swap(l,r);
        return ask(r) - ask(l-1);
    }
    int kth(T k){
        int pos = 1<<_log;
        for(int i=pos>>1;pos-i<node.size()&&i;i>>=1)
            if(node[pos-i]>=k) pos-=i;
            else k-=node[pos-i];
        return pos;
    }
};

template<typename T=int>
struct BIT_range{
    BIT<T> di,dn;
    BIT_range(int N=16):di(N),dn(N){}
    void add_d(int k,int v){
        di.add(k,k*v);
        dn.add(k,v);
    }
    void add(int l,int r,T v){
        if(l>r)std::swap(l,r);
        add_d(l,v);add_d(r+1,-v);
    }
    void add(int k,T v){add(k,k,v);}
    T ask(int k){
        return (k+1)*dn.ask(k) - di.ask(k);
    }
    T ask(int l,int r){
        if(l>r)std::swap(l,r);
        return ask(r) - ask(l-1);
    }
};

#endif