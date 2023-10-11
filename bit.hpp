#ifndef BIT_CPP
#define BIT_CPP
// 树状数组

int lowbit(int x){return x&-x;}

template<typename T=int>
struct BIT{
    int n;
    std::vector<T> node;
    BIT(int n=0){init(n);}
    BIT(const std::vector<T>&initarr){
        init(initarr);
    }
    void init(const std::vector<T>&initarr){
        init(initarr.size());
        for(int i=1;i<=n;i++){
            node[i] += initarr[i-1];
            if(i+lowbit(i)<=n)node[i+lowbit(i)] += node[i];
        }
    }
    void init(int n){
        this->n = n;
        node.assign(n+1,T());
    }
    void add(int k,T v){
        if(!k)return;
        for(;k<=n;k+=lowbit(k))node[k]+=v;
    }
    T __ask(int k){
        T ret = 0;
        for(;k>0;k-=lowbit(k))ret+=node[k];
        return ret;
    }
    T ask(int l,int r){
        if(l>r)std::swap(l,r);
        return query(r) - query(l-1);
    }
    T ask(int p){return ask(p,p);}
    int kth(T k){
        int x = 0;
        for(int i=1<<std::__lg(n);i;i>>=2){
            if(x+i<=n&&k>=node[x+i-1]){
                x += i;
                k -= node[x-1];
            }
        }
        return x;
    }
};

template<typename T=int>
struct BIT_range{
    int n;
    BIT<T> di,dn;
    BIT_range(int n=0):di(n),dn(n){}
    BIT_range(const std::vector<T>&initarr){
        init(initarr);
    }
    void init(const std::vector<T>&initarr){
        init(initarr.size());
        for(int i=1;i<=n;i++){
            di.node[i] += i*initarr[i-1];
            if(i<n)di.node[i+1] -= (i+1)*initarr[i-1];
            if(i+lowbit(i)<=n)di.node[i+lowbit(i)] += di.node[i];
            dn.node[i] += initarr[i-1];
            if(i<n)dn.node[i+1] -= initarr[i-1];
            if(i+lowbit(i)<=n)dn.node[i+lowbit(i)] += dn.node[i];
        }
    }
    void init(int n){
        this->n = n;
        di.init(n);dn.init(n);
    }
    void __add(int k,int v){
        di.add(k,k*v);
        dn.add(k,v);
    }
    void add(int l,int r,T v){
        if(l>r)std::swap(l,r);
        __add(l,v);__add(r+1,-v);
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