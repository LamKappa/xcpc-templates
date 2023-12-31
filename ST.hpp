#ifndef SEPARATE_TABLE
#define SEPARATE_TABLE

template<typename T>
struct ST{
    unsigned N, B;
    T IE;
    std::function<T(T,T)> calc;
    std::vector<std::vector<T>> a, b;
    std::vector<T> pre, suf;

    ST(){}
    template<typename ... Args>
    ST(Args ... args){
        init(args ...);
    }
    template<typename ... Args>
    T Calc(const T&x, const T&y, Args ... args){
        return Calc(calc(x, y), args...);
    }
    T Calc(const T&x, const T&y){return calc(x, y);}
    void init(const std::vector<T>&v, 
        const std::function<T(T,T)>&calc = [](T x,T y){return std::max(x, y);},
        const T&IE = T()){
        this->N = v.size();
        if(N <= 0) return;
        this->B = sqrt(N) + 1;
        this->calc = calc;
        this->IE = IE;
        pre = suf = v;
        const int M = (N + B - 1) / B;
        const int lgM = std::__lg(M);
        const int lgB = std::__lg(B);
        a.resize(lgM+1); a[0].resize(M);
        b.resize(lgB+1); b[0] = v;
        for(int i=0; i<M; i++){
            a[0][i] = v[i * B];
            for(int j=1; j<B && i*B + j < N; j++){
                a[0][i] = Calc(a[0][i], v[i*B + j]);
            }
        }
        for(int j=0; j<lgM; j++){
            a[j+1].resize(M - (2<<j) + 1);
            for(int i=0; i+(2<<j)<=M; i++){
                a[j+1][i] = Calc(a[j][i], a[j][i + (1<<j)]);
            }
        }
        for(int j=0; j<lgB && (2<<j)<=N; j++){
            b[j+1].resize(N - (2<<j) + 1);
            for(int i=0; i+(2<<j)<=N; i++){
                b[j+1][i] = b[j][i];
                if((i+(1<<j))/B == i/B){
                    b[j+1][i] = Calc(b[j+1][i], b[j][i + (1<<j)]);
                }
            }
        }
        for(int i=1; i<N; i++){
            if(i%B != 0){
                pre[i] = Calc(pre[i-1], pre[i]);
            }
        }
        for(int i=N-2; i>=0; i--){
            if(i%B != B-1){
                suf[i] = Calc(suf[i], suf[i+1]);
            }
        }
    }
    T ask_O1(int l, int r){
        if(l > r) return IE;
        assert(0<=l && r < N);
        if(l/B != r/B){
            T ans = Calc(suf[l], pre[r]);
            l = l/B + 1;
            r = r/B - 1;
            if(l <= r){
                int k = std::__lg(r-l+1);
                ans = Calc(ans, a[k][l], a[k][r + 1 - (1<<k)]);
            }
            return ans;
        }else{
            int k = std::__lg(r-l+1);
            return Calc(b[k][l], b[k][r + 1 - (1<<k)]);
        }
    }
    T ask(int l, int r){
        if(l > r) return IE;
        assert(0<=l && r < N);
        T ans = IE;
        if(l/B != r/B){
            ans = suf[l];
            T ansr = pre[r];
            l = l/B + 1; r = r/B - 1;
            while(l <= r){
                int k = std::__lg(r-l+1);
                ans = Calc(ans, a[k][l]);
                l += (1<<k);
            }
            ans = Calc(ans, ansr);
        }else{
            while(l <= r){
                int k = std::__lg(r-l+1);
                ans = Calc(ans, b[k][l]);
                l += (1<<k);
            }
        }
        return ans;
    }
};

#endif