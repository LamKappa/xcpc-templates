#ifndef SIEVE_CPP
#define SIEVE_CPP
// 线性筛

namespace Seive{
    std::vector<int> minp,primes,phi;
    void init(int N){
        minp.assign(N+1,0);
        phi.assign(N+1,0);
        primes.clear();
        
        for(int n=2;n<=N;n++){
            if(0==minp[n]){
                minp[n] = n;
                phi[n] = n-1;
                primes.push_back(n);
            }
            for(int p:primes){
                if(n*p>N)break;
                minp[n*p] = p;
                if(p==minp[n]){
                    phi[n*p] = phi[n] * p;
                    break;
                }
                phi[n*p] = phi[n] * (p-1);
            }
        }
    }
}


#endif