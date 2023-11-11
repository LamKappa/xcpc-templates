#ifndef LINEAR_PROGRAMMING
#define LINEAR_PROGRAMMING
// 线性规划 支持无解、无穷大判定 常数大 缓存命中差

enum Relation{
    LESS_EQ,
    EQUAL,
    GREATER_EQ
};

struct LP{
    constexpr static double EPS = 1e-12;
    std::size_t n = 0, m = 0;
    std::vector<int> base, varmp;
    std::vector<double> b, c;
    std::vector<std::vector<double>> A;

    bool exist = true;

    // stdandard: Ax <= b maximize cx
    // assert xi >= 0
    // initialize: Ajx + x{n+j} == bj

    void add_restrict(std::vector<double> aj, Relation r, double bj){
        n = std::max(n, aj.size());
        if(r != GREATER_EQ){
            A.push_back(aj); b.push_back(bj); m++;
        }
        if(r != LESS_EQ){
            A.push_back(aj); b.push_back(-bj); m++;
            for(auto&x : A.back()) x = -x;
        }
    }

    void pivot(int x, int j){
        double Ajx = A[j][x], cx = c[varmp[x]];
        b[j] /= Ajx;
        c[n + m + 1] += cx * b[j];
        for(int i=0; i<n+1; i++){
            c[varmp[i]] -= cx * (A[j][i] /= Ajx);
        }
        c[base[j]] -= cx / Ajx;
        A[j][x] = 1. / Ajx;

        for(int jj=0; jj<m; jj++){
            if(jj == j) continue;
            double C = A[jj][x];
            for(int i=0; i<n+1; i++){
                A[jj][i] -= C * A[j][i];
            }
            b[jj] -= C * b[j];
            A[jj][x] = -C / Ajx;
        }
        std::swap(varmp[x], base[j]);
    }

    void simplex(){
        for(;;){
            int max_i = -1;
            for(int i=0; i<n+1; i++){
                if(c[varmp[i]] <= 0) continue;
                if(max_i < 0 || c[varmp[max_i]] < c[varmp[i]]){
                    max_i = i;
                }
            }
            if(max_i < 0) return;
            int min_j = -1;
            for(int j=0; j<m; j++){
                if(A[j][max_i] <= EPS) continue;
                if(min_j < 0 ||
                    b[j]*A[min_j][max_i] < b[min_j]*A[j][max_i]){
                    min_j = j;
                }
            }
            if(min_j < 0) return (void)(exist = false);
            pivot(max_i, min_j);
        }
    }

    void initialize(){
        base.resize(m);
        varmp.resize(n+1);
        std::iota(varmp.begin(),varmp.end(),0);
        varmp[n] = n + m;
        int min_j = 0;
        for(int j=0; j<m; j++){
            A[j].resize(n+1, 0.);
            base[j] = n + j;

            if(b[j] < b[min_j]){
                min_j = j;
            }
        }
        if(b[min_j] < 0){
            for(int j=0; j<m; j++){
                A[j][n] = -1.;
            }
            decltype(c) _c; std::swap(_c, c);
            c.assign(n+m+2, 0.); c[n + m] = -1.;
            pivot(n + m, min_j); simplex(); 
            if(std::abs(c[n + m] + 1) > EPS){
                return (void)(exist = false);
            }
            std::swap(c, _c);
            for(int i=0; i<n; i++){
                for(int j=0; j<m; j++){
                    c[varmp[i]] -= c[base[j]] * A[j][i];
                }
            }
            for(int j=0; j<m; j++){
                c[n + m + 1] += c[base[j]] * b[j];
                A[j][n] = c[base[j]] = 0.;
            }
        }
    }

    double maximize(std::vector<double> _c){
        c = _c; c.resize(n+m+2, 0.);
        c[n + m + 1] = c[n]; c[n] = 0.;
        initialize();
        if(!exist) return NAN;
        simplex();
        if(!exist) return 1. / 0.;
        double res = c[n + m + 1];
        for(int j=0; j<m; j++){
            res += b[j] * c[base[j]];
        }
        return res;
    }
    
    double minimize(std::vector<double> _c){
        c = _c; for(auto&x : c) x = -x;
        return -maximize(c);
    }
};

#endif

// e.g.

// LP lp1;
// lp1.add_restrict({1.,1.,3.}, LESS_EQ, 30.);
// lp1.add_restrict({2.,2.,5.}, LESS_EQ, 24.);
// lp1.add_restrict({4.,1.,2.}, LESS_EQ, 36.);
// cout<<lp1.maximize({3.,1.,2.})<<endl;    // 28

// LP lp2;
// lp2.add_restrict({2.,-1.}, LESS_EQ, 2.);
// lp2.add_restrict({1.,-5.}, LESS_EQ, -4.);
// cout<<lp2.maximize({2.,-1.})<<endl;      // 2

// LP lp3;
// lp3.add_restrict({1,-1}, LESS_EQ, -2);
// cout<<lp3.maximize({1})<<endl;              // inf

// LP lp4;
// lp4.add_restrict({1}, LESS_EQ, -1);
// lp4.add_restrict({1}, GREATER_EQ, 1);
// cout<<lp4.maximize({1})<<endl;              // nan
