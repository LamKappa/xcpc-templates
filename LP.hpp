#ifndef LINEAR_PROGRAMMING_HPP
#define LINEAR_PROGRAMMING_HPP
// 线性规划

enum Relation{
    LESS_THAN,
    EQUAL_TO,
    GREATER_THAN
};

struct LP{
    constexpr static double EPS = 1e-12;
    std::size_t n = 0, m = 0;
    std::vector<int> base;
    std::vector<double> b, c;
    std::vector<std::vector<double>> A;

    // stdandard: Ax <= b maximize cx
    // assert xi >= 0
    // initialize: Ajx + x{n+j} == bj

    void add_restrict(auto&&aj, Relation r, double bj){
        n = std::max(n, aj.size());
        if(r==LESS_THAN || r==EQUAL_TO){
            A.push_back(aj); b.push_back(bj); m++;
        }
        if(r==GREATER_THAN || r== EQUAL_TO){
            for(auto&aji : aj) aji = -aji;
            A.push_back(aj); b.push_back(-bj); m++;
        }
    }

    void pivot(int x, int j){
        double Ajx = A[j][x], cx = c[x];
        base[j] = x; b[j] /= Ajx;
        for(int i=0; i<n+m+1; i++){
            A[j][i] /= Ajx;
            c[i] -= cx * A[j][i];
        }
        c[n + m + 1] += cx * b[j];
        for(int jj=0; jj<m; jj++){
            if(jj == j) continue;
            double C = A[jj][x];
            for(int i=0; i<n+m+1; i++){
                A[jj][i] -= C * A[j][i];
            }
            b[jj] -= C * b[j];
        }
    }

    void simplex(){
        bool fl = false;
        for(int i=0; i<n+m+1; i++){
            if(c[i] <= 0) continue;
            int min_j = -1;
            for(int j=0; j<m; j++){
                if(A[j][i] <= 0) continue;
                if(min_j < 0 ||
                    b[j]*A[min_j][i] < b[min_j]*A[j][i]){
                    min_j = j;
                }
            }
            assert(min_j >= 0);
            pivot(i, min_j);
            fl = true;
            break;
        }
        if(fl) return (void)simplex();
    }

    void initialize(){
        base.resize(m);
        double min_j = 0;
        for(int j=0; j<m; j++){
            A[j].resize(n+m+1, 0.);
            base[j] = n + j;
            A[j][base[j]] = 1.;

            if(b[j] < b[min_j]){
                min_j = j;
            }
        }
        if(b[min_j] < 0){
            for(int j=0; j<m; j++){
                A[j][n + m] = -1.;
            }
            auto _c = c;
            c.assign(n+m+2, 0.); c[n + m] = -1.;
            pivot(n + m, min_j);
            simplex();
            std::swap(c, _c);
            for(int j=0; j<m; j++){
                A[j][n + m] = 0.;
                if(std::abs(c[base[j]]) > EPS){
                    double cx = c[base[j]];
                    for(int i=0; i<n+m; i++){
                        c[i] -= cx * A[j][i];
                    }
                    c[n + m + 1] += cx * b[j];
                }
            }
        }
    }

    double maximize(auto&&_c){
        c = _c; c.resize(n+m+2, 0.);
        c[n + m + 1] = c[n]; c[n] = 0.;
        initialize();
        simplex();
        double res = c[n + m + 1];
        for(int j=0; j<m; j++){
            res += b[j] * c[base[j]];
        }
        return res;
    }
    
    double minimize(auto&&_c){
        c = _c; for(auto&x : c) x = -x;
        return -maximize(c);
    }
};

// e.g.

// LP lp;
// lp.add_restrict(vector{1.,1.,3.}, LESS_THAN, 30.);
// lp.add_restrict(vector{2.,2.,5.}, LESS_THAN, 24.);
// lp.add_restrict(vector{4.,1.,2.}, LESS_THAN, 36.);
// cout<<lp.maximize(vector{3.,1.,2.})<<endl;   // 28

// LP lp;
// lp.add_restrict(vector{2.,-1.}, LESS_THAN, 2.);
// lp.add_restrict(vector{1.,-5.}, LESS_THAN, -4.);
// cout<<lp.maximize(vector{2.,-1.})<<endl;     // 2

#endif