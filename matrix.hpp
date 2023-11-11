#ifndef MATRIX
#define MATRIX
// 矩阵

template<std::size_t N, std::size_t M, typename I=long long>
struct Matrix{
    I val[N][M];

    I* operator[](int i){
        return val[i];
    }
    const I* operator[](int i) const {
        return val[i];
    }
    Matrix<N,M,I> operator+(const Matrix<N,M,I>&o)const{
        Matrix<N,M,I> res;
        for(int i=0;i<N;i++){
            for(int j=0;j<M;j++){
                res[i][j] = val[i][j] + o[i][j];
            }
        }
        return res;
    }
    template<std::size_t K>
    Matrix<N,K,I> operator*(const Matrix<M,K,I>&o)const{
        Matrix<N,K,I> res;
        for(int i=0;i<N;i++){
            for(int k=0;k<K;k++){
                res[i][k] = 0;
                for(int j=0;j<M;j++){
                    res[i][k] += val[i][j] * o[j][k];
                }
            }
        }
        return res;
    }
    friend std::ostream& operator<<(std::ostream&out, const Matrix<N,M,I>&x){
        for(int i=0;i<N;i++){
            for(int j=0;j<M;j++){
                out<<x[i][j]<<" \n"[j+1==M];
            }
        }
        return out;
    }
};

#endif