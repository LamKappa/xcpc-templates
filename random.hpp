#ifndef RANDOM
#define RANDOM
// 随机

namespace Random{
    std::mt19937 eng(std::random_device{}());
    double double_uniform(double l,double r){
        std::uniform_real_distribution<> dis(l,r);
        return dis(eng);
    }
    int int_uniform(int l,int r){
        std::uniform_int_distribution<> dis(l,r);
        return dis(eng);
    }
};

#endif