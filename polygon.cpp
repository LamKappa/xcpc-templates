#ifndef POLYGON
#define POLYGON
// 平面凸包

namespace Polygon{
    constexpr double PI = acosl(-1.);
    constexpr double INF = 1e20;
    constexpr double EPS = 1e-12;
    struct Point{
        double x,y;
        Point(double _x=0, double _y=0):x(_x),y(_y){}
        Point& operator=(const Point&o){
            x = o.x; y = o.y;
            return *this;
        }
        friend std::ostream& operator<<(std::ostream&out, const Point&p){
            return out<<'('<<p.x<<','<<p.y<<')';
        }
        bool operator<(const Point&o)const{
            return x==o.x?y<o.y:x<o.x;
        }
        Point operator+(const Point&o)const{
            return {x+o.x,y+o.y};
        }
        Point operator-()const{
            return {-x,-y};
        }
        Point operator-(const Point&o)const{
            return (*this) + (-o);
        }
        double operator*(const Point&o)const{
            return dot(o);
        }
        bool operator==(const Point&o)const{
            return between(o,o);
        }
        double norm()const{
            return sqrt((*this) * (*this));
        }
        double dot(const Point&o)const{
            return x*o.x+y*o.y;
        }
        double cross(const Point&o)const{
            return x*o.y-y*o.x;
        }
        bool between(Point a,Point b)const{
            if(abs((a.x-x)*(b.y-y)-(b.x-x)*(a.y-y))>EPS)return false;
            if(a.x>b.x)swap(a.x,b.x);
            if(a.y>b.y)swap(a.y,b.y);
            return a.x-EPS<=x&&x<=b.x+EPS &&
                a.y-EPS<=y&&y<=b.y+EPS;
        }
        double theta()const{
            return x+y?atan2(y,x):-INF;
        }
    };
    using Points = std::vector<Point>;
    bool check_inside(Point&p, const Points&dots){
        static std::mt19937 eng(std::random_device{}());
        static std::uniform_real_distribution<> dis(0,PI/2.);
        double K = dis(eng);
        double B = p.y - K*p.x;

        int cnt = 0;
        for(int i=0;i<dots.size();i++){
            int j = i+1==dots.size()?0:i+1;
            Point pi{dots[i]}, pj{dots[j]}, pk;

            double k = (pi.y-pj.y) / (pi.x-pj.x);
            double b = pi.y - k*pi.x;

            pk.x = (b-B) / (K-k);
            if(!finite(k)) pk.x = pi.x;
            pk.y = K*pk.x + B;

            if(pk.x>=p.x && pk.between(pi,pj))cnt++;
        }
        return cnt&1;
    }
    namespace Graham{
        Points chull(const Points&_dots){
            Points dots{_dots};
            Point c = *min_element(dots.begin(), dots.end());
            sort(dots.begin(), dots.end(), [&c](auto a,auto b){
                auto ac = a - c, bc = b - c;
                double delta = ac.theta() - bc.theta();
                return delta?delta<0:ac.norm()<bc.norm();
            });
            Points H{dots[0],dots[1]};
            for(int i=2;i<dots.size();i++){
                while(H.size()>1&&(H.back()-H[H.size()-2]).cross(dots[i]-H.back())<=0){
                    H.pop_back();
                }
                H.push_back(dots[i]);
            }
            return H;
        }
    }
    namespace Andrew{
        Points chull(const Points&_dots){
            Points dots{_dots};
            sort(dots.begin(),dots.end());
            std::vector<int> I{0};
            std::vector<bool> used(dots.size(), true);
            for(int i=1;i<dots.size();i++){
                while(I.size()>1&&(dots[I.back()]-dots[I[I.size()-2]]).cross(dots[i]-dots[I.back()])<=0){
                    used[I.back()] = false; I.pop_back();
                }
                I.push_back(i);
            }
            used[0] = false;
            std::size_t s = I.size();
            for(int i=dots.size()-2;i>=0;i--){
                if(used[i]) continue;
                while(I.size()>s&&(dots[I.back()]-dots[I[I.size()-2]]).cross(dots[i]-dots[I.back()])<=0){
                    I.pop_back();
                }
                I.push_back(i);
            }
            Points H(I.size()-1);
            for(int i=0;i<I.size()-1;i++) H[i] = dots[I[i]];
            return H;
        }
    }
}

#endif