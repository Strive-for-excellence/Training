#include <cstdio>//C语言io
#include <cstring>//以下是c语言常用头文件
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <cctype>
#include <cstring>
#include <cmath>
#include <iostream>//c++IO
#include <sstream>
#include <string>
#include <list>//c++常用容器
#include <vector>
#include <set>
#include <map>
#include <queue>
#include <stack>
#include <algorithm>//c++泛型的一些函数
#include <functional>//用来提供一些模版
#define fo0(i,n) for(int i = 0;i < n; ++i)
#define fo1(i,n) for(int i = 1;i <= n; ++i)
#define mem(ar,num) memset(ar,num,sizeof(ar))
#define me(ar) memset(ar,0,sizeof(ar))
#define lowbit(x) (x&(-x))
using namespace std;
typedef long long LL;
typedef unsigned long long ULL;
const int    prime = 999983;
const int    INF = 0x7FFFFFFF;
const LL     INFF =0x7FFFFFFFFFFFFFFF;
const double pi = acos(-1.0);
const double inf = 1e18;
const double eps = 5 * 1e-13;
const LL     mod = 1e9 + 7;
//......................................................

struct Point
{
//    Point() = default;
    double x,y;
    Point(double x = 0,double y = 0):x(x),y(y) {}

};
typedef Point Vector;
Vector operator + (Vector A,Vector B)
{
    return Vector(A.x + B.x,A.y + B.y);
}
Vector operator - (Vector A,Vector B)
{
    return Vector(A.x-B.x,A.y-B.y);
}
Vector operator / (Vector A,double p)
{
    return Vector(A.x/p,A.y/p);
}
Vector operator * (Vector A,double p)
{
    return Vector(A.x*p,A.y*p);
}
int dcmp(double x)
{
    if(fabs(x)<eps)
        return 0;
    else
        return x < 0?-1:1;
}
bool operator < (const Point &a,const Point &b)
{
    if(dcmp(a.x-b.x)==0)
        return a.y<b.y;
    else
        return a.x<b.x;
}

bool operator == (const Point &a,const Point &b)
{
    return !dcmp(a.x-b.x)&&!dcmp(a.y-b.y);
}
double Dot(Vector A,Vector B)
{
    return A.x*B.x+A.y*B.y;
}
double Length(Vector A)
{
    return sqrt(A.x*A.x+A.y*A.y);
}
double Angle(Vector A,Vector B)
{
    return acos(Dot(A,B)/Length(A)/Length(B));
}
double Cross(Vector A,Vector B)
{
    return A.x*B.y - A.y*B.x;
}
double Area2(Point A,Point B,Point C)
{
    return Cross(B-A,C-A);
}
Vector Rotate(Vector A,double rad)
{
    return Vector (A.x*cos(rad)-A.y*sin(rad),A.x*sin(rad)+A.y*cos(rad));
}
Vector Normal(Vector A)//单位法线
{
    double L = Length(A);
    return Vector(-A.y/L,A.x/L);
}
//调用前确保直线有唯一交点，当且仅当Cross(v,w)非0
Point Get_Line_Intersection(Point P,Vector v,Point Q,Vector w)
{
    Vector u = P - Q;
    double t = Cross(w,u)/Cross(v,w);
    return P+v*t;
}
double angle(Vector v)
{
    return atan2(v.y,v.x);
}

double Distance_To_Line(Point P,Point A,Point B)//点到直线的距离
{
    Vector v1 = B-A,v2 = P-A;
    return fabs(Cross(v1,v2)/Length(v1));
}
double Distance_To_Segment(Point P,Point A,Point B)
{
    if(A==B)
        return Length(P-A);
    Vector v1 = B-A,v2 = P-A,v3 = P-B;
    if(dcmp(Dot(v1,v2))<0)
        return Length(v1);
    else if(dcmp(Dot(v1,v3))>0)
        return Length(v3);
    else
        return fabs(Cross(v1,v2))/Length(v1);
}
Point Get_Line_Projection(Point P,Point A,Point B)//求投影点
{
    Vector v = B- A;
    return A + v*(Dot(v,P-A)/Dot(v,v));
}
//线段相交判定 相交不在线段的端点
bool Segment_Proper_Intersection(Point a1,Point a2,Point b1,Point b2)
{
    double c1 =  Cross(a2-a1,b1-a1),c2 = Cross(a2-a1,b2-a1),
           c3 =  Cross(b2-b1,a2-b1),c4 = Cross(b2-b1,a1-b1);
    return dcmp(c1)*dcmp(c2)<0&&dcmp(c3)*dcmp(c4)<0;
}
//判断点是否在线段上(不包括端点）
bool Onsegment(Point p,Point a1,Point a2)
{
    return dcmp(Cross(a1-p,a2-p))==0&&dcmp(Dot(a1-p,a2-p))<0;
}
//多边形的有向面积
//1
double PolygonArea (Point * p,int n)
{
    double area = 0;
    for(int i = 1; i < n - 1; ++i)
    {
        area += Cross(p[i]-p[0],p[i+1]-p[0]);
    }
    return area/2;
}
//2
double PolygonArea2(Point *p,int n)
{
    double area = 0;
    for(int i = 0; i <= n-1; ++i)
    {
        if(i!=n-1)
            area += Cross(p[i],p[i+1]);
        else
            area += Cross(p[n-1],p[0]);
    }
    return area/2;
}
Point Get(Point A,Point B,double rad1,double rad2)
{
    return Get_Line_Intersection(A,Rotate(B-A,rad1/3),B,Rotate(A-B,-rad2/3));
}
struct Circle
{
    Circle() =  default;
    Point c;
    double r;
    Circle(Point c,double r):c(c),r(r) {}
    Point point(double a)
    {
        return Point (c.x+cos(a)*r,c.y+sin(a)*r);
    }
};

Point read_point(void)
{
    Point p;
    scanf("%lf%lf",&p.x,&p.y);
    return p;
}
void get_Circle_Cricle_Intersection(Circle C1,Circle C2,vector<double> &sol)
{
    double d =  Length(C1.c-C2.c);
    if(dcmp(d)==0)
        return ;
    if(dcmp(C1.r+C2.r-d)<0)
        return ;
    if(dcmp(fabs(C1.r-C2.r)-d)>0)
        return ;
    double a = angle(C2.c-C1.c);
    double da = acos((C1.r*C1.r+d*d-C2.r*C2.r)/(2*C1.r*d));
    sol.push_back(a-da);
    sol.push_back(a+da);
}
const int maxn = 100+12;
Circle cir[maxn];
int n;
bool vis[maxn];
int judge(Point p)
{
    for(int i = n-1; i >= 0; --i)
        if(Length(p-cir[i].c)<cir[i].r)
            return i;
    return -1;
}
int main(void)
{
    while(cin>>n&&n)
    {
        me(vis);
        for(int i = 0; i < n; ++i)
            cin>>cir[i].c.x>>cir[i].c.y>>cir[i].r;
        for(int i = 0; i < n ; ++i)
        {
            vector<double> rad;
            rad.push_back(0);
            rad.push_back(2*pi);
            for(int j = 0; j < n; ++j)
                get_Circle_Cricle_Intersection(cir[i],cir[j],rad);
            sort(rad.begin(),rad.end());
            for(int j = 0; j < rad.size()-1; ++ j)
            {

                double    mid = (rad[j] + rad[j+1]) / 2.0;

                for(int k = -1; k <= 1; k += 2)
                {
                    double rr =  cir[i].r + k*eps;
                    int t = judge(Point(cir[i].c.x + rr*cos(mid),cir[i].c.y+rr*sin(mid)));
                    if(t >= 0)
                        vis[t] = true;
                }
            }
        }
        int num =  0;
        for(int i = 0; i < n; ++i)
            if(vis[i])
                ++num;
//        cout<<endl;
        cout<<num<<endl;
//        cout<<endl;

    }

    return 0;
}
