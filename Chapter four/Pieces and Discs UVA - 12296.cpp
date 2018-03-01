
//重点是在求圆与直线的交点，然后判断交点是否在线段上这一段，浮点误差

#include <bits/stdc++.h>
#define mem(ar,num) memset(ar,num,sizeof(ar))
#define me(ar) memset(ar,0,sizeof(ar))
#define lowbit(x) (x&(-x))
#define forn(i,n) for(int i = 0;i < n; ++i)
using namespace std;
typedef long long LL;
typedef unsigned long long ULL;
const int    prime = 999983;
const int    INF = 0x7FFFFFFF;
const LL     INFF =0x7FFFFFFFFFFFFFFF;
const double pi = acos(-1.0);
const double inf = 1e18;
const double eps = 1e-12;
const LL     mod = 1e9 + 7;
struct Point
{
    double x,y;
    Point(double x = 0,double y = 0):x(x),y(y) {}

};
typedef Point Vector;
typedef vector<Point> Polygon;
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
    if( fabs(x) < eps)
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
double Cross(Vector A,Vector B)
{
    return A.x*B.y - A.y*B.x;
}
double angle(Vector v)//求向量的角度从0到2*pi
{
    return atan2(v.y,v.x);
}
double Angle(Vector A,Vector B)
{
    return acos(Dot(A,B)/Length(A)/Length(B));
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
    return dcmp(Cross(a1-p,a2-p)) == 0&&dcmp(Dot(a1-p,a2-p)) < 0;
}
int isPointInPolygon(Point p,Polygon poly)
{
    int n  = poly.size();
    int wn = 0;
    for(int i = 0; i < n; ++i)
    {
        if(Onsegment(p,poly[i],poly[(i+1)%n]))
            return -1;
        int k = dcmp(Cross(poly[(i+1)%n]-poly[i],p-poly[i]));
        int d1 = dcmp(poly[i].y-p.y);
        int d2 = dcmp(poly[(i+1)%n].y-p.y);
        if(k>0 && d1 <= 0&&d2 > 0)
            wn ++;
        if(k<0 && d2 <= 0&&d1 > 0)
            wn --;
    }
    if(wn != 0)
        return 1;
    return 0;
}
double PolygonArea(Polygon poly)
{
    double area = 0;
    int n = poly.size();
    for(int i = 1; i < n-1; i++)
        area += Cross(poly[i]-poly[0], poly[(i+1)%n]-poly[0]);
    return area/2;
}

Polygon CutPolygon(Polygon poly,Point A,Point B)
{
    Polygon newpoly;
    int n = poly.size();
    for(int i = 0; i < n; ++i)
    {
        Point C = poly[i];
        Point D = poly[(i+1)%n];
        if(dcmp(Cross(B-A,C-A)) >= 0)
            newpoly.push_back(C);
        if(dcmp((Cross(B-A,D-C))) != 0)
        {
            Point ip = Get_Line_Intersection(A,B-A,C,D-C);
            if(Onsegment(ip,C,D))
                newpoly.push_back(ip);
        }
    }

    return newpoly;
}

Polygon Init(int L,int W)//初始的矩形
{
    Polygon poly;
    poly.push_back(Point(0,0));
    poly.push_back(Point(L,0));
    poly.push_back(Point(L,W));
    poly.push_back(Point(0,W));
    return poly;
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
struct Line
{
    Line() = default;
    Vector v;
    Point p;
    Line(Vector V,Point P):v(V),p(P) {}
    Point point(double t)
    {
        return p + v*t;
    }
};
int Get_Line_Circle_Intersection(Line L,Circle C,vector<Point> &sol,double &t1,double &t2)
{
    double a = L.v.x,b = L.p.x-C.c.x,c = L.v.y,d = L.p.y - C.c.y;
    double e = a*a + c*c,f = 2*(a*b+c*d),g = b*b+d*d-C.r*C.r;
    double delta = f*f-4*e*g;
    if(dcmp(delta)<0)
        return 0;//相离
    if(dcmp(delta)==0)//相切
    {
        t1 = t2 = -f/(2*e);
        sol.push_back(L.point(t1));
        return 1;
    }
    t1 = (-f+sqrt(delta))/(2*e);
    t2 = (-f-sqrt(delta))/(2*e);
    sol.push_back(L.point(t1));
    sol.push_back(L.point(t2));
    return 2;
}
bool Intersect(int x,int y,int r,Polygon poly)
{
    int n = poly.size();
    Point c(x,y);
    if(isPointInPolygon(c,poly) != 0)
        return true;
    for(int i =  0; i < n; ++i)
    {
        if(dcmp(Length(poly[i] - c) - r) < 0)
            return true;
        Point mid = (poly[i]+poly[(i+1)%n])/2;
        if(dcmp(Length(mid     - c) - r) < 0)
            return true;
        vector<Point> sol;
        double t1,t2;
        int tmp = Get_Line_Circle_Intersection(Line(poly[(i+1)%n]-poly[i],poly[i]),Circle(c,r),sol,t1,t2);
        if(tmp==2)
        {
           /* for(int i = 0; i < sol.size(); ++i)
            {
                if(Onsegment(sol[i],poly[i],poly[(i+1) % n]))
                    return true;
            }*/
           if(dcmp(t1) > 0 && dcmp(t1-1) < 0)
              return true; // 端点在圆上
           if(dcmp(t2) > 0 && dcmp(t2-1) < 0)
               return true;
        }
    }
    return false;
}
int main(void)
{
    int n,m,L,W;
    while(cin>>n>>m>>L>>W&&L)
    {
        vector<Polygon> ans;
        ans.push_back(Init(L,W));
        while(n--)
        {
            vector<Polygon> tmp;
            int len = ans.size();
            int x1,y1,x2,y2;
            scanf("%d%d%d%d",&x1,&y1,&x2,&y2);
            for(int i  = 0; i < len; ++i)
            {
                Polygon p1  = CutPolygon(ans[i],Point(x1,y1),Point(x2,y2));
                Polygon p2  = CutPolygon(ans[i],Point(x2,y2),Point(x1,y1));
                if(dcmp(PolygonArea(p1)) > 0)
                    tmp.push_back(p1);
                if(dcmp(PolygonArea(p2)) > 0)
                    tmp.push_back(p2);
            }
            ans = tmp;
        }
        int len = ans.size();
        vector<double> Area(len+1);
        for(int i = 0; i < len; ++i)
        {
            Area[i] = fabs(PolygonArea(ans[i]));
        }
        while(m--)
        {
            vector<double> vec;
            int x,y,r;
            scanf("%d%d%d",&x,&y,&r);
            for(int i = 0; i < len; ++i)
                if(Intersect(x,y,r,ans[i]))
                    vec.push_back(Area[i]);
            sort(vec.begin(),vec.end());
            int t = vec.size();
            printf("%d",t);
            for(int i = 0; i < t; ++i)
                printf(" %.2f",vec[i]);

            puts("");
        }
        puts("");
    }
    return 0;
}
