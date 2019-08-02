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
const LL     INFF = 0x7FFFFFFFFFFFFFFF;
const double pi = acos(-1.0);
const double inf = 1e18;
const double eps = 1e-10;
const LL     mod = 1e9 + 7;
struct Point
{
    double x, y;

    Point(double x = 0, double y = 0): x(x), y(y) {}

};
typedef Point Vector;
Vector operator + (Vector A, Vector B)
{
    return Vector(A.x + B.x, A.y + B.y);
}
Vector operator - (Vector A, Vector B)
{
    return Vector(A.x - B.x, A.y - B.y);
}
Vector operator / (Vector A, double p)
{
    return Vector(A.x / p, A.y / p);
}
Vector operator * (Vector A, double p)
{
    return Vector(A.x * p, A.y * p);
}
double angle(Vector v)//求向量的角度从0到2*pi
{
    return atan2(v.y, v.x);
}
int dcmp(double x)
{
    if (fabs(x) < eps)
        return 0;
    else
        return x < 0 ? -1 : 1;
}
bool operator < (const Point &a, const Point &b)
{
    if (dcmp(a.x - b.x) == 0)
        return a.y < b.y;
    else
        return a.x < b.x;
}


bool operator == (const Point &a, const Point &b)
{
    return !dcmp(a.x - b.x) && !dcmp(a.y - b.y);
}
double Dot(Vector A, Vector B)
{
    return A.x * B.x + A.y * B.y;
}
double Length(Vector A)
{
    return sqrt(A.x * A.x + A.y * A.y);
}
double Angle(Vector A, Vector B)
{
    return acos(Dot(A, B) / Length(A) / Length(B));
}
double Cross(Vector A, Vector B)
{
    return A.x * B.y - A.y * B.x;
}
double Area2(Point A, Point B, Point C)
{
    return Cross(B - A, C - A);
}
Vector Rotate(Vector A, double rad)
{
    return Vector (A.x * cos(rad) - A.y * sin(rad), A.x * sin(rad) + A.y * cos(rad));
}
Vector Normal(Vector A)//单位法线
{
    double L = Length(A);
    return Vector(-A.y / L, A.x / L);
}
//调用前确保直线有唯一交点，当且仅当Cross(v,w)非0
Point Get_Line_Intersection(Point P, Vector v, Point Q, Vector w)
{
    Vector u = P - Q;
    double t = Cross(w, u) / Cross(v, w);
    return P + v * t;
}
double Distance_To_Line(Point P, Point A, Point B) //点到直线的距离
{
    Vector v1 = B - A, v2 = P - A;
    return fabs(Cross(v1, v2) / Length(v1));
}
double Distance_To_Segment(Point P, Point A, Point B)
{
    if (A == B)
        return Length(P - A);
    Vector v1 = B - A, v2 = P - A, v3 = P - B;
    if (dcmp(Dot(v1, v2)) < 0)
        return Length(v1);
    else if (dcmp(Dot(v1, v3)) > 0)
        return Length(v3);
    else
        return fabs(Cross(v1, v2)) / Length(v1);
}
Point Get_Line_Projection(Point P, Point A, Point B) //求投影点
{
    Vector v = B - A;
    return A + v * (Dot(v, P - A) / Dot(v, v));
}
//线段相交判定 相交不在线段的端点
bool Segment_Proper_Intersection(Point a1, Point a2, Point b1, Point b2)
{
    double c1 =  Cross(a2 - a1, b1 - a1), c2 = Cross(a2 - a1, b2 - a1),
           c3 =  Cross(b2 - b1, a2 - b1), c4 = Cross(b2 - b1, a1 - b1);
    return dcmp(c1) * dcmp(c2) < 0 && dcmp(c3) * dcmp(c4) < 0;
}
//判断点是否在线段上(不包括端点）
bool Onsegment(Point p, Point a1, Point a2)
{
    return dcmp(Cross(a1 - p, a2 - p)) == 0 && dcmp(Dot(a1 - p, a2 - p)) < 0;
}
int ConvexHull(Point *p, int n , Point *ch)
{
    sort(p, p + n);
    int m = 0;
    for (int i = 0; i < n; ++i)
    {
        while (m > 1 && Cross(ch[m - 1] - ch[m - 2], p[i] - ch[m - 2]) < 0) m--;
        ch[m++] = p[i];

    }
    int k = m;
    for (int i = n - 2; i >= 0; --i)
    {
        while (m > k && Cross(ch[m - 1] - ch[m - 2], p[i] - ch[m - 2]) < 0) m--;
        ch[m++] = p[i];
    }
    if (n > 1) m--;
    return m;
}
const int maxn = 100 + 10;
Point P1[maxn], P2[maxn], p1[maxn], p2[maxn];
int x[maxn], y[maxn], v[maxn];
int isPointInPolygon(Point p, Point * poly, int n)
{
    // int n  = poly.size();
    int wn = 0;
    for (int i = 0; i < n; ++i)
    {
        if (Onsegment(p, poly[i], poly[(i + 1) % n])) return -1;
        int k = dcmp(Cross(poly[(i + 1) % n] - poly[i], p - poly[i]));
        int d1 = dcmp(poly[i].y - p.y);
        int d2 = dcmp(poly[(i + 1) % n].y - p.y);
        if (k > 0 && d1 <= 0 && d2 > 0) wn ++;
        if (k < 0 && d2 <= 0 && d1 > 0) wn --;
    }
    if (wn != 0)  return 1;
    return 0;
}
bool solve() {
    int n;
    cin >> n;
    int cnt1, cnt2;
    cnt1 = cnt2 = 0;
    for (int i = 1; i <= n; ++i) {
        scanf("%d%d%d", &x[i], &y[i], &v[i]);
        for (int j = 1; j < i; ++j)
        {
            if (x[i] == x[j] && y[i] == y[j]) {
                if (v[i] != v[j])
                    return false;
            }
            else
                continue;
        }
        if (v[i] == -1) {
            P1[cnt1].x = x[i];
            P1[cnt1++].y = y[i];

        }
        else {
            P2[cnt2].x = x[i];
            P2[cnt2++].y = y[i];
        }

    }
    if (cnt1 <= 1 && cnt2 <= 1) return true;
    int m1 = ConvexHull(P1, cnt1, p1);
    int m2 = ConvexHull(P2, cnt2, p2);
    if (m1 == 2) {
        for (int i = 0; i < m2; ++i)
            if (Onsegment(p2[i], p1[0], p1[1]))
                return false;
    }
    if (m2 == 2) {
        for (int i = 0; i < m1; ++i)
            if (Onsegment(p1[i], p2[0], p2[1]))
                return false;
    }
    if (m1 >= 3) {
        for (int i = 0; i < m2; ++i)
            if (isPointInPolygon(p2[i], p1, m1))
                return false;
    }
    if (m2 >= 3) {
        for (int i = 0; i < m1; ++i) {
            if (isPointInPolygon(p1[i], p2, m2))
                return false;
        }
    }
    return true;
}
int main(void) {
    int T;
    cin >> T;
    while (T--) {
        if (solve())
            puts("Successful!");
        else
            puts("Infinite loop!");


    }




    return 0;
}
