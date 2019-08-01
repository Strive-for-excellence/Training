#include <bits/stdc++.h>
#define mem(ar,num) memset(ar,num,sizeof(ar))
#define me(ar) memset(ar,0,sizeof(ar))
#define lowbit(x) (x&(-x))
#define Pb push_back
#define  FI first
#define  SE second
#define rep(i,a,n) for (int i=a;i<n;i++)
#define per(i,a,n) for (int i=n-1;i>=a;i--)
#define IOS ios::sync_with_stdio(false)
#define DEBUG cout<<endl<<"DEBUG"<<endl;
using namespace std;
typedef long long LL;
typedef unsigned long long ULL;
const int    prime = 999983;
const LL    INF = 1e14;
const LL     INFF = 0x7FFFFFFFFFFFFFFF;
const double pi = acos(-1.0);
const double inf = 1e18;
const double eps = 1e-6;
const LL     mod = 1e9 + 7;
LL qpow(LL a, LL b) {LL s = 1; while (b > 0) {if (b & 1)s = s * a % mod; a = a * a % mod; b >>= 1;} return s;}
LL gcd(LL a, LL b) {return b ? gcd(b, a % b) : a;}
int dr[2][4] = {1, -1, 0, 0, 0, 0, -1, 1};
typedef pair<LL, int> P;
const int maxn = 10000 + 10;
typedef pair<int, int> PII;
struct Edge2 {
    int to;
    LL w;
};
vector<Edge2> G[maxn];
int n, m;
LL dis[maxn];
bool done[maxn];
void Dij() {
    for (int i = 1; i <= n; ++i)
        dis[i] = INF, done[i] = 0;
    priority_queue<P, vector<P>, greater<P>> Q;
    dis[1] = 0;
    Q.push(P(dis[1], 1));

    while (!Q.empty()) {
        P p = Q.top(); Q.pop();
        // cout << u << endl;
        int u = p.second;
        if (done[u]) continue;
        done[u] = true;
        for (int i = 0; i < G[u].size(); ++i) {
            int v = G[u][i].to;
            if (!done[v] && dis[v] > dis[u] + G[u][i].w) {
                dis[v] = dis[u] + G[u][i].w;
                Q.push(P(dis[v], v));
            }
        }
    }
}

const int LEN = maxn;

struct Edge {
    int from, to;
    LL  cap, flow;
    Edge(int u, int v, int w, int f): from(u), to(v), cap(w), flow(f) {}
};
struct Dinic {
    LL n, m, s, t;
    vector<Edge> edges;
    vector<LL> G[LEN];
    LL a[LEN];
    LL vis[LEN];
    LL d[LEN];
    LL cur[LEN];//好吧就是点, 代表该点在一次求增广的过程中搜索到了那条边,意思就是从这条边往下肯定搜索不到结果了
    void init(LL n)
    {
        this->n  = n;
        for (LL i = 0; i < n; ++i)
            G[i].clear();
        edges.clear();
    }
    void Add(LL u, LL v, LL w)
    {
        edges.push_back(Edge(u, v, w, 0));
        edges.push_back(Edge(v, u, 0, 0));
        m = edges.size();
        G[u].push_back(m - 2);
        G[v].push_back(m - 1);
    }
    bool Bfs(void)//分层
    {
        me(d);
        me(vis);
        d[s] = 0;
        vis[s] = 1;

        queue<LL> Q;
        Q.push(s);
        while (!Q.empty())
        {
            LL q = Q.front(); Q.pop();

            for (size_t i = 0; i < G[q].size(); ++i)
            {
                Edge &tmp = edges[G[q][i]];
                if (!vis[tmp.to] && tmp.cap > tmp.flow)
                {
                    vis[tmp.to] = 1;
                    d[tmp.to] = d[q] + 1;
                    Q.push(tmp.to);
                }
            }
        }
        return vis[t];
    }
    LL Dfs(LL node, LL a)
    {

        if (node == t || a == 0)
            return a;
        LL flow =  0, f;
        for (LL &i = cur[node]; i < G[node].size(); ++i)
        {
            Edge &tmp = edges[G[node][i]];
            if (d[tmp.to] == d[node] + 1 && (f = Dfs(tmp.to, min(a, tmp.cap - tmp.flow))) > 0)
            {
                flow += f;
                tmp.flow += f;
                edges[G[node][i] ^ 1].flow -= f;
                a -= f;
                if (a == 0)
                    break;
            }
        }
        return flow;
    }
    LL MaxFlow(LL s, LL t)
    {
        this->s = s;
        this->t = t;
        LL flow = 0;
        while (Bfs())
        {
            me(cur);
            flow += Dfs(s, maxn);
        }
        return flow;

    }


};
Dinic dinic;
int  main(void)
{

    LL T;
    cin >> T;
    while (T--) {

        cin >> n >> m;
        for (int i = 1; i <=  n; ++i) G[i].clear();
        for (int i = 1; i <=  m; ++i) {
            int u, v, w;
            scanf("%d%d%d", &u, &v, &w);
            G[u].Pb(Edge2{v, w});
        }
        Dij();
        if (dis[n] == INF) {
            cout << 0 << endl;
            continue;
        }
        dinic.init(n + 2);
        for (int u = 1; u <= n; ++u) {
            for (int j = 0; j < G[u].size(); ++j) {
                auto e = G[u][j];
                if (dis[e.to]  == dis[u] + e.w) {
                    dinic.Add(u, e.to, e.w);
                    // cout << u << " " << e.to << " " << e.w << endl;
                }
            }
        }
        cout << dinic.MaxFlow(1, n) << endl;

    }


    return 0;
}
