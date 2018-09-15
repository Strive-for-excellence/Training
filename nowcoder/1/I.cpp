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
const int    INF = 0x7FFFFFFF;
const LL     INFF =0x7FFFFFFFFFFFFFFF;
const double pi = acos(-1.0);
const double inf = 1e18;
const double eps = 1e-6;
const LL     mod = 1e9 + 7;
LL qpow(LL a,LL b){LL s=1;while(b>0){if(b&1)s=s*a%mod;a=a*a%mod;b>>=1;}return s;}
LL gcd(LL a,LL b) {return b?gcd(b,a%b):a;}
int dr[2][4] = {1,-1,0,0,0,0,-1,1};
typedef pair<int,int> P;
const int maxn = 1e6 + 10;
 
struct SuffixArray {
  int s[maxn];      // 原始字符数组（最后一个字符应必须是0，而前面的字符必须非0）
  int sa[maxn];     // 后缀数组
  int rank[maxn];   // 名次数组. rank[0]一定是n-1，即最后一个字符
  int height[maxn]; // height数组
  int t[maxn], t2[maxn], c[maxn]; // 辅助数组
  int n; // 字符个数
 
  void clear() { n = 0; memset(sa, 0, sizeof(sa)); }
 
  // m为最大字符值加1。调用之前需设置好s和n
  void build_sa(int m) {
    int i, *x = t, *y = t2;
    for(i = 0; i < m; i++) c[i] = 0;
    for(i = 0; i < n; i++) c[x[i] = s[i]]++;
    for(i = 1; i < m; i++) c[i] += c[i-1];
    for(i = n-1; i >= 0; i--) sa[--c[x[i]]] = i;
    for(int k = 1; k <= n; k <<= 1) {
      int p = 0;
      for(i = n-k; i < n; i++) y[p++] = i;
      for(i = 0; i < n; i++) if(sa[i] >= k) y[p++] = sa[i]-k;
      for(i = 0; i < m; i++) c[i] = 0;
      for(i = 0; i < n; i++) c[x[y[i]]]++;
      for(i = 0; i < m; i++) c[i] += c[i-1];
      for(i = n-1; i >= 0; i--) sa[--c[x[y[i]]]] = y[i];
      swap(x, y);
      p = 1; x[sa[0]] = 0;
      for(i = 1; i < n; i++)
        x[sa[i]] = y[sa[i-1]]==y[sa[i]] && y[sa[i-1]+k]==y[sa[i]+k] ? p-1 : p++;
      if(p >= n) break;
      m = p;
    }
  }
 
  void build_height() {
    int i, j, k = 0;
    for(i = 0; i < n; i++) rank[sa[i]] = i;
    for(i = 0; i < n; i++) {
      if(k) k--;
      int j = sa[rank[i]-1];
      while(s[i+k] == s[j+k]) k++;
      height[rank[i]] = k;
    }
  }
};
 
 
//本题相关
char ar[maxn];
SuffixArray  SA;
int n;
LL cal1(){
    LL ans = 0;
    rep(i,6,SA.n){
        LL Left = n-SA.sa[i]%(n+1);
        ans += max(Left-SA.height[i],0LL);
    }
    return ans;
}
LL cal2(){
   int maxlen = 1,len = 1;
   rep(i,0,n){
     if(ar[i]==ar[i-1])
        len++;
     else
        len = 1;
     maxlen = max(len,maxlen);
   }   
   return maxlen;
}
int main(void){
 
    while(cin>>n){
        SA.clear();
        scanf("%s",ar);
        // cout<<ar<<endl;
        // rep(i,0,n) ar[i;
        int a[4] = {1,2,3};
        do{
           rep(i,0,n) SA.s[SA.n++] = a[ar[i]-'a'];
           SA.s[SA.n++] = 0;
        } while(next_permutation(a,a+3));
        SA.build_sa(10);
        SA.build_height();
        if(1){
            LL a,b;
            a = cal1();
            b = cal2();
            LL ans = (a+3*b)/6;
            printf("%lld\n",ans);
        }
         
    }
    return 0;
}
