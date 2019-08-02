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
// const int    prime = 999983;
const int    INF = 0x7FFFFFFF;
const LL     INFF =0x7FFFFFFFFFFFFFFF;
const double pi = acos(-1.0);
const double inf = 1e18;
const double eps = 1e-6;
const LL     mod = 998244353;
LL qpow(LL a,LL b){LL s=1;while(b>0){if(b&1)s=s*a%mod;a=a*a%mod;b>>=1;}return s;}
LL gcd(LL a,LL b) {return b?gcd(b,a%b):a;}
int dr[2][4] = {1,-1,0,0,0,0,-1,1};
typedef pair<int,int> P;
const int maxn = 1e7+1000;
int g[maxn];
int prime[maxn/10];
 int phi[maxn];
bool vis[maxn];
void init(){
    int tot = 0;
    int N = 1e7+10;
    g[1] = 1;
    phi[1] = 1;
    for(int i = 2;i <= N; ++i){
        if(!vis[i]){
            prime[tot++] = i;
            g[i] = 2*i-1;
             phi[i] = i-1;
        }
        for(int j = 0;j < tot; ++j){
            if(1ll*prime[j]*i > N) break;
            vis[prime[j]*i] = true;
            
            if(i%prime[j]==0){
                 phi[i*prime[j]] = phi[i]*prime[j]%mod;
                LL t = prime[j];
                LL ii = i;
                LL cnt = 1;
                while(ii %prime[j] == 0)
                    ii /= prime[j],t *= prime[j],++cnt;
                // if(prime[j]*i == 9){
                //     cout<<cnt<<" "<<( g[prime[j]*i] =g[ii]*((1ll*(cnt+1)*t-1ll*cnt*t/prime[j])%mod)%mod)<<endl;
                // }
                 g[prime[j]*i] =g[ii]*((1ll*(cnt+1)*t-1ll*cnt*t/prime[j])%mod)%mod;
                 break;
            }
            else
                {
                     phi[i*prime[j]] = phi[i]*(prime[j]-1);
                    g[prime[j]*i] = 1ll*g[i]*g[prime[j]]%mod;
                }
        }
    }
}
// template<
inline void Add(int &a,int b){
    a += b;
    if(a >= mod)
        a -= mod;
}
inline void Add(LL &a,LL b){
    a += b;
    if(a >= mod)
        a -= mod;
}
LL Get(int p,int n){
    int ans =0;
    for(int i = 1;i <= n; ++i){
        Add(ans,gcd(i,p));
    }
    return ans;
}
LL Get2(int p,int n){
    int ans = 0;
    for(int i = 1;i *i <= p; ++i){
        if(p%i) continue;
        Add(ans,n/i*phi[i]);
        if(p/i != i)
        Add(ans,n/(p/i)*phi[p/i]);
    }
    return ans;
}
template <class T>
void read(T &x) {
    static char ch;static bool neg;
    for(ch=neg=0;ch<'0' || '9'<ch;neg|=ch=='-',ch=getchar());
    for(x=0;'0'<=ch && ch<='9';(x*=10)+=ch-'0',ch=getchar());
    x=neg?-x:x;
}
int main(void)
{
    // cout<<maxn<<endl;
    init();
//    for(int i =  1;i <= 100; ++i){
//        for(int j =1;j <= i; ++j){
//            cout<<i<<" "<<j<<" "<<Get(i,j)<<" "<<Get2(i,j)<<endl; 
//            assert(Get(i,j) == Get2(i,j)); 
//        }
//    }
    int T;cin>>T;
    while(T--){
        __int128 n;
        read(n);
        __int128 i;
        LL  ans = 0;
        for(i = 1;1ll*(i+1)*(i+1)*(i+1)-1 <= n; ++i){
            // assert(i < maxn);

            Add(ans,(1ll*g[i]*(3*i+3)+i)%mod);
        }
        
        // cout<<i<<endl;
        Add(ans,(n-i*i*i)/i*g[i]%mod);
        // aasert()
        if(n >= i*i*i)
            Add(ans,i);
        LL t = (n-i*i*i)%i;
        Add(ans,Get2(i,t));
        cout<<ans<<endl;
    }
    

   return 0;
}
