const int maxn = 1e5+100;
LL dp[maxn];
int main(void)
{
    LL N,M;
    dp[0] = 1;
    dp[1] = 0;
    dp[2] = 1;
    while(scanf("%lld %lld",&N,&M) == 2){
        if(M == 1) puts("1");
        rep(n,3,N+1)
        dp[n] = 1ll*(n-1)*dp[n-1]%M+1ll*(n-1)*dp[n-2]%M-1ll*(n-1)*(n-2)/2%M*dp[n-3]%M;
        printf("%lld\n",(dp[N]%M+M)%M);
    }


   return 0;
}
