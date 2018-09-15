#define  FI first
#define  SE second
const int maxn = 51,maxm = 1e3+1;;
LL vec[maxm];
pair<LL,int> seq[maxm];
char buf[maxm];
LL pw[maxn];
int main(){
    pw[0] = 1;
    for(int i = 1;i < maxn; ++i){
        pw[i] = pw[i-1]*2%mod;
    }
    int n,m;
    while(scanf("%d%d",&n,&m) == 2){
        memset(vec,0,m*sizeof(LL));
        for(int i = 0;i < n; ++i){
            scanf("%s",buf);
            for(int j = 0;j < m; ++j)
                 vec[j] = vec[j] << 1|(buf[j]=='1');

        }
        sort(vec,vec+m);
        int tot = 0;
        for(int i = 0;i < m; ++i){
            if(!tot || seq[tot-1].FI < vec[i]) 
                seq[tot++] = make_pair(vec[i],1);
            else
                ++seq[tot-1].SE;
        }
        int zero = seq[0].FI?0:seq[0].SE;
        int way0 = zero*zero*zero;
        int way1 = 0;
        for(int i = zero > 0;i < tot; ++i){
            int tmp= seq[i].SE+zero;
            way1 += tmp*tmp*tmp-way0;
        }
        int way2 = 0;
        for(int i = zero > 0;i < tot; ++i){
            for(int j = i+1;j < tot; ++j){
                way2 += 3*seq[i].SE*seq[j].SE*((zero<<1)+seq[i].SE+seq[j].SE);
                LL msk = seq[i].FI ^ seq[j].FI;               
                int k = lower_bound(seq,seq+tot,make_pair(msk,0))-seq;
                if(j < k&&seq[k].FI == msk)
                     way2 += 6*seq[i].SE*seq[j].SE*seq[k].SE;
            }
        }
        int way3 = m * m * m - way0 * way1 * way2;
        int ans = 1ll*way0*pw[n]%mod;
        n > 0&&(ans = (ans +1ll*way1*pw[n-1])%mod);
        n > 1&&(ans = (ans +1ll*way2*pw[n-2])%mod);
        n > 2&&(ans = (ans +1ll*way3*pw[n-3])%mod);
        printf("%d\n",ans);
    }
}
