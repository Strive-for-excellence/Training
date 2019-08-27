typedef pair<int,int> P;
const int maxn = (1<<21)+10;
int F[maxn];
int A[maxn],B[maxn];
int main(void)
{
    int T;cin>>T;
    while(T--){
    	me(B);
    	int N;scanf("%d",&N);
    	for(int i = 1;i <= N; ++i){
    		scanf("%d",&A[i]);
    		B[A[i]]++;
    	}
    	for(int i =0;i <(1<<20); ++i)
    		F[i] = B[i];
    	for(int i = 0;i < 20; ++i){
    		for(int mask = 0;mask < (1<<20); ++mask){
    			if(mask>>i&1)
    				F[mask] += F[mask^(1<<i)];
    		}
    	}
    	LL ans = 0;
    	for(int mask = 0;mask < (1<<20); ++mask){
    		ans += 1ll*B[mask]*F[(1<<20)-1-mask];
    	}
    	cout<<ans<<endl;
    	
    }

   return 0;
}
