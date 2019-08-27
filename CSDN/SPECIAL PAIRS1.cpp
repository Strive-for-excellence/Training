/*
FWT 求解
*/
const int maxn = 2e6+10;
LL A[maxn],B[maxn];
 
void FWT_and(LL *a,int N,int opt)
{
    for(int i=1;i<N;i<<=1)
        for(int p=i<<1,j=0;j<N;j+=p)
            for(int k=0;k<i;++k)
                if(opt==1)a[j+k]=(a[j+k]+a[i+j+k]);
                else a[j+k]=(a[j+k]-a[i+j+k]);
}
 
int main(void)
{
	int T;cin>>T;
	int N = 1<<20;
	while(T--){
		me(B);
		int n;cin>>n;
		for(int i = 1;i <= n; ++i){
			scanf("%lld",&A[i]);
			B[A[i]]++;
		}
		FWT_and(B,N,1);
		for(int i = 0;i < N; ++i)
			B[i] = B[i]*B[i];
		FWT_and(B,N,-1);
	
		cout<<B[0]<<endl;
 
	}
   return 0;
}
