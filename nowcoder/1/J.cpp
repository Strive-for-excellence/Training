```

const int maxn = 30000+1000;
const int maxq = 200000+10;
struct Query{
 int l,r,id;
};
bool operator < (const Query &a,const Query &b){
 return a.r < b.r;
}
Query query[maxq];
int ans[maxq];
map<int,int> mp;
int Tree[maxn],a[maxn];
void Add(int x,int p){
   while(x < maxn){
       Tree[x] += p;
       x += lowbit(x);
   }
}
int Sum(int i){
  int ans = 0;
  while(i > 0){
    ans += Tree[i];
    i -= lowbit(i);
  }
  return ans;
}
int main(void)
{
   int n,m;
//   cout<<Sum(0)<<endl;
   while(cin>>n>>m){
    mp.clear();
    me(Tree);
//     cin>>n;
   for(int i = 1;i <= n; ++i)
     scanf("%d",&a[i]);
   for(int i = n+1;i <= 2*n; ++i)
      a[i] = a[i-n];
//   cin>>m;

   for(int i = 1;i <= m; ++i){
        int a,b;
    scanf("%d %d",&a,&b),query[i].id = i;
    query[i].l = b;
    query[i].r = a+n;
   }
   n *= 2;

   sort(query+1,query+1+m);
   int pre = 1;

   for(int i = 1;i <= m; ++i){
    for(int j = pre;j <= query[i].r; ++j){
     if(mp[a[j]] != 0){
        Add(mp[a[j]],-1);
     }
     Add(j,1);
     mp[a[j]] = j;
    }
    pre = query[i].r + 1;
    ans[query[i].id] = Sum(query[i].r)-Sum(query[i].l-1);
   }
   for(int i = 1;i <= m; ++i)
      printf("%d\n",ans[i]);

   }
   return 0;
}




```

