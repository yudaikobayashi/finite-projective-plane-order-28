// Exhaustive search for planar difference families {A,B,C} in Z_P with |A|=a, |B|=b, |C|=c:
// every nonzero element of Z_P is x-y for exactly one ordered pair x!=y in a common block.
// Normalisation (paper, Section 3): multipliers and translation of each block. A contains no 3-term
// arithmetic progression; the remaining triple classes l_1, l_2, ... (classes of l under l -> 1/l, 1-l) are
// processed in stages: in stage s, A contains {0,1,l_s} and no triple of class l_1..l_{s-1}.
// B (or C) contains {0,d}, d the smallest difference not covered by A; the other block contains {0,e} likewise.
// usage: search P a b c [sampleMod sampleRes]   (only jobs q with q % sampleMod == sampleRes; see run.sh)
// Families are printed to stdout, per-stage statistics to stderr; the last line of stderr is "DONE ...".
#include <bits/stdc++.h>
#include <omp.h>
using namespace std; typedef unsigned long long u64; const int W=5;
struct S{ u64 w[W]; };
int P,a,b,c; S FULLM;
static inline bool tst(const S&s,int i){ return s.w[i>>6]>>(i&63)&1; } static inline void st(S&s,int i){ s.w[i>>6]|=1ull<<(i&63);} static inline void cl(S&s,int i){ s.w[i>>6]&=~(1ull<<(i&63)); }
static inline int cnt(const S&s){ int r=0; for(int i=0;i<W;i++) r+=__builtin_popcountll(s.w[i]); return r; }
static inline S shl(const S&s,int k){ S r; int q=k>>6,b=k&63; for(int i=W-1;i>=0;i--){ u64 v=0; if(i-q>=0){ v=s.w[i-q]<<b; if(b&&i-q-1>=0) v|=s.w[i-q-1]>>(64-b);} r.w[i]=v; } return r; }
static inline S shr(const S&s,int k){ S r; int q=k>>6,b=k&63; for(int i=0;i<W;i++){ u64 v=0; if(i+q<W){ v=s.w[i+q]>>b; if(b&&i+q+1<W) v|=s.w[i+q+1]<<(64-b);} r.w[i]=v; } return r; }
static inline S rot(const S&s,int k){ if(!k) return s; S x=shl(s,k),y=shr(s,P-k); for(int i=0;i<W;i++) x.w[i]=(x.w[i]|y.w[i])&FULLM.w[i]; return x; }
struct Ctx{ u64 lev[40]={0}; u64 leaves=0,ab=0,fam=0; vector<string> out; };
template<class F> void ext(int*cur,int k,int K,const S&cand,const S&al,F&&fn){ // al: allowed (unused) differences
  if(k==K){ fn(al); return; } S cd=cand; int left=cnt(cd);
  for(int w=0;w<W;w++) while(cd.w[w]){ if(left--<K-k) return; int x=w*64+__builtin_ctzll(cd.w[w]); cd.w[w]&=cd.w[w]-1; S s=al; bool ok=true;
    for(int q=0;q<k;q++){ int d=x-cur[q]; if(d<0)d+=P; if(!tst(s,d)){ok=false;break;} cl(s,d); cl(s,P-d);} if(!ok) continue; cur[k]=x; S nc=cd; for(int t=0;t<=k;t++){ S r=rot(s,cur[t]); for(int i=0;i<W;i++) nc.w[i]&=r.w[i]; } ext(cur,k+1,K,nc,s,fn); } }
void finish(const int*A,const S&al,Ctx&C){ int d=1; while(!tst(al,d)) d++;
  for(int wh=0;wh<2;wh++){ int K=wh?c:b,K2=wh?b:c; if(K<2) continue; S s=al; int B[40]; B[0]=0;B[1]=d; cl(s,d); cl(s,P-d); S cand=s,r=rot(s,d); for(int i=0;i<W;i++) cand.w[i]&=r.w[i];
    ext(B,2,K,cand,s,[&](const S&r){ C.ab++; int D[40];
      auto emit=[&](){ C.fam++; string o; for(int i=0;i<a;i++) o+=to_string(A[i])+" "; o+="| "; for(int i=0;i<K;i++) o+=to_string(B[i])+" "; o+="| "; for(int i=0;i<K2;i++) o+=to_string(D[i])+" "; if(C.out.size()<100000) C.out.push_back(o); };
      int n=cnt(r); if(K2==1){ if(!n){ D[0]=0; emit(); } return; }
      int e=1; while(e<P&&!tst(r,e)) e++; if(e==P) return; S z=r; D[0]=0; D[1]=e; cl(z,e); cl(z,P-e); S c2=z,q=rot(z,e); for(int i=0;i<W;i++) c2.w[i]&=q.w[i];
      ext(D,2,K2,c2,z,[&](const S&zz){ if(!cnt(zz)) emit(); }); }); } }
struct Forb{ vector<S> F; bool any; }; // F[d]: {d*l : l in forbidden classes}
// add point x to A[0..k-1]; updates al (unused diffs) and cand
static inline void addpt(const int*A,int k,int x,S&al,S&cand,const Forb&fb){ int nd[40]; S ND; memset(&ND,0,sizeof ND);
  for(int q=0;q<k;q++){ int d=x-A[q]; if(d<0)d+=P; nd[q]=d; cl(al,d); cl(al,P-d); st(ND,d); st(ND,P-d); }
  S r=rot(al,x); for(int i=0;i<W;i++) cand.w[i]&=r.w[i];
  for(int t=0;t<k;t++){ S f=rot(ND,A[t]); for(int i=0;i<W;i++) cand.w[i]&=~f.w[i]; int m=A[t]+x; if(m&1) m+=P; m/=2; if(m>=P)m-=P; cl(cand,m);
    if(fb.any){ S g=rot(fb.F[nd[t]],A[t]); for(int i=0;i<W;i++) cand.w[i]&=~g.w[i]; } } }
void dfsA(int*A,int k,const S&cand,const S&al,const Forb&fb,Ctx&C){ C.lev[k]++; S cd=cand; int left=cnt(cd);
  if(k==a-1){ for(int w=0;w<W;w++) while(cd.w[w]){ int x=w*64+__builtin_ctzll(cd.w[w]); cd.w[w]&=cd.w[w]-1; S s=al; for(int q=0;q<k;q++){ int d=x-A[q]; if(d<0)d+=P; cl(s,d); cl(s,P-d);} A[k]=x; C.leaves++; finish(A,s,C);} return; }
  for(int w=0;w<W;w++) while(cd.w[w]){ if(left--<a-k) return; int x=w*64+__builtin_ctzll(cd.w[w]); cd.w[w]&=cd.w[w]-1; S s=al,nc=cd; addpt(A,k,x,s,nc,fb); if(cnt(nc)<a-k-1) continue; A[k]=x; dfsA(A,k+1,nc,s,fb,C); } }
int inv(int x){ int r=1,e=P-2,bb=x; while(e){ if(e&1) r=r*bb%P; bb=bb*bb%P; e>>=1;} return r; }
int main(int argc,char**argv){ P=atoi(argv[1]); a=atoi(argv[2]); b=atoi(argv[3]); c=atoi(argv[4]); int mod=argc>5?atoi(argv[5]):1,res=argc>6?atoi(argv[6]):0;
  memset(&FULLM,0,sizeof FULLM); for(int i=0;i<P;i++) st(FULLM,i);
  vector<vector<int>> cls; vector<char> seen(P,0); for(int l=2;l<P;l++) if(!seen[l]){ set<int> s={l}; bool ch=true; while(ch){ ch=false; for(int x:set<int>(s)) for(int y:{inv(x),(1-x+P)%P}) if(s.insert(y).second) ch=true; } for(int x:s) seen[x]=1; if(s.count(P-1)) continue; cls.push_back(vector<int>(s.begin(),s.end())); }
  // larger classes first (a class of size 2 is rarer)
  stable_sort(cls.begin(),cls.end(),[](const vector<int>&x,const vector<int>&y){return x.size()>y.size();});
  int s0=getenv("S0")?atoi(getenv("S0")):0,s1=getenv("S1")?atoi(getenv("S1")):(int)cls.size(); fprintf(stderr,"P=%d (%d,%d,%d) classes %zu\n",P,a,b,c,cls.size());
  u64 TF=0,TL=0; double t00=omp_get_wtime();
  for(int sg=s0;sg<s1&&sg<(int)cls.size();sg++){ Forb fb; fb.any=sg>0; fb.F.resize(P); for(auto&f:fb.F) memset(&f,0,sizeof f); for(int d=1;d<P;d++) for(int q=0;q<sg;q++) for(int l:cls[q]) st(fb.F[d],d*l%P);
    int base[3]={0,1,cls[sg][0]}; int A0[40]; S al=FULLM,cand=FULLM; cl(al,0); cl(cand,0); A0[0]=0; addpt(A0,1,1,al,cand,fb); cl(cand,1); A0[1]=1; if(!tst(cand,base[2])) continue; addpt(A0,2,base[2],al,cand,fb); A0[2]=base[2]; cl(cand,base[2]); cl(cand,0); cl(cand,1);
    vector<pair<int,int>> jobs; for(int x=2;x<P;x++) if(tst(cand,x)) for(int y=x+1;y<P;y++) if(tst(cand,y)) jobs.push_back({x,y});
    mt19937 rng(5); shuffle(jobs.begin(),jobs.end(),rng); Ctx T; double t0=omp_get_wtime();
    #pragma omp parallel for schedule(dynamic,1)
    for(size_t q=0;q<jobs.size();q++){ if((int)(q%mod)!=res) continue; int A[40]={A0[0],A0[1],A0[2]}; S s=al,nc=cand; int x=jobs[q].first,y=jobs[q].second; addpt(A,3,x,s,nc,fb); A[3]=x; if(!tst(nc,y)) continue; addpt(A,4,y,s,nc,fb); A[4]=y; for(int i=0;i<=y;i++) cl(nc,i); Ctx C; if(a==5){ C.leaves++; finish(A,s,C);} else dfsA(A,5,nc,s,fb,C);
      #pragma omp critical
      { for(int i=0;i<40;i++) T.lev[i]+=C.lev[i]; T.leaves+=C.leaves; T.ab+=C.ab; T.fam+=C.fam; for(auto&o:C.out) puts(o.c_str()); } }
    fflush(stdout); TF+=T.fam; TL+=T.leaves; fprintf(stderr,"stage %d l=%d nodes",sg,base[2]); for(int i=5;i<a;i++) fprintf(stderr," %llu",T.lev[i]); fprintf(stderr," leaves %llu AB %llu fam %llu t %.1f\n",T.leaves,T.ab,T.fam,omp_get_wtime()-t0); }
  fprintf(stderr,"DONE leaves %llu families %llu time %.1f\n",TL,TF,omp_get_wtime()-t00); }
