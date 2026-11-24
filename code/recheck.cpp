// Independent re-implementation of the search for planar difference families {A,B,C} in Z_P,
// |A|=a > |B|=b > |C|=c, with the normalisation of Section 3 of the paper.
// Written from scratch (std::bitset, per-job output); does not share code with search.cpp.
// usage: recheck P a b c outfile      (resumable: jobs already in outfile are skipped)
// Each output line: "s x y leaves families"  (stage s; x < y the two least elements of A outside {0,1,l_s}).
#include <bitset>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <set>
#include <array>
#include <algorithm>
#include <omp.h>
using namespace std;
const int MAXP = 271;
typedef bitset<MAXP> BS;
int P, a, b, c;
BS ALL;

static BS rotl(const BS& s, int k) { // {x+k : x in s}
  k %= P; if (k < 0) k += P; if (!k) return s;
  BS r = (s << k) | (s >> (P - k));
  return r & ALL;
}
static int md(long long x) { x %= P; if (x < 0) x += P; return (int)x; }

struct Block { int el[40]; int k; BS cand; BS rem; };

// add z to block X (X.cand must contain z); F (may be null): F[d] = {d*l : l in forbidden classes}
static void addElem(Block& X, int z, const vector<BS>* F) {
  BS nd;
  for (int i = 0; i < X.k; i++) { int d = md(z - X.el[i]); X.rem.reset(d); X.rem.reset(P - d); nd.set(d); nd.set(P - d); }
  X.cand &= rotl(X.rem, z);                       // z' - z must be an unused difference
  for (int i = 0; i < X.k; i++) {
    X.cand &= ~rotl(nd, X.el[i]);                 // z' - el[i] must not be a new difference
    long long m = X.el[i] + z; if (m & 1) m += P; X.cand.reset(md(m / 2)); // midpoints
    if (F) X.cand &= ~rotl((*F)[md(z - X.el[i])], X.el[i]);               // forbidden triples
  }
  X.cand.reset(z);
  X.el[X.k++] = z;
}

static Block start(const BS& rem) { Block X; X.k = 0; X.rem = rem; X.cand = rem; X.cand.reset(0); X.el[X.k++] = 0; return X; }

static int firstBit(const BS& s) { for (int i = 1; i < P; i++) if (s[i]) return i; return -1; }

// enumerate all completions of X to size K (new elements in increasing order); call f(X) on each
template <class Fn> static void grow(const Block& X, int K, const vector<BS>* F, Fn&& f) {
  if (X.k == K) { f(X); return; }
  if ((int)X.cand.count() < K - X.k) return;
  BS it = X.cand;
  if (X.k == K - 1) { // last element: every candidate is admissible; only the differences are needed
    for (int z = 1; z < P; z++) {
      if (!it[z]) continue;
      Block Y; Y.k = X.k; Y.rem = X.rem;
      for (int i = 0; i < X.k; i++) { int d = md(z - X.el[i]); Y.rem.reset(d); Y.rem.reset(P - d); Y.el[i] = X.el[i]; }
      Y.el[Y.k++] = z;
      f(Y);
    }
    return;
  }
  for (int z = 1; z < P; z++) {
    if (!it[z]) continue;
    it.reset(z);
    if ((int)it.count() < K - X.k - 1) return;
    Block Y = X; Y.cand = it; Y.cand.set(z);   // only elements > z remain (plus z itself, removed in addElem)
    addElem(Y, z, F);
    grow(Y, K, F, f);
  }
}

static long long completions(const BS& rem) { // number of {B,C} completing A, normalised
  long long fam = 0; int d = firstBit(rem);
  for (int w = 0; w < 2; w++) {
    int K = w ? c : b, K2 = w ? b : c;
    Block X = start(rem); if (!X.cand[d]) continue; addElem(X, d, nullptr);
    grow(X, K, nullptr, [&](const Block& Xd) {
      BS r = Xd.rem; if ((int)r.count() != K2 * (K2 - 1)) return;
      int e = firstBit(r); if (e < 0) return;
      Block Y = start(r); if (!Y.cand[e]) return; addElem(Y, e, nullptr);
      grow(Y, K2, nullptr, [&](const Block& Yd) { if (Yd.rem.none()) { fam++;
        #pragma omp critical
        { fprintf(stderr, "FAMILY found\n"); } } });
    });
  }
  return fam;
}

int inv(int x) { for (int y = 1; y < P; y++) if ((long long)x * y % P == 1) return y; return -1; }

int main(int argc, char** argv) {
  if (argc < 6) { fprintf(stderr, "usage: recheck P a b c outfile\n"); return 1; }
  P = atoi(argv[1]); a = atoi(argv[2]); b = atoi(argv[3]); c = atoi(argv[4]); const char* out = argv[5];
  if (P > MAXP) return 1;
  for (int i = 0; i < P; i++) ALL.set(i);
  // classes of l (l != 0,1) under l -> 1/l, 1-l; skip the class of -1 (3-term progressions)
  vector<vector<int>> cls; vector<int> done(P, 0);
  for (int l = 2; l < P; l++) {
    if (done[l]) continue;
    int i = inv(l), o = md(1 - l), oi = inv(o);
    set<int> s = {l, i, o, oi, md((long long)l * inv(md(l - 1))), md((long long)md(l - 1) * i)};
    for (int x : s) done[x] = 1;
    if (s.count(P - 1)) continue;
    cls.push_back(vector<int>(s.begin(), s.end()));
  }
  stable_sort(cls.begin(), cls.end(), [](const vector<int>& x, const vector<int>& y) { return x.size() > y.size(); });
  fprintf(stderr, "classes %zu\n", cls.size());
  // resume
  set<array<int, 3>> have;
  if (FILE* f = fopen(out, "r")) { int s, x, y; long long L, Fm; while (fscanf(f, "%d %d %d %lld %lld", &s, &x, &y, &L, &Fm) == 5) have.insert({s, x, y}); fclose(f); }
  FILE* fo = fopen(out, "a");
  for (int s = 0; s < (int)cls.size(); s++) {
    vector<BS> F(P);
    for (int d = 1; d < P; d++) for (int q = 0; q < s; q++) for (int l : cls[q]) F[d].set(md((long long)d * l));
    int lam = cls[s][0];
    BS full = ALL; full.reset(0);
    Block A = start(full);
    addElem(A, 1, &F);
    if (!A.cand[lam]) { fprintf(stderr, "stage %d: l not admissible\n", s); continue; }
    addElem(A, lam, &F);
    // jobs: x < y, the two least elements of A outside {0,1,lam}
    vector<pair<Block, array<int, 2>>> jobs;
    for (int x = 2; x < P; x++) {
      if (!A.cand[x]) continue;
      Block B = A; BS keep = A.cand; for (int t = 0; t <= x; t++) keep.reset(t);
      B.cand = keep; B.cand.set(x); addElem(B, x, &F);
      for (int y = x + 1; y < P; y++) {
        if (!B.cand[y]) continue;
        Block C = B; BS k2 = B.cand; for (int t = 0; t <= y; t++) k2.reset(t);
        C.cand = k2; C.cand.set(y); addElem(C, y, &F);
        if (!have.count({s, x, y})) jobs.push_back({C, {x, y}});
      }
    }
    fprintf(stderr, "stage %d jobs %zu\n", s, jobs.size());
    #pragma omp parallel for schedule(dynamic, 1)
    for (size_t j = 0; j < jobs.size(); j++) {
      long long leaves = 0, fam = 0;
      grow(jobs[j].first, a, &F, [&](const Block& Ad) { leaves++; fam += completions(Ad.rem); });
      #pragma omp critical
      { fprintf(fo, "%d %d %d %lld %lld\n", s, jobs[j].second[0], jobs[j].second[1], leaves, fam); fflush(fo); }
    }
  }
  fclose(fo);
  fprintf(stderr, "DONE\n");
}
