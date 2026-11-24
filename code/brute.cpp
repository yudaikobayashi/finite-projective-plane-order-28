// Brute force, no normalisation except: each block contains 0 (translation).
// Lists all planar difference families {A,B,C} in Z_P with |A|=a,|B|=b,|C|=c (a>b>c, so blocks are distinguishable).
// usage: brute P a b c > out.txt   (lines "A | B | C", same format as search.cpp)
#include <cstdio>
#include <cstdlib>
int P, sz[3], el[3][40], used[400];
long long cnt = 0;
void block(int w, int k, int from);
void emit() {
  cnt++;
  for (int w = 0; w < 3; w++) { for (int i = 0; i < sz[w]; i++) printf("%d ", el[w][i]); if (w < 2) printf("| "); }
  printf("\n");
}
// choose element k of block w, from 'from' upward
void block(int w, int k, int from) {
  if (w == 3) { for (int d = 1; d < P; d++) if (!used[d]) return; emit(); return; }
  if (k == sz[w]) { block(w + 1, 1, 1); return; }
  for (int x = from; x < P; x++) {
    bool ok = true; int i;
    for (i = 0; i < k; i++) { int d = (x - el[w][i] + P) % P; if (used[d] || used[P - d] || d == P - d) { ok = false; break; } used[d] = used[P - d] = 1; }
    if (ok) { el[w][k] = x; block(w, k + 1, x + 1); }
    for (int j = 0; j < i; j++) { int d = (x - el[w][j] + P) % P; used[d] = used[P - d] = 0; }
  }
}
int main(int argc, char** argv) {
  P = atoi(argv[1]); for (int i = 0; i < 3; i++) sz[i] = atoi(argv[2 + i]);
  for (int w = 0; w < 3; w++) el[w][0] = 0;
  block(0, 1, 1);
  fprintf(stderr, "families %lld\n", cnt);
}
