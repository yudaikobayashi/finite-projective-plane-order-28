#!/bin/bash
# Full search for n = 28, p = 271 (blocks 13, 9, 7) in 40 chunks; resumable (finished chunks are skipped).
g++ -O3 -march=x86-64-v2 -fopenmp search.cpp -o search || exit 1
mkdir -p run
for r in $(seq 0 39); do
  [ -f run/$r.done ] && continue
  ./search 271 13 9 7 40 $r > run/$r.fam 2> run/$r.log && grep -q DONE run/$r.log && touch run/$r.done
done
ls run/*.done | wc -l      # 40
cat run/*.fam | wc -l      # 0 (no family)
