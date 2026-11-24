# Search code

- `search.cpp`: exhaustive search for planar difference families in Z/P with given block sizes, normalised as in Section 3 of the paper.
- `run.sh`: the full run for n = 28, P = 271 (blocks 13, 9, 7) in 40 chunks; resumable.
- `results.txt`: the last line of each chunk's log from the run reported in the paper.
- `recheck.cpp`: a second, separately written implementation of the same search; `stage_counts.txt` compares the number of admissible sets A per stage.
- `brute.cpp`: all families without normalisation (only 0 in each block), for small P; `brute_vs_search.txt` compares its classes with those of `search.cpp`.
- `classes.py`: classes of families under multipliers and translations.

Requirements: g++ with OpenMP. Set `OMP_NUM_THREADS` to limit the threads.

    g++ -O3 -march=x86-64-v2 -fopenmp search.cpp -o search

    # controls (seconds): PG(2,13) and PG(2,19)
    ./search 61 7 4 3 > f61.txt;   python classes.py 61 f61.txt     # 28 families, 12 classes
    ./search 127 9 7 4 > f127.txt; python classes.py 127 f127.txt   # 8 families, 4 classes

    # brute force vs search (example)
    g++ -O2 brute.cpp -o brute
    ./brute 29 5 3 2 > b.txt; ./search 29 5 3 2 > s.txt; python classes.py 29 b.txt s.txt   # 8 classes each, equal True

    # n = 28 (about 11 h on 4 threads)
    bash run.sh                    # 40 chunks, 0 families

    # re-check (resumable; about 16 h on 4-12 threads)
    g++ -O3 -march=x86-64-v2 -fopenmp recheck.cpp -o recheck
    ./recheck 271 13 9 7 out271.txt
    awk '{L+=$4; F+=$5} END {print L, F}' out271.txt   # 6871065175 0

Source licence: MIT (`LICENSE` in this directory).
