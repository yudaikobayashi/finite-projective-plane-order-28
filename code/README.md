# Search code

- `search.cpp`: exhaustive search for planar difference families in Z/P with given block sizes, normalised as in Section 3 of the paper.
- `run.sh`: the full run for n = 28, P = 271 (blocks 13, 9, 7) in 40 chunks; resumable.
- `results.txt`: the last line of each chunk's log from the run reported in the paper.
- `classes.py`: classes of the families found (used for the controls).

Requirements: g++ with OpenMP. Set `OMP_NUM_THREADS` to limit the threads.

    # controls (seconds): PG(2,13) and PG(2,19)
    g++ -O3 -march=x86-64-v2 -fopenmp search.cpp -o search
    ./search 61 7 4 3 > f61.txt;   python classes.py 61 f61.txt     # 28 families, 12 classes
    ./search 127 9 7 4 > f127.txt; python classes.py 127 f127.txt   # 8 families, 4 classes

    # n = 28 (about 11 h on 4 threads)
    bash run.sh                    # 40 chunks, 0 families

Source licence: MIT (`LICENSE` in this directory).
