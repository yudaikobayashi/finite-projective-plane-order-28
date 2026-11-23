# A projective plane of order 28 admits no collineation of order 271

> **Status.** Computer-assisted result, **not peer reviewed**; novelty not verified against the literature.
> The search was run with a single program and has not been repeated with an independent implementation.
> Most derivations and programs were made by AI (Claude, Anthropic) under the author's direction.
> This does **not** settle whether a projective plane of order 28 exists.

**Main result.** A projective plane of order 28 has no collineation of order 271. Such a collineation would give a planar difference family in Z/271 with blocks of sizes 13, 9, 7; an exhaustive search finds none.

**Corollary.** For every order n < 36 that is not a prime power, a projective plane of order n has no collineation of prime order p > n + 1.

- Paper: [`paper/main.pdf`](paper/main.pdf) (source `paper/main.tex`; build: `make -C paper`).
- Code: [`code/`](code/) (`search.cpp`, `run.sh`; per-chunk output in `code/results.txt`).
- Licence: paper and documents CC BY 4.0 (`LICENSE`), source code MIT (`code/LICENSE`). Citation: `CITATION.cff`.
