# G2: which primes empty the four-class set (Roadmap G, run 6)

Written by an AI model (Claude) in a separate session; reviewed by a separate AI agent (x6, `reviews/review_x6.md`, its
own checks in `../../referee_x6/`); not read by a person. **Experiment, not counted; tier D of `LEDGER.md`.**
`REPORT_G2.md` is the session's report, unchanged; the notes below take precedence over it.

## What it claims

- **Test 6 of `COMPLETE.md` §7 answered (computed): the bad primes alone do not always empty the four-class set.** With
  every bad prime (every prime of `2mn(m² − n²)`, of `m² + n²` and of the octic `m⁸ + 14m⁴n⁴ + n⁸`; the local test at odd
  primes by `../g1/exact_local.py`, at 2 the repository's), they leave tuples for **101 of the 997** and **836 of the
  8,739** of `class_selmer_15000_list.txt`. Each tuple left has an explicit 2-adic point
  (`g2_badalone_witness*.txt`). The review found, by separate code, explicit local points **at 2 and at every odd bad
  prime** for every tuple left (`../../referee_x6/badalone_witness_all.txt`), so this answer does not rest on the
  exactness of `exact_local.py`.
- Bad primes and good primes up to 101 (the 997) and 181 (the 8,739), all `≡ 1 (mod 4)`, empty all 9,736. The kill sets
  in `g2_*_all.txt` are irredundant (no prime can be dropped), not proved of least size. That good primes `≡ 3 (mod 4)`
  never help is forced: there every candidate class has the same label.
- No residue-symbol law found; the feature scan (`g2_features*.txt`) lists counts for every feature tried and states no
  conjecture.
- G1's 63 kept ratios (`g2_g1_kept_0.txt`): 57 emptied by bad primes alone, 3 need good primes (largest 61), 3 not run
  (out of memory at 8 GB).

## Review

- `g2_check.py`, `g2_features.py` (both lists), `g2_check.py --negative`: reproduced on PARI 2.17.4.
- 20 kill sets recomputed (`../../referee_x6/killsets_witness.txt`): empty with the repository's `tuple_descent.Local`
  (19 of the 20 need nothing else), and irredundant **with a proof for each drop**: after dropping any prime, one tuple
  left has explicit local points at every remaining prime.
- 10 "bad primes alone fail" ratios rechecked in the same file, then all 937 in `badalone_witness_all.txt`: the tuples
  left number as here for every ratio (6,904 in all), and each has explicit local points at 2 and at every odd bad prime.

## Changes when committed

The two input lists were byte-identical copies of the repository's `scripts/cover_descent/tuple_selmer_997.txt` and
`structure/class_selmer_15000_list.txt`; they are left out and `run_g2.sh` now reads the repository's files. The four
per-worker outputs `g2_tuple_selmer_997_{0,1}.txt`, `g2_class_selmer_15000_list_{0,1}.txt` are left out: sorted, they
are the lines of `g2_997_all.txt` and `g2_15000_all.txt` (checked).

## Rerunning

From this folder: `python3 g2_check.py`, `python3 g2_check.py --negative`, `python3 g2_features.py g2_997_all.txt`
(and `g2_15000_all.txt`), `python3 g2_badalone_witness.py g2_997_all.txt` (about 40 s), each asserting its committed
output; `bash run_g2.sh i N` (the 997: about 3 minutes per worker on two workers; the 8,739: about 57 minutes each) and
`bash run_g1kept.sh` (G1's 63, one process each under 8 GB) remake the kill sets. The scripts find the repository
through `MSS_REPO`, or else `../../..`, and `exact_local.py` in `../g1`.
