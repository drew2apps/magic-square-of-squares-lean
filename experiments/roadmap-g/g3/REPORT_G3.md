# G3: census along the degenerate-type family, 11000 < m ≤ 13000

> **Read `NOTES.md` first.** This is the session's report, committed unchanged; where they differ, `NOTES.md` (after review x6) holds.

Written by an AI model (Claude); not read by a person. Experiment, not counted. Nothing here changes the census, which stops at m ≤ 1000.

Run at 062c54f (0.6.0 release candidate); the programs used are identical at 24d6ef6. Real sympy 1.14.0 is installed. Machine: PARI 2.19.1 and cypari2, not CI's PARI 2.17.4.

## Bottom line

- **Two blocks are done, 11000 < m ≤ 13000: 611 ratios, all closed except 4, with 3 more conditional on a rank being 1.** CENSUS_TABLE_G3.tsv gives, for each of the 611 ratios of square or twice-square type that pass the congruences, the first step that closes it and that step's label. Format and labels are those of F0's CENSUS_TABLE_1000_11000.tsv (1000 < m ≤ 11000, run 5), which this continues for these two types. *(computed; PARI rank bounds; labels per row)*
- **Open:** (11165,7548), (12201,7360), (12204,9853) and (12703,2767). The last is the twice-square partner of (7735,4968), the one open ratio of F0's table, so it was expected.
- **Conditional on a rank being 1:** (11441,1519), (11969,4529) and (12007,1799), each closed only through a partner that F0's table already labels conditional.
- **Scope:** only the two degenerate types. Ratios of other type with m > 11000 would need the class step first, which G3 does not run. The run goes on toward 30,000 (fam_11000_30000.json has 5,796 ratios, about 300 per block); each block takes 2.6 to 2.9 h of second descent on two cores plus about 11 min of later stages.

## The steps (same programs as CENSUS9000.md, unchanged)

1. **Family (gen_family.py):** coprime m : n with m² + n² a square or twice a square, 11000 < m ≤ 30000, enumerated from Pythagorean parameters, that pass the congruences of ../../census2000/gen_pairs.py (the 484 classes; the 5 mod 8 condition always holds for these types). 5,796 ratios.
2. **Second descent** (../../census-wide/wide_pari.py, unchanged; run_g3.sh), by blocks of 1000 in m.
3. **Rank-0 quotients** (the repository's scripts/g1_full/global.py, unchanged) and **Demjanenko–Manin** (../../census2000/dm_pari.py, generators saturated at primes ≤ 1000), on each block's open ratios together with their partners (pipe_g3.py, run_g3_stages.sh).
4. **Partners:** a ratio closes when its partner does, as in CENSUS9000.md. The partner of m : n is (m + n) : (m − n), halved when both are even, in lowest terms. A twice-square ratio's partner is of square type with a smaller m. When that partner is in F0's table (m ≤ 11000), its closing row is taken, with F0's label.
5. **Rank-1 route** (../../census-wide/route/route_pari.py, unchanged; route/run_route.sh): counted for square-type ratios only. CENSUS9000.md: the route is built for square type, and on a twice-square ratio its fibre step reports "no" rather than "degenerate". Run on a square-type leftover whose '+13' class d gives E1 rank bounds [1,1], found by ../../census-wide/route/scope_pari.py (route/scope_g3_<lo>.jsonl).

g3_table.py builds the table from these records and runs no curve computation.

## Results (g3_table.txt)

```
blocks: 11000-12000, 12000-13000 (11000 < m <= 13000); ratios: 611
square: second descent 324, rank 0 18, Demjanenko-Manin 18, special 1, open 3
twice: second descent 209, rank 0 21, Demjanenko-Manin 12, special 4, open 1
OPEN: [(11165, 7548), (12201, 7360), (12204, 9853), (12703, 2767)]
CONDITIONAL on a rank being 1: [(11441, 1519), (11969, 4529), (12007, 1799)]
closed by the route (conditions (a), (b)): [(12416, 5313), (12487, 4241)]
```

Per block (stage_<lo>.json):

| block | ratios | after second descent | after rank 0 | after DM | DM alerts |
|---|---|---|---|---|---|
| 11000–12000 | 308 | 34 | 16 | 3 | none |
| 12000–13000 | 303 | 44 | 23 | 6 | none |

## The ratios not closed outright

- **(11165,7548), square, OPEN.** Every step leaves it. Its '+13' class d = 135369842 gives E1 rank bounds [0,2], and −d gives [2,2] (route/scope_g3_11000.jsonl), so the route doesn't apply. Its partner (18713,3617) is in the family and is reached in block 18000.
- **(12201,7360) and (12204,9853), square, OPEN.** The route needs E1 of rank 1. For (12201,7360) the bounds are [1,3] for d and [0,2] for −d; for (12204,9853), [1,3] for both (route/scope_g3_12000.jsonl). Their partners, (19561,4841) and (22057,2351), are reached in blocks 19000 and 22000.
- **(12703,2767), twice, OPEN.** Its partner is (7735,4968), the one open ratio of F0's table, and CENSUS9000.md shows every usable curve of both has rank 2 or more.
- **(12416,5313), square, closed by the route:** d = 251858174, E1 rank [1,1] (both bounds, and the smallest upper bound in its isogeny class is 1), the generator saturated at primes ≤ 2000. The sieve at four primes leaves only t = ±1, which the fibre step reports as degenerate (route/log_12416.txt, 9 s). *(computed; conditional on (a) G generating E1(ℚ) modulo torsion, (b) the reading of Checcoli–Veneziano–Viada's constants, as in CENSUS9000.md)*
- **(12487,4241), twice:** closes through its partner (8364,4123), which CENSUS9000.md closes by the route, with that label.
- **Conditional on a rank being 1, through their partners in F0's table:** (11441,1519) via (6480,4961), (11969,4529) via (8249,3720), and (12007,1799) via (6903,5104). Each partner is conditional there (root number −1, PARI bounds [1,3]).

## What holds and how it was checked

- **Table (g3_table.py):** the default run checks g3_table.txt and CENSUS_TABLE_G3.tsv byte for byte.
- **Negative mode (g3_table.py --negative, g3_table_negative.txt):** for one ratio of each step, every record of it and of its partner is removed, and it must come out OPEN. 4 of 4 come out OPEN, as they must.
- **Cross-checks with run 5 and CENSUS9000.md:** (11969,4529) and (12703,2767) are the partners of (8249,3720) and (7735,4968). CENSUS9000.md already ran the route on both and set those runs aside as twice type. Here (11969,4529) is again closed only through (8249,3720), and (12703,2767) is open, as (7735,4968) is.
- **The route on a twice-square ratio:** route/log_11441.txt is a run on (11441,1519) with d = 410, E1 rank [1,1], which ends closed. It is **not counted**, for the reason in step 5. It is kept as a record, and route/scope_g3_11000.jsonl gives the class d and rank bounds for the leftovers of block 11000.

## What it rests on

The same tiers as CENSUS9000.md and F0's table: second descent with PARI rank upper bounds (not proved); global.py; dm_pari.py with PARI generators saturated at primes ≤ 1000; for rows taken from F0's table, F0's label. A row labeled CONDITIONAL rests on PARI's rank bounds [1,3] and the root number, as in CENSUS9000.md.

## What failed

- The first stage run for block 11000 stopped because MSS_REPO was not set, so global.py was not found. run_g3_stages.sh now sets MSS_REPO (default ../../..), and the block was rerun.
- census-wide's sd_pari.py, which scope_pari.py and route_pari.py import, finds the repository by an absolute path into my workspace (its line 9, as committed in 0.6.0). I ran them with that clone present and did not change them. To rerun the route elsewhere, that line has to point at a clone.
- I first read the route's "closed" on (11441,1519) as a closure. It is not one, for the reason in step 5.

## Next step

- Let the run continue to 30,000; the table and summary update by `python3 g3_table.py --write`.
- (11165,7548): its partner (18713,3617) is in the family and is reached in block 18000.
- For the conditional ratios, the open question is the same as for (8249,3720) and (6480,4961): a rank of 1 that PARI bounds only in [1,3].

## Proposed LEDGER / CHANGES text

> **Experiment, not counted.** Census of the two degenerate types (m² + n² a square or twice a square) beyond F0's table (experiments/roadmap-g/g3): 611 ratios with 11000 < m ≤ 13000 that pass the congruences, each with its closing step and label in CENSUS_TABLE_G3.tsv. Open: (11165,7548), (12201,7360), (12204,9853), (12703,2767). Conditional on a rank being 1, through partners already conditional: (11441,1519), (11969,4529), (12007,1799). (12416,5313) is closed by the rank-1 route, conditional as in CENSUS9000.md. The same programs as CENSUS9000.md; PARI rank bounds throughout. Only these two types; other ratios with m > 11000 are not covered. Computed. Written by an AI model; not read by a person.

## Proposed CI lines

```
cd experiments/roadmap-g/g3        # needs F0's table: experiments/census-wide/table/CENSUS_TABLE_1000_11000.tsv, or F0_TABLE=...
python3 g3_table.py | tail -1 | grep -q '^# matches'                # under 1 s
python3 g3_table.py --negative | tail -1 | grep -q '^# matches'     # under 1 s
# not for CI: bash run_g3.sh and bash run_g3_stages.sh (2.6 to 2.9 h of second descent per block of 1000 on two cores,
#             plus about 11 min of stages; times.txt)
```

pipe_g3.py and gen_family.py find the repository through MSS_REPO, or else `../../..`.
