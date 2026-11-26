# Referee report, slice 5 (SIGMA_TRANSFER.md, scripts/sigma_*), branch slice5, commits ff3f815..393bbdf

Hostile referee; scratch in /home/claude/conductor/review5/. Files under review not edited. Notes appended as I go.

## 0. Started 05:45 UTC
Read: SIGMA_TRANSFER.md, README.md, base_census.py, sigma_lvalues.py, sigma_lcert.py, outputs.
arXiv abstract 2604.28072 (WebFetch): title "A torsion-intersection proof of perfect-cuboid nonexistence on 1,072 explicit
master-tuple fibers"; abstract mentions a genus-3 reduction C_A: w^2 = lambda^8 + A lambda^4 + 1 (palindromic), elliptic
quotients E_A', E_A'', modular symbols L(E,1)/Omega via Sage, Kolyvagin. The abstract itself does not mention t -> 1/t,
E_PQ/E_uV/E_3 or Kani-Rosen (those are attributed to his section 4 / PDF; checked below if time).

## 4. The 17 "new" bases vs Peschmann's list (05:55 UTC)
Cloned renpe/euler-brick-obstructions (depth 1) into review5/pesch_repo. paper3/data/proven_fibers.csv: 1,072 rows.
- None of the 17 [(19,16) (40,33) (52,5) (56,41) (58,7) (63,40) (72,7) (77,26) (79,40) (80,17) (81,64) (83,68) (88,9)
  (89,40) (89,60) (96,49) (96,65)] is in it (checked both orders (m,n),(n,m)).
- All 14 README bases are in it, quotient E_3 with rank_resolution=modular_symbol (13) and E_uV (26,23): matches
  SIGMA_TRANSFER.md:66 and :90.
- His scan_full_M100.jsonl (2,040 rows) lists the 17 with his three quotients at rank lower bound >= 1 or ambiguous
  ([0,2]/[0,4] for E_uV on 9 of them), no L-ratio computed: consistent with "not closed by him", not with "closed by
  another means". STANDS as a list statement.

## 2/3. N, w and L(E,1): independent recomputation (05:52 UTC)
Own code, review5/myL.py (no PARI/Sage): a_p by numpy point counts at good odd p; at odd bad p (each checked to be
multiplicative: exactly one pair of the integer roots collides) a_p = Legendre(r_i - r_k) from the tangent cone; a_{2^k}
by hypothesis (f2, a2); N_odd = rad of the odd primes of the root differences; f2 and w CHOSEN by the functional equation
L(1) = G(x) + w G(1/x), G(x) = sum a_n/n e^{-2 pi n x/sqrt N}, compared at x = 1 and x = 1.25 (all f2 in 1..8, w = +-1).
| curve | my N (f2) | my w | my L(E,1) | their enclosure | FE residual (best / next) |
| (41,40) E_bd | 4035630 (1) | +1 | 3.9225616777 | [3.922561369, 3.922561986] | 4e-11 / 7e-4 |
| (19,16) E_ad | 19694640 (4) | +1 | 10.5322996312 | [10.532299318, 10.532299945] | 3e-11 / 1.5e-3 |
| (61,60) E_bd | 147371730 (1) | +1 | 3.5155968920 | [3.515596575, 3.515597209] | 4e-11 / 6e-4 |
| (40,23) E_bd | 174812190 (1) | +1 | 12.9857100582 | [12.985709741, 12.985710375] | 1e-11 / 3e-4 |
| (81,64) E_ad | 1260936240 (4) | +1 | 2.2105954882 | [2.210595170, 2.210595806] | 5e-11 / 4.5e-3 |
| (89,80) E_bd | 497081910 (1) | +1 | 2.7032644621 | [2.703264145, 2.703264780] | 3e-11 / 1.1e-3 |
All six agree with PARI's N and w and lie inside the sigma_lcert enclosures (agreement with the midpoint ~1e-12).
Root numbers from local data, review5/myw.py (no PARI): for the semistable curves (E_ab, E_bd, E_2 are Frey-shaped at 2:
roots 0, -alpha, beta with alpha = 1 mod 4... model Y^2+XY = X^3 + ((be-al-1)/4)X^2 - (al be/16)X; a_2 by counting mod 2;
good at 2 when v2(b^2) = 4), w = -prod(-a_p). Validated against PARI on all 945 such census curves with m < 40
(480 with w = -1): 945/945 agree in N and w. On the 15 semistable closing curves (13 E_bd, (56,41) E_ab, (26,23) E_2):
w = +1 and N = PARI's N, all 15. Conductors of (19,16),(72,7),(96,49) E_ad, (41,40) E_bd, (56,41) E_ab, (26,23) E_2 also by
Sage's own Tate implementation (algorithm='generic') and eclib (algorithm='mwrank'): equal to PARI's, all six.
For all 16 E_ad (additive at 2, f2 = 4): FE test of PARI's (N, w) with PARI's a_n (review5/fe_pari.py, float64):
residual <= 2e-14 for PARI's (N, w); flipping w gives >= 0.03, N -> 2N or N/2 gives >= 7e-4. So N, w are right
(numerically decisive) for all 31 closing curves: 15 by independent local root numbers, 16 by the FE test.

## 1. Rank 0 => closure (05:55 UTC)
- Fibre = census fibre: t = c != 0 rational with t^2+A, t^2+B, t^2+D squares (README:9-13, base_census.py:4-8);
  SIGMA_TRANSFER.md:7 uses the same. Maps C -> E: X = t^2, Y = t*y_x*y_y (E_ab, E_ad, E_bd), Y = t*y_a*y_b*y_d (E_2,
  quartic; cubic model x = ABD/X, checked by hand: (Y x^2/ABD)^2 = (x+BD)(x+AD)(x+AB)). Trivial identities, correct.
  A cuboid point has Y != 0 and X a nonzero square, so if E(Q) is finite and no point has X a nonzero square, closed.
- Torsion, own check (review5/mytors.py, own point counts, odd good p < 600): gcd #E(F_p) = 8 on all 34 closing curves
  except (26,23) E_2 (16), equal to the repo's orders. For order 8 the group is explicit: O, X in {0, -x^2, -y^2} and
  X = +-xy (checked on the curve), and xy (= m^4-n^4, 2mn(m^2+n^2), 2mn(m^2-n^2)) is not a square in any case (also
  Fermat); for (26,23) E_2 all 16 points re-checked on the curve, X = ABD/x: none a nonzero square. STANDS.
- sigma_lvalues.py:57-81 (torsion list + gcd bound, original model) is sound: L != 0 => E(Q) finite, gcd bound >= #E(Q),
  and gcd-many distinct rational points listed => list complete. sigma_lvalues.py controls are adequate.
- L(E,1) != 0 => rank 0: correct as quoted (modularity; GZ + Kolyvagin need an auxiliary imaginary quadratic twist with
  L' != 0, Bump-Friedberg-Hoffstein / Murty-Murty; or Kato directly). Suggest citing one of these.

## 2. sigma_lcert.py rigour
Formula (1+w) sum a_n/n e^{-2 pi n/sqrt N} for w = +1: correct (Cremona 2.11). q enclosure via mpmath iv at 200 bits;
mpf(q.a) of the degenerate endpoint is exact (checked). Fixed-point recursion rounds outward correctly (lines 50-56:
floor for lo, ceil for hi, sign handled). Tail: |a_n| <= d(n) sqrt n <= 2n, sum_{n>M} <= 2 q^{M+1}/(1-q); lines 58-60
give an upper bound in units 2^-K (note `>> K + 0` parses as `>> (K+0)`, fine); assert tail < 1e-6 (line 73). The
`certified` flag is integer-exact. Given (a_n, N, w) the enclosure is RIGOROUS (modulo trusting mpmath's iv exp/pi).
Issues:
 (i) NOT fail-safe in N or w (lines 65-66 take PARI's ellglobalred/ellrootno, nothing checks them): PARI w = +1 while the
     true w = -1, or a wrong N, would yield a confident "certified" nonsense interval. Acknowledged at
     SIGMA_TRANSFER.md:79 and sigma_lcert.py:18-19, but the table (SIGMA_TRANSFER.md:76) still reads as a certificate.
     Fix: enclose also G(x) + w G(1/x) at x = 5/4 (same integer arithmetic) and require overlap with x = 1, and
     disjointness for -w and for N' = N*p^{+-1}; or check w = -prod w_p from local data (done here, review5/myw.py).
     My checks (sec. 2/3 above) show N and w are in fact right for all 31.
 (ii) w = -1 case: returns not certified (line 66): fail-closed, correct.
 (iii) a_n from PARI up to 2.2e7, spot-checked only p < 1000 in the repo. My own a_n == PARI's ellan exactly up to the
     certificate's M on (41,40) E_bd, (19,16) E_ad, (61,60) E_bd, (40,23) E_bd [M up to 50509] (review5/cmp_an.py).
 (iv) L_lo/L_hi in the jsonl are float-printed (line 75, f'{lo/one:.12f}'), not rounded outward; display only.
 (v) SIGMA_TRANSFER.md:62 "lower ends from 0.50": certified minimum is 0.7105 ((96,49) E_ad); 0.504 is (80,27) E_ad,
     which sigma_lcert does not certify (--closing keeps only the smallest-N closing curve per base, line 99-102). So
     SIGMA_TRANSFER.md:66 "also by E_ad for (73,16),(80,27),(86,3),(97,80)" rests on Sage/PARI numerics only.

## 5. Claim (a) and Peschmann's involution
- Multiplicity one: sigma_multiplicity.py checks (no CM at inert primes for the nine class-number-one fields; |a_p|
  differ at a prime good for both => not Q-isogenous to any twist; a_p != 0 at p = 3 mod 4 separates E^(D), E^(-D)).
  The logic is right; a_p for the quartic E_2 uses -1-S (two points at infinity), right. Hom_Q(J,E) has rank = multiplicity
  (End_Q(E) = Z for every E/Q, so the no-CM clause is only needed over K, where it is needed: CM over K would double the
  rank of Hom_K). DM needs k > rank E(K) independent maps, so with multiplicity one it needs rank 0; over K containing i
  the twins give k = 2 but rank_K >= rank E + rank E^(-1). Argument STANDS within its stated scope (fibre and the
  first-descent covers; higher covers, Chabauty-type and elliptic-Chabauty arguments not covered, as the text says).
  Not re-verified by me: the genus-17 Kani-Rosen description of D_delta (SIGMA_TRANSFER.md:39-43).
  J(C) ~ E_ab x E_ad x E_bd x E_2 x E_5: checked (the 15 quadratic subcovers of the (Z/2)^4 cover of the X-line;
  genus 1 exactly for the five with 3 or 4 of X, X+A, X+B, X+D; 5 = g(C)). Twin identities checked by hand.
- Peschmann: arXiv PDF (WebFetch) and his repo (paper3/analysis/jh_full_decomp.sage:6-14): sigma_1: t -> -t,
  sigma_2: (t,v) -> (1/t, v/t^4), u = t+1/t for E_uV, w = t-1/t for E_3 (sigma_1 sigma_2), s = t^2 for E_PQ;
  J(H) ~ product cited as Kani-Rosen [KR89]; modular symbols via Sage L_ratio (restricted to semistable curves),
  Kolyvagin; hard fibres = all three ranks >= 1; Chabauty-Coleman said to apply only for total rank < 3, quadratic
  Chabauty named as next tool (sec. 5); Demjanenko-Manin and isogenies between quotients not discussed.
  SIGMA_TRANSFER.md:87 is accurate except: the involutions are defined in his sec. 3.2 (recalled from Paper 1) and used
  in sec. 4; "Chabauty only as future work" slightly understates sec. 5.1 (he states when Chabauty-Coleman applies).
- 17 new bases: for each, the closing curve (E_ad, or E_ab for (56,41)) has conductor != that of each of his E_PQ,
  E_uV, E_3 and an a_p differing at some p <= 41 (review5/pesch17.py, PARI): not isogenous to his quotients. STANDS.
- a_n check completed: own a_n == PARI ellan exactly (0 differences) up to the certificate's M on all six:
  (41,40) E_bd, (19,16) E_ad, (61,60) E_bd, (40,23) E_bd, (89,80) E_bd, (81,64) E_ad [M up to 141220].
- Merged records: 35 closing curves, all with both Sage modsym and PARI lfun; 6 zero curves, 2 of them Sage only
  ((65,12), (85,6) E_ab): consistent with SIGMA_TRANSFER.md:68.

## Verdict (05:57 UTC)
(a) STANDS, within its stated scope (fibre + first-descent covers); small result, multiplicity-one certificate sound.
(b) STANDS: all 31 closures re-verified (N, w independently, L recomputed on 6, torsion and square test independently,
    17 bases absent from Peschmann's CSV). Tier "computed, with quoted theorems" is right.
(c) GAP: enclosure rigorous given (a_n, N, w), but not fail-safe in N, w; misstatements below.
Must change before publishing:
 1. sigma_lcert.py:65-66 - add an internal check of N and w (interval FE test at x = 5/4 vs x = 1, and/or local root
    numbers), so a wrong PARI value cannot certify. (They are right for all 31; review5/myw.py, fe_pari.py.)
 2. SIGMA_TRANSFER.md:62 - "lower ends from 0.50": certified minimum is 0.7105, (96,49) E_ad; 0.504 is (80,27) E_ad,
    which is not certified.
 3. SIGMA_TRANSFER.md:66 - "also by E_ad for (73,16), (80,27), (86,3), (97,80)": not interval-certified
    (sigma_lcert.py:99-102 keeps one curve per base); say so or certify them.
 4. sigma_lcert.py:75 - endpoints float-printed, not rounded outward; print exact rationals.
 Minor: cite BFH/Murty-Murty or Kato (SIGMA_TRANSFER.md:60); Peschmann's involutions are his sec. 3.2, and sec. 5.1
 does say when Chabauty-Coleman applies (SIGMA_TRANSFER.md:87); "without CM" unnecessary over Q (:35).
 Not checked: genus-17 D_delta decomposition (:39-43), sigma_mobius.py.
