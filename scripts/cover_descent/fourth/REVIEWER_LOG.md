# REPORT (tuple) - independent re-implementation of the four-class local test
Started. Findings are appended as I go.

## Step 1 (from prose only: LAST_PAIR 11,12,14; RATIONAL_DISTANCE 2) - definition I implement
Forms g_j(u,v)=a_j u+b_j v over Z[i], j=1..4.  A tuple (d_1..d_4) of squarefree positive integers is "possible at p" iff
there is a PRIMITIVE pair (u,v) in Z_p^2 with g_j(u,v) invertible in Q_p(x)Q(i) and d_j*g_j(u,v) a square there, all j.
Only the image of d_j in H_p = Q_p^*/(Q_p^* cap squares of Q_p(x)Q(i)) matters:
 inert p: H_p = Z/2 (v_p mod 2);  split p: (Z/2)^2 (v_p mod 2, Legendre of unit part);  p=2: (Z/2)^2 (v_2 mod 2, odd part = +-1 or +-3 mod 8), since -1=i^2 and 2,5,10 are not squares in Q_2(i).
Because all primitive pairs are allowed, the set is automatically closed under (c_1..c_4) -> (lam c_1..lam c_4), lam a unit class
(non-trivial only for split p [non-residue] and p=2 [class of 5]).  NOT closed under multiplying by p (pair must stay primitive).
Ambiguities in the prose:
 A1 sec.11 writes gamma_L=(1+i)((mu-nv)+i(mv-nu)), Lemma S3 writes (mu-nv)+i(nu-mv) (the conjugate). Harmless for rational classes (conj swaps components) but the text is inconsistent.
 A2 sec.12 "with gamma_L,gamma_U multiplied by 1+i" - sec.11 already has the factor; one reading doubles it. (1+i) is NOT a rational class times a square, so this matters. I take N(gamma_L)=GL, i.e. one factor.
 A3 unit closure is stated in sec.12 and sec.14 ("closed under the classes of units") but nowhere is the unit class group written down (which units are non-trivial at inert p / p=2); R3 says rational units are squares at inert p; for p=2 nothing is said about -1,5.
 A4 whether g_j(u,v)=0 (or zero in one component at split p) counts as "square" is not stated. I exclude it (a global point has non-zero norm).
 A5 rational distance: classes of the forms are D_j times sqf(s*v) ("up to the common factor sv"); supports stated as primes of 2rs / 2s(r-s).

## Step 2 (my implementation: mytuple.py, drive.py)
- U_p: exhaustion mod p^N with a table of ALL squares of Z[i]/p^N (no valuations), charts (t,1),(1,pt) times unit reps
  ({1,non-residue} odd p; {1,3,5,7} for p=2).  Proof that reps suffice: lam = rep*mu^2 mod p^N and mu^2 preserves squares.
  Checked against exhaustion over ALL primitive pairs: 360 random (system, prime power) cases, 0 mismatches.
- L_p: exact class in H_p of g_j(u,v) at explicit integer pairs (witnesses of the exhaustion + digit walk to zeros).
- 2-adic: unit squares mod 16,32,64 in Z[i] are exactly the lifts of the 4 unit squares mod 8 (verified by enumeration).
- second upper bound without the table (valuations/CRT) for deeper N: 2.4M elements + 72 systems vs table: 0 mismatches.
- (627,364), primes of support + all <=61: U leaves exactly (1,2,263,991),(1,2,1982,526) = prose claim "two tuples".

## Step 3 (compare.py -> compare_out.txt): repo tuple_descent.py Local.allowed() used as a black box, per prime
Systems: 17 magic ratios (7 named + 10 random from pairs157.txt, seed 11), 19 rational-distance lines, 3 cuboid bases
(third form doubled to make four). Primes: supports + all primes <= 61.  723 (system, prime) cases.
 L_p inside repo set R_p: 723/723  (NO rejected certified tuple = no bug found)
 R_p inside U_p:          723/723  (repo never weaker than my exhaustion)
 U_p inside R_p:          722/723 ; L=R=U in 715/723 ; L=U in 715.
Final global sets from my U_p equal repo descend(pmax=61) in all 39 systems (e.g. 627,364: the two known tuples;
820,799 / 889,320 / 899,724: empty; 26/225: 12 tuples, 8 non-constant; cuboid 15,8: 12; 3,2: 4; 2,1: 4).
The 8 cases with L != U were all weaknesses of MY witness search (e.g. (627,364) at p=991: the exhaustion mod 991^2 is
ambiguous at u=-v mod 991 since 991*g is 0 mod 991^2, and my walk only went deeper, not sideways; the missing tuple is
that of the known point t=-1).  After fixing the walk (start one level up), the gap cases are re-run below.

## Step 3 rerun after the walk fix (compare_out2.txt)
723 cases: L in R 723, R in U 723, U in R 722, L=R=U 722.  The one remaining case has L subset R, R strictly inside U
(my exhaustion at its N weaker than the repo there; not identified individually).  No rejected certified tuple.

## Step 4 (attack.py -> attack_out.txt; mut.py)
- 120 random systems (coefficients in [-5,5], zeros and rational zeros allowed) x 10 primes (2..41): 1200 cases, L in R, R in U, U in R all 1200; L=R=U 1185.
- 72 designed systems (exact rational zero of one form, second form vanishing mod p at the same point, p | coefficients,
  powers (1+i)^k as content) at p in {2,3,5,7,13,17}: 252 cases, no disagreement.
- 60 planted rational points (g_j(u0,v0) = c_j z_j^2): tuple accepted by repo at all 25 primes < 100: 1500/1500; plus 300 set comparisons.
- total step 4: 1752 cases, L in R 1752, R in U 1752, U in R 1752, all equal 1732 (the 20 others: my L smaller, U = R).
- mutations of the repo program (copy, never the original): unit closure removed at split p -> a certified tuple is
  rejected in 21/28 cases (first: (627,364), p=5, point (50,2)); removed at p=2 -> 3/28 (point (3,3));
  hint replaced by the wrong class -> 0/28 detected (the hint is used only at depth K and shallower discs already give the tuples).
- Not done: primes 67..150 with N=1 for all systems (only support primes > 70, at N=2 by the table-free exhaustion).

## Verdict
No disagreement found between tuple_descent.py and an independent implementation with certified lower bounds.
