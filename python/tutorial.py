"""
tutorial.py — Hochschild cohomology algebras & Baues–Wirsching cohomology,
explained by computation.  Run:  python3 tutorial.py

Every section states a piece of theory, computes it with hochschild.py,
and `check(...)`s the computation against the theorem.

© Computational Classnotes · CC BY-NC-SA 4.0
"""
import itertools
import numpy as np
from hochschild import *
from hochschild import _as_function

P0 = 32003            # stand-in for characteristic 0


def section(t):
    print("\n" + "=" * 78 + "\n" + t + "\n" + "=" * 78)


def say(*lines):
    for l in lines:
        print("  " + l)


def check(cond, msg):
    print(("  ✔ " if cond else "  ✘ FAILED: ") + msg)
    assert cond, msg


def signed(v, p):
    return [int(x) if x <= p // 2 else int(x) - p for x in v]


# ---------------------------------------------------------------------------
section("1. The Hochschild complex, written out for the dual numbers A = k[x]/(x²)")
# ---------------------------------------------------------------------------
say("C^n(A,A) = Hom_k(A^{⊗n}, A).  Normalizing (cochains vanish if any input is 1)",
    "lets us feed only the 'radical' basis {x}.  So a basis of C^n is",
    "(chain of n x's, output ∈ {1, x}):  dim C^n = 2 for every n.")
A = quiver_algebra("k[x]/x^2", [0], {'x': (0, 0)}, forbidden=[('x', 'x')])
cx = hochschild_complex(A, 5, P0)
for n in range(3):
    say(f"basis of C^{n}: {cx.basis[n]}")
say("", "The differential  (df)(a1..a_{n+1}) = a1 f(a2..) + Σ(-1)^i f(..a_i a_{i+1}..)"
        " + (-1)^{n+1} f(..a_n) a_{n+1}.")
say("Since x·x = 0, the middle terms vanish here, and only the outer terms survive:")
for n in range(3):
    say(f"d^{n} =\n      " + str(signed(cx.d[n].flatten(), P0)) + f"  (as a {cx.d[n].shape} matrix, rows=C^{n+1})")
say("Row 2 of d^1 says: for f(x) = 1,  (df)(x,x) = x·1 - f(x·x) + 1·x = 2x.  So d^1 ≠ 0 exactly",
    "when 2 ≠ 0 in k — this single entry '2' is the whole difference between char 0 and char 2.")
check(all(cx.check_d_squared()), "d∘d = 0 in every computed degree")
dims = [cx.betti(n) for n in range(6)]
say(f"dim HH^n, n=0..5 :  {dims}")
check(dims == [2, 1, 1, 1, 1, 1], "char 0:  HH^0 = A (commutative), HH^n ≅ k for n ≥ 1")

say("", "Reading the low degrees:")
_, H0 = cx.cohomology(0)
say(f"HH^0 = centre Z(A): cocycles {[dict((cx.basis[0][i][1], int(v)) for i, v in enumerate(c) if v) for c in H0.T]}")
_, H1 = cx.cohomology(1)
say(f"HH^1 = Der(A)/Inn(A): representative derivation {_as_function(cx, 1, H1[:, 0])}",
    "     i.e. D(x) = x — the Euler derivation x d/dx.  (Inner derivations are 0: A commutative.)")
_, H2 = cx.cohomology(2)
f2 = _as_function(cx, 2, H2[:, 0])
say(f"HH^2 = infinitesimal deformations: representative {f2}",
    "     i.e. μ_t(x, x) = x·x + t·f(x,x) = t·1 — the deformation  k[x]/(x² - t).",
    "     For t ≠ 0 that algebra is k×k (semisimple): the dual numbers are NOT rigid.")

# ---------------------------------------------------------------------------
section("2. The characteristic matters:  truncated polynomial rings k[x]/(x^n)")
# ---------------------------------------------------------------------------
say("Theorem: for A = k[x]/(x^n),  HH^0 = A,  and for i ≥ 1",
    "    HH^{2i-1} ≅ HH^{2i} ≅ k^{n-1}   if char k ∤ n,",
    "    HH^{2i-1} ≅ HH^{2i} ≅ k^{n}     if char k | n.",
    "Reason: the periodic resolution has maps 'multiply by x' and 'multiply by n x^{n-1}';",
    "the second vanishes exactly when p | n.")
for n, p in [(2, P0), (2, 2), (3, P0), (3, 3), (3, 2), (4, 2), (4, 3)]:
    T = quiver_algebra(f"k[x]/x^{n}", [0], {'x': (0, 0)}, forbidden=[('x',) * n])
    d, _ = HH_dims(T, 4, p)
    expect = [n] + [n if p != P0 and n % p == 0 else n - 1] * 4
    check(d == expect, f"k[x]/x^{n} over F_{p if p != P0 else '(char 0)'}: HH^0..4 = {d}")

# ---------------------------------------------------------------------------
section("3. HH*(A) is an algebra: the cup product")
# ---------------------------------------------------------------------------
say("(f ∪ g)(a1..a_{m+n}) = f(a1..a_m)·g(a_{m+1}..a_{m+n}).  It satisfies the Leibniz rule",
    "d(f∪g) = df∪g + (-1)^m f∪dg, so it descends to cohomology, where Gerstenhaber proved it is",
    "GRADED COMMUTATIVE:  [f]∪[g] = (-1)^{mn} [g]∪[f].")
cx = hochschild_complex(A, 6, P0)
# graded commutativity check on representatives
for m, n in [(1, 1), (1, 2), (2, 2), (0, 1), (1, 3)]:
    _, Hm = cx.cohomology(m)
    _, Hn = cx.cohomology(n)
    ok = True
    for i in range(Hm.shape[1]):
        for j in range(Hn.shape[1]):
            fg = cx.class_of(m + n, cup(cx, m, Hm[:, i], n, Hn[:, j]))
            gf = cx.class_of(m + n, cup(cx, n, Hn[:, j], m, Hm[:, i]))
            ok &= np.array_equal(fg % P0, ((-1) ** (m * n) * gf) % P0)
    check(ok, f"[f]∪[g] = (-1)^{m*n}[g]∪[f] on HH^{m} × HH^{n}")

# Name generators: HH^0 basis {1, x}; u ∈ HH^1; v ∈ HH^2.
say("", "Structure constants for k[x]/x² in char 0  (1, x ∈ HH^0;  u ∈ HH^1;  v ∈ HH^2):")
for m, n in [(0, 1), (0, 2), (1, 1), (1, 2), (2, 2), (1, 4), (2, 3)]:
    say(f"  HH^{m} × HH^{n} → HH^{m+n}:  {cup_table(cx, m, n).tolist()}")
say("Reading off:  x·u = 0, x·v = 0, u² = 0, u·v ≠ 0 spans HH^3, v² ≠ 0 spans HH^4, ...",
    "Hence   HH*(k[x]/x²) ≅ k[x, u, v] / (x², xu, xv, u²),   |x|=0, |u|=1, |v|=2   (char ≠ 2).")
t11 = cup_table(cx, 1, 1)
check(t11.sum() == 0, "u ∪ u = 0, as graded commutativity forces when 2 is invertible")
check(all(cup_table(cx, 2, k).any() for k in (2, 3, 4)), "v is not nilpotent: v∪(HH^k) ≠ 0 — HH* is infinite and finitely generated")

say("", "In characteristic 2 graded commutativity no longer kills squares of odd classes:")
cx2 = hochschild_complex(A, 4, 2)
say(f"  dims HH^0..4 over F_2: {[cx2.betti(n) for n in range(5)]}")
say(f"  HH^1 × HH^1 → HH^2 structure constants: {cup_table(cx2, 1, 1).tolist()}")
check(cup_table(cx2, 1, 1).any(), "char 2: some odd-degree class squares to a non-zero class")

# ---------------------------------------------------------------------------
section("4. The Gerstenhaber bracket:  HH^1 is a Lie algebra")
# ---------------------------------------------------------------------------
say("f ∘ g = Σ_i ±f(.., g(..), ..);  [f, g] = f∘g - (-1)^{(m-1)(n-1)} g∘f  makes HH^{*+1}",
    "a graded Lie algebra; on HH^1 it is the commutator of (outer) derivations.",
    "Example: the Kronecker algebra  1 ⇉ 2  (arrows a, b).  Its outer derivations act on the",
    "2-dimensional space span(a, b) by 2×2 matrices; scalars act innerly, so HH^1 ≅ gl_2/k ≅ pgl_2 ≅ sl_2 (char ≠ 2).")
K = quiver_algebra("Kronecker", [1, 2], {'a': (1, 2), 'b': (1, 2)})
cxK = hochschild_complex(K, 2, P0)
_, H1 = cxK.cohomology(1)
k = H1.shape[1]
ad = np.zeros((k, k, k), dtype=np.int64)
for i in range(k):
    for j in range(k):
        ad[i, j] = cxK.class_of(1, gerstenhaber_bracket(cxK, 1, H1[:, i], 1, H1[:, j]))
for i in range(k):
    say(f"  basis derivation h{i}: {_as_function(cxK, 1, H1[:, i])}")
for i in range(k):
    for j in range(i + 1, k):
        say(f"  [h{i}, h{j}] = {signed(ad[i, j], P0)}")
adm = [ad[i].T for i in range(k)]
Kf = np.array([[int(np.trace(adm[i] @ adm[j])) % P0 for j in range(k)] for i in range(k)])
check(k == 3, "dim HH^1(Kronecker) = 3")
check(rank_mod(Kf, P0) == 3, "Killing form nondegenerate ⇒ HH^1 is semisimple, 3-dimensional: sl_2")
check(rank_mod(ad.reshape(k * k, k), P0) == 3, "[HH^1, HH^1] = HH^1 (perfect Lie algebra)")

# ---------------------------------------------------------------------------
section("5. Quivers: Happel's formula for hereditary algebras")
# ---------------------------------------------------------------------------
say("For a path algebra kQ of a finite connected acyclic quiver:  HH^n = 0 for n ≥ 2 and",
    "   dim HH^1 = 1 - |Q_0| + Σ_{arrows α: i→j} dim e_j kQ e_i      (paths parallel to α).")
quivers = {
    "A3 linear 1→2→3": ([1, 2, 3], {'a': (1, 2), 'b': (2, 3)}),
    "A3 sink 1→2←3": ([1, 2, 3], {'a': (1, 2), 'b': (3, 2)}),
    "Kronecker 1⇉2": ([1, 2], {'a': (1, 2), 'b': (1, 2)}),
    "1⇉2→3": ([1, 2, 3], {'a': (1, 2), 'b': (1, 2), 'c': (2, 3)}),
    "triangle 1→2→3, 1→3": ([1, 2, 3], {'a': (1, 2), 'b': (2, 3), 'c': (1, 3)}),
}
for nm, (V, arr) in quivers.items():
    Q = quiver_algebra(nm, V, arr)
    d, _ = HH_dims(Q, 3, P0)
    happel = 1 - len(V) + sum(len(Q.homs(s, t)) for (s, t) in arr.values())
    check(d == [1, happel, 0, 0], f"{nm:22s} HH^0..3 = {d}   (Happel predicts HH^1 = {happel})")
say("Trees are rigid (HH^1 = 0); multiple arrows and 'shortcut' paths create outer derivations.")

say("", "Adding relations can create HH^2 — e.g. the oriented 2-cycle with both length-2 paths zero:")
Z = quiver_algebra("2-cycle, rad^2=0", [1, 2], {'a': (1, 2), 'b': (2, 1)}, forbidden=[('a', 'b'), ('b', 'a')])
d, _ = HH_dims(Z, 5, P0)
say(f"  HH^0..5 = {d}   (not hereditary: infinitely many non-zero HH^n)")

# ---------------------------------------------------------------------------
section("6. Small categories and their category algebras kC")
# ---------------------------------------------------------------------------
say("A finite category C gives the algebra kC: basis = morphisms, product = composition",
    "(0 if not composable), identity = Σ 1_x.  Groups, monoids, posets, quivers are all special cases.")
C2 = group_category("C2", ['1', 'g'], lambda a, b: '1' if a == b else 'g', '1')
for p in (P0, 3, 2):
    d, _ = HH_dims(C2.algebra(), 4, p)
    say(f"  k[C2] over F_{p if p != P0 else '(0)'}: HH^0..4 = {d}")
check(HH_dims(C2.algebra(), 4, 2)[0] == HH_dims(A, 4, 2)[0], "F_2[C2] ≅ F_2[x]/x² (x = 1+g): same HH")
check(HH_dims(C2.algebra(), 4, 3)[0] == [2, 0, 0, 0, 0], "Maschke: p ∤ |G| ⇒ kG separable ⇒ HH^{>0} = 0")

say("", "Centralizer decomposition:  HH*(kG) ≅ ⊕_{[g] conj. class} H*(C_G(g); k).",
    "For S3 (classes: 1, transpositions, 3-cycles; centralizers S3, C2, C3):")
els = list(itertools.permutations(range(3)))
nm_ = {g: ''.join(map(str, g)) for g in els}
inv_ = {v: k for k, v in nm_.items()}
S3 = group_category("S3", [nm_[g] for g in els],
                    lambda f, g: nm_[tuple(inv_[f][inv_[g][i]] for i in range(3))], '012')
C3 = group_category("C3", ['0', '1', '2'], lambda a, b: str((int(a) + int(b)) % 3), '0')
for p in (3, 2):
    hh, _ = HH_dims(S3.algebra(), 3, p)
    gS3, _ = BW_dims(S3, constant_system(S3), 3, p)
    gC2, _ = BW_dims(C2, constant_system(C2), 3, p)
    gC3, _ = BW_dims(C3, constant_system(C3), 3, p)
    pred = [a + b + c for a, b, c in zip(gS3, gC2, gC3)]
    say(f"  F_{p}:  H*(S3)={gS3}  H*(C2)={gC2}  H*(C3)={gC3}")
    check(hh == pred, f"F_{p}[S3]: HH^0..3 = {hh} = sum of centralizer cohomologies")

# ---------------------------------------------------------------------------
section("7. Baues–Wirsching cohomology of a small category")
# ---------------------------------------------------------------------------
say("The factorization category F(C): objects are morphisms f of C; a map f → f' is a pair (u, v)",
    "with f' = u∘f∘v.  A natural system is a functor D: F(C) → k-mod  (u_* = D(u,1), v^* = D(1,v)).",
    "Cochains: F^n(C;D) = Π over chains x0 ←f1— x1 ←f2— ... ←fn— xn of D(f1∘...∘fn), with",
    "  (δF)(f1..f_{n+1}) = f1_* F(f2..) + Σ(-1)^i F(.., f_i∘f_{i+1}, ..) + (-1)^{n+1} f_{n+1}^* F(f1..fn).",
    "Special natural systems:",
    "  • constant D = k           →  H_BW(C; k) = H*(BC; k), cohomology of the nerve;",
    "  • D(f) = M(cod f), M a functor C → k-mod   →  H*(C; M) = Ext_{kC}(k, M)  (functor cohomology);",
    "  • D(f) = k·Hom_C(dom f, cod f),  u_*g = u∘g, v^*g = g∘v   →  H_BW(C; D) ≅ HH*(kC).")
say("", "(a) Constant coefficients see the homotopy type of the nerve:")
crown = poset_category("crown", ['a', 'b', 'c', 'd'], lambda x, y: x == y or (x in 'ab' and y in 'cd'))
square = poset_category("square", ['0', 'x', 'y', '1'], lambda u, v: u == v or u == '0' or v == '1')
KronC = FiniteCategory("Kron", [1, 2], ['i1', 'i2', 'a', 'b'],
                       {'i1': 1, 'i2': 2, 'a': 1, 'b': 1}, {'i1': 1, 'i2': 2, 'a': 2, 'b': 2},
                       lambda f, g: g if f in ('i1', 'i2') else f, {1: 'i1', 2: 'i2'})
S2 = poset_category("S² poset", ['a', 'b', 'c', 'd', 'e', 'f'],
                    lambda x, y: x == y or ('ab'.find(x) >= 0 and y in 'cdef') or (x in 'cd' and y in 'ef'))
for Cc, exp, what in [(square, [1, 0, 0, 0], "square has a minimum: contractible"),
                      (crown, [1, 1, 0, 0], "crown's order complex is a 4-cycle ≃ S¹"),
                      (KronC, [1, 1, 0, 0], "two parallel arrows: nerve ≃ S¹"),
                      (S2, [1, 0, 1, 0], "three levels of two points: join S⁰*S⁰*S⁰ ≃ S²"),
                      (C2, None, None)]:
    if exp is None:
        d0, _ = BW_dims(C2, constant_system(C2), 4, 2)
        check(d0 == [1] * 5, f"one-object category C2: H*(BC2;F_2)=H*(RP^∞;F_2) = {d0}")
        continue
    d, _ = BW_dims(Cc, constant_system(Cc), 3, P0)
    check(d == exp, f"{Cc.name:9s} H_BW(C; k) = {d}   — {what}")

say("", "(b) Functor coefficients: C2 acting on k by the sign representation (p odd):")
sign = {'dim': {'*': 1}, 'mat': {'1': np.array([[1]]), 'g': np.array([[P0 - 1]])}}
d, _ = BW_dims(C2, functor_system(C2, sign), 4, P0)
check(d == [0, 0, 0, 0, 0], f"H*(C2; sign) = {d}  (no invariants, and |C2| invertible)")

say("", "(c) The comparison theorem  HH*(kC) ≅ H*_BW(C; k Hom(-,-)):")
say("    Normalized Hochschild cochains on composable non-identity morphisms (f1..fn) take values in",
    "    k·Hom(dom fn, cod f1) — exactly D(f1∘...∘fn).  The Hochschild and BW differentials agree term",
    "    by term; normalization doesn't change cohomology.  Numerically:")
for Cc in [C2, C3, crown, square, KronC, S2, S3]:
    for p in ([P0, 2] if Cc in (C2, S3) else [P0]):
        N = 3 if Cc is not S3 else 2
        hh, _ = HH_dims(Cc.algebra(), N, p)
        bw, cxb = BW_dims(Cc, parallel_system(Cc), N, p)
        check(hh == bw and all(cxb.check_d_squared()),
              f"{Cc.name:9s} F_{p if p != P0 else '(0)'}: HH = {hh} = H_BW = {bw}"
              f"   (BW complex dims {[cxb.dim(n) for n in range(N + 1)]})")

# ---------------------------------------------------------------------------
section("8. Gerstenhaber–Schack: incidence algebras of posets")
# ---------------------------------------------------------------------------
say("For a finite poset P with incidence algebra kP (= category algebra of P):",
    "   HH*(kP) ≅ H*(Δ(P); k)   — simplicial cohomology of the order complex.",
    "So for posets Hochschild cohomology is a homotopy invariant of the nerve, and both of the",
    "natural systems from §7 give the same answer:")
for Pp in [square, crown, S2]:
    hh, _ = HH_dims(Pp.algebra(), 3, P0)
    nerve, _ = BW_dims(Pp, constant_system(Pp), 3, P0)
    check(hh == nerve, f"{Pp.name:9s} HH*(kP) = {hh} = H*(ΔP) = {nerve}")
say("Every finite simplicial complex is (up to subdivision) an order complex, so every finite",
    "simplicial cohomology arises as a Hochschild cohomology of a finite-dimensional algebra.")

# cup product on the S² poset: H^2 generator squared must vanish (degree 4 > dim)
cxs = hochschild_complex(S2.algebra(), 4, P0)
check(cxs.betti(4) == 0, "S² poset: HH^4 = 0, so the cup square of the HH^2 class is 0 (as in H*(S²))")

print("\nAll checks passed.  © Computational Classnotes · CC BY-NC-SA 4.0")
