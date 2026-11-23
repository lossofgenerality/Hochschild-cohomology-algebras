"""
hochschild.py — Hochschild cohomology of algebras with a multiplicative basis,
and Baues–Wirsching cohomology of finite categories.

© Computational Classnotes · CC BY-NC-SA 4.0

Everything is exact linear algebra over the prime field F_p.
  * p = 32003 stands in for characteristic 0 (the answers below agree with
    the char-0 theorems; a large prime only fails to detect torsion that
    happens to divide p).
  * p = 2, 3, ... are used to see genuinely modular phenomena.

Conventions
-----------
A basis element b is drawn as an arrow  b : s(b) -> t(b).
The product a*b means "first b, then a" (composition a∘b); it is non-zero
only when s(a) = t(b).  Every algebra here has a *multiplicative basis*:
the product of two basis elements is another basis element or 0.
This covers
  - path algebras of quivers modulo monomial relations  (k[x]/x^n, A_n, Kronecker, ...)
  - category algebras kC of finite categories           (group algebras, monoid algebras,
                                                         incidence algebras of posets, ...)
"""
from __future__ import annotations
from itertools import product as iproduct
import numpy as np

# ---------------------------------------------------------------------------
# 1.  Linear algebra over F_p
# ---------------------------------------------------------------------------

def rref_mod(M, p):
    """Reduced row echelon form of M over F_p.  Returns (R, pivot_columns)."""
    R = np.array(M, dtype=np.int64) % p
    rows, cols = R.shape
    pivots, r = [], 0
    for c in range(cols):
        if r == rows:
            break
        nz = np.nonzero(R[r:, c])[0]
        if nz.size == 0:
            continue
        k = r + nz[0]
        if k != r:
            R[[r, k]] = R[[k, r]]
        R[r] = (R[r] * pow(int(R[r, c]), p - 2, p)) % p
        others = np.nonzero(R[:, c])[0]
        others = others[others != r]
        if others.size:
            R[others] = (R[others] - np.outer(R[others, c], R[r])) % p
        pivots.append(c)
        r += 1
    return R, pivots


def rank_mod(M, p):
    if M.size == 0:
        return 0
    return len(rref_mod(M, p)[1])


def nullspace_mod(M, p):
    """Columns spanning {x : Mx = 0} over F_p."""
    rows, cols = M.shape
    if rows == 0:
        return np.eye(cols, dtype=np.int64)
    R, piv = rref_mod(M, p)
    free = [c for c in range(cols) if c not in piv]
    N = np.zeros((cols, len(free)), dtype=np.int64)
    for j, f in enumerate(free):
        N[f, j] = 1
        for i, pc in enumerate(piv):
            N[pc, j] = (-R[i, f]) % p
    return N


def solve_mod(M, b, p):
    """One solution x of Mx = b over F_p, or None."""
    rows, cols = M.shape
    R, piv = rref_mod(np.hstack([M, b.reshape(-1, 1)]), p)
    if cols in piv:
        return None
    x = np.zeros(cols, dtype=np.int64)
    for i, pc in enumerate(piv):
        x[pc] = R[i, -1]
    return x


# ---------------------------------------------------------------------------
# 2.  Cochain complexes and their cohomology
# ---------------------------------------------------------------------------

class CochainComplex:
    """
    A finite stretch C^0 -> C^1 -> ... -> C^N of a cochain complex over F_p.
    `basis[n]` lists labels of a basis of C^n, `d[n]` is the matrix of
    d^n : C^n -> C^{n+1}  (shape  dim C^{n+1} x dim C^n).
    """

    def __init__(self, basis, d, p):
        self.basis, self.d, self.p = basis, d, p
        self.index = [{lab: i for i, lab in enumerate(B)} for B in basis]
        self._H = {}

    def check_d_squared(self):
        """d^{n+1} d^n = 0 for every n we have.  Returns list of booleans."""
        return [not np.any((self.d[n + 1] @ self.d[n]) % self.p)
                for n in range(len(self.d) - 1)]

    def dim(self, n):
        return len(self.basis[n])

    def betti(self, n):
        """dim H^n = dim C^n - rank d^n - rank d^{n-1}."""
        rk_out = rank_mod(self.d[n], self.p)
        rk_in = rank_mod(self.d[n - 1], self.p) if n > 0 else 0
        return self.dim(n) - rk_out - rk_in

    def cohomology(self, n):
        """
        Returns (Bmat, Hmat):  columns of Bmat span the coboundaries B^n,
        columns of Hmat are cocycles whose classes form a basis of H^n.
        """
        if n in self._H:
            return self._H[n]
        p = self.p
        Z = nullspace_mod(self.d[n], p)
        if n > 0:
            Bfull = self.d[n - 1]
            _, piv = rref_mod(Bfull, p)              # independent coboundaries
            Bmat = Bfull[:, piv]
        else:
            Bmat = np.zeros((self.dim(n), 0), dtype=np.int64)
        both = np.hstack([Bmat, Z])
        _, piv = rref_mod(both, p)
        nb = Bmat.shape[1]
        Hmat = both[:, [c for c in piv if c >= nb]]
        self._H[n] = (Bmat, Hmat)
        return self._H[n]

    def class_of(self, n, cocycle):
        """Coordinates of the class [cocycle] in the basis returned by cohomology(n)."""
        Bmat, Hmat = self.cohomology(n)
        x = solve_mod(np.hstack([Bmat, Hmat]), np.asarray(cocycle) % self.p, self.p)
        if x is None:
            raise ValueError("not a cocycle (or not in Z^n)")
        return x[Bmat.shape[1]:] % self.p


# ---------------------------------------------------------------------------
# 3.  Algebras with a multiplicative basis
# ---------------------------------------------------------------------------

class MultAlgebra:
    """
    labels : list of basis-element names
    src, tgt : dict label -> object (vertex)
    mul : function (a, b) -> label of a∘b, or None for 0   [only called when s(a)=t(b)]
    idem : dict object -> label of its idempotent / identity e_x
    """

    def __init__(self, name, labels, src, tgt, mul, idem):
        self.name, self.labels, self.src, self.tgt = name, labels, src, tgt
        self._mul, self.idem = mul, idem
        self.objects = list(idem.keys())
        self.E = set(idem.values())                     # the separable subalgebra kE
        self.rad = [b for b in labels if b not in self.E]   # basis of A/kE
        # hom-sets: (x, y) -> basis elements x -> y
        self.hom = {}
        for b in labels:
            self.hom.setdefault((src[b], tgt[b]), []).append(b)
        self._mtab = {}
        for a in labels:
            for b in labels:
                if src[a] == tgt[b]:
                    self._mtab[(a, b)] = mul(a, b)

    def mul(self, a, b):
        return self._mtab.get((a, b))

    def dim(self):
        return len(self.labels)

    def homs(self, x, y):
        return self.hom.get((x, y), [])

    def __repr__(self):
        return f"<{self.name}: dim {self.dim()}, {len(self.objects)} object(s)>"


# ---- constructors ---------------------------------------------------------

def quiver_algebra(name, vertices, arrows, forbidden=(), max_len=12):
    """
    Path algebra kQ/I, I generated by the monomial paths in `forbidden`.
    arrows: dict name -> (source, target).  A path is a tuple of arrow names
    written right-to-left like composition: ('b','a') means first a then b.
    """
    forb = [tuple(f) for f in forbidden]

    def ok(path):
        for f in forb:
            L = len(f)
            for i in range(len(path) - L + 1):
                if path[i:i + L] == f:
                    return False
        return True

    paths = {('e', v): (v, v) for v in vertices}      # trivial paths
    frontier = [((a,), arrows[a][0], arrows[a][1]) for a in arrows]
    while frontier:
        new = []
        for path, s, t in frontier:
            if not ok(path):
                continue
            paths[path] = (s, t)
            if len(path) < max_len:
                for a, (sa, ta) in arrows.items():
                    if sa == t:
                        new.append(((a,) + path, s, ta))
        frontier = new
        if len(paths) > 2000:
            raise ValueError("algebra infinite-dimensional (or too big); add relations")

    def lab(path):
        if path[0] == 'e':
            return 'e' + str(path[1])
        return ''.join(path)

    labels = [lab(P) for P in paths]
    back = {lab(P): P for P in paths}
    src = {lab(P): st[0] for P, st in paths.items()}
    tgt = {lab(P): st[1] for P, st in paths.items()}
    idem = {v: 'e' + str(v) for v in vertices}

    def mul(a, b):
        A, B = back[a], back[b]
        if A[0] == 'e':
            return b
        if B[0] == 'e':
            return a
        C = A + B
        return lab(C) if C in paths else None

    return MultAlgebra(name, labels, src, tgt, mul, idem)


class FiniteCategory:
    """A finite category: objects, morphisms with dom/cod, composition table."""

    def __init__(self, name, objects, morphisms, dom, cod, compose, ident):
        self.name, self.objects, self.morphisms = name, objects, morphisms
        self.dom, self.cod, self.ident = dom, cod, ident
        self.comp = {}
        for f in morphisms:
            for g in morphisms:
                if dom[f] == cod[g]:
                    h = compose(f, g)          # f ∘ g
                    assert dom[h] == dom[g] and cod[h] == cod[f]
                    self.comp[(f, g)] = h
        self._check_axioms()

    def _check_axioms(self):
        for f in self.morphisms:
            assert self.comp[(f, self.ident[self.dom[f]])] == f
            assert self.comp[(self.ident[self.cod[f]], f)] == f
        for (f, g), fg in self.comp.items():
            for h in self.morphisms:
                if self.dom[g] == self.cod[h]:
                    assert self.comp[(fg, h)] == self.comp[(f, self.comp[(g, h)])], "not associative"

    def hom(self, x, y):
        return [f for f in self.morphisms if self.dom[f] == x and self.cod[f] == y]

    def algebra(self):
        """The category algebra kC: basis = morphisms, product = composition."""
        return MultAlgebra("k[" + self.name + "]", list(self.morphisms),
                           dict(self.dom), dict(self.cod),
                           lambda a, b: self.comp[(a, b)], dict(self.ident))

    def __repr__(self):
        return f"<category {self.name}: {len(self.objects)} objects, {len(self.morphisms)} morphisms>"


def group_category(name, elements, mult, e):
    """A group (or monoid) as a one-object category."""
    return FiniteCategory(name, ['*'], list(elements),
                          {g: '*' for g in elements}, {g: '*' for g in elements},
                          lambda f, g: mult(f, g), {'*': e})


def poset_category(name, elements, leq):
    """A finite poset as a category: one arrow x->y iff x <= y."""
    mors = [(x, y) for x in elements for y in elements if leq(x, y)]
    labs = {m: f"{m[0]}<{m[1]}" if m[0] != m[1] else f"1_{m[0]}" for m in mors}
    inv = {v: k for k, v in labs.items()}
    return FiniteCategory(name, list(elements), list(labs.values()),
                          {labs[m]: m[0] for m in mors}, {labs[m]: m[1] for m in mors},
                          lambda f, g: labs[(inv[g][0], inv[f][1])],
                          {x: labs[(x, x)] for x in elements})


# ---------------------------------------------------------------------------
# 4.  The (normalized, kE-relative) Hochschild cochain complex
# ---------------------------------------------------------------------------
#
#   C^n(A, A) = Hom_{E-E}( (A/E)^{⊗_E n}, A )
#
# A basis of C^n is indexed by pairs (chain, out) where
#   chain = (a1, ..., an) composable radical basis elements  a1∘a2∘...∘an,
#   out   = a basis element parallel to that composite (same source & target).
# For n = 0 the "chain" is an object x and out ∈ e_x A e_x.
#
#   (df)(a1..a_{n+1}) = a1·f(a2..a_{n+1})
#                       + Σ_i (-1)^i f(.., a_i a_{i+1}, ..)
#                       + (-1)^{n+1} f(a1..a_n)·a_{n+1}
#
# with products landing in E counting as 0 (normalization).

def _chains(A, n):
    """Composable n-tuples of radical basis elements, as (tuple, source, target)."""
    if n == 0:
        return [(('@', x), x, x) for x in A.objects]
    out = [((a,), A.src[a], A.tgt[a]) for a in A.rad]
    for _ in range(n - 1):
        out = [(ch + (b,), A.src[b], t) for ch, s, t in out for b in A.rad if A.tgt[b] == s]
    return out


def hochschild_complex(A, N, p=32003):
    """Hochschild cochains C^0..C^{N+1} of A with differentials d^0..d^N."""
    basis, chains = [], []
    for n in range(N + 2):
        ch = _chains(A, n)
        chains.append(ch)
        basis.append([(c, o) for c, s, t in ch for o in A.homs(s, t)])
    idx = [{lab: i for i, lab in enumerate(B)} for B in basis]

    def tail(ch):            # remove a1  -> chain starting at its source
        return ch[1:] if len(ch) > 1 else ('@', A.src[ch[0]])

    def head(ch):            # remove a_{n+1}
        return ch[:-1] if len(ch) > 1 else ('@', A.tgt[ch[0]])

    d = []
    for n in range(N + 1):
        M = np.zeros((len(basis[n + 1]), len(basis[n])), dtype=np.int64)
        for row, (ch, c) in enumerate(basis[n + 1]):
            a1, al = ch[0], ch[-1]
            # term a1 · f(tail)
            T = tail(ch)
            for b in A.homs(A.src[al], A.src[a1]):
                if A.mul(a1, b) == c:
                    M[row, idx[n][(T, b)]] += 1
            # middle terms
            for i in range(len(ch) - 1):
                ab = A.mul(ch[i], ch[i + 1])
                if ab is None or ab in A.E:
                    continue
                key = (ch[:i] + (ab,) + ch[i + 2:], c)
                M[row, idx[n][key]] += (-1) ** (i + 1)
            # term (-1)^{n+1} f(head) · a_{n+1}
            H = head(ch)
            for b in A.homs(A.tgt[al], A.tgt[a1]):
                if A.mul(b, al) == c:
                    M[row, idx[n][(H, b)]] += (-1) ** (n + 1)
        d.append(M % p)
    cx = CochainComplex(basis, d, p)
    cx.algebra = A
    return cx


def HH_dims(A, N, p=32003):
    cx = hochschild_complex(A, N, p)
    return [cx.betti(n) for n in range(N + 1)], cx


# ---- multiplicative structure: cup product and Gerstenhaber bracket ------

def _as_function(cx, n, vec):
    """Turn a cochain vector into  dict chain -> {out: coeff}."""
    f = {}
    for (ch, o), v in zip(cx.basis[n], vec):
        if v % cx.p:
            f.setdefault(ch, {})[o] = int(v) % cx.p
    return f


def _from_function(cx, n, f):
    v = np.zeros(len(cx.basis[n]), dtype=np.int64)
    for ch, vals in f.items():
        for o, c in vals.items():
            v[cx.index[n][(ch, o)]] = (v[cx.index[n][(ch, o)]] + c) % cx.p
    return v


def cup(cx, m, f, n, g):
    """(f ∪ g)(a1..a_{m+n}) = f(a1..a_m) · g(a_{m+1}..a_{m+n})."""
    A, p = cx.algebra, cx.p
    F, G = _as_function(cx, m, f), _as_function(cx, n, g)
    out = np.zeros(len(cx.basis[m + n]), dtype=np.int64)
    for row, (ch, c) in enumerate(cx.basis[m + n]):
        if m + n == 0:
            x = ch[1]
            left, right = ('@', x), ('@', x)
        else:
            seq = ch
            left = seq[:m] if m else ('@', A.tgt[seq[0]])
            right = seq[m:] if n else ('@', A.src[seq[-1]])
        s = 0
        for u, cu in F.get(left, {}).items():
            for w, cw in G.get(right, {}).items():
                if A.mul(u, w) == c:
                    s += cu * cw
        out[row] = s % p
    return out


def _circ_i(cx, m, F, n, G, i):
    """(f ∘_i g)(a1..a_{m+n-1}) = f(a1..a_{i-1}, g(a_i..a_{i+n-1}), ..);
    the output of g is projected to A/E (normalization)."""
    A, p = cx.algebra, cx.p
    k = m + n - 1
    out = {}
    for ch, c in cx.basis[k]:
        seq = () if k == 0 else ch
        inner = seq[i:i + n] if n else None
        if n == 0:
            # g is an element of e_x A e_x inserted at position i
            x = A.tgt[seq[i]] if i < len(seq) else (A.src[seq[-1]] if seq else ch[1])
            gvals = G.get(('@', x), {})
        else:
            gvals = G.get(inner, {})
        s = 0
        for u, cu in gvals.items():
            if u in A.E:
                continue
            new = seq[:i] + (u,) + seq[i + n:] if n else seq[:i] + (u,) + seq[i:]
            s += cu * F.get(new, {}).get(c, 0)
        if s % p:
            out.setdefault(ch, {})[c] = s % p
    return out


def gerstenhaber_bracket(cx, m, f, n, g):
    """[f, g] = f∘g - (-1)^{(m-1)(n-1)} g∘f,  f∘g = Σ_i (-1)^{(i)(n-1)} f∘_i g.
    Implemented for m, n >= 1 (degree-1 classes give the Lie algebra HH^1)."""
    assert m >= 1 and n >= 1
    p = cx.p
    F, G = _as_function(cx, m, f), _as_function(cx, n, g)

    def circ(m, F, n, G):
        tot = {}
        for i in range(m):
            sign = (-1) ** (i * (n - 1))
            for ch, vals in _circ_i(cx, m, F, n, G, i).items():
                for o, c in vals.items():
                    tot.setdefault(ch, {})
                    tot[ch][o] = (tot[ch].get(o, 0) + sign * c) % p
        return tot

    fg = _from_function(cx, m + n - 1, circ(m, F, n, G))
    gf = _from_function(cx, m + n - 1, circ(n, G, m, F))
    return (fg - (-1) ** ((m - 1) * (n - 1)) * gf) % p


def cup_table(cx, m, n):
    """Structure constants  [h_i^m] ∪ [h_j^n] = Σ_k c_ijk [h_k^{m+n}]."""
    _, Hm = cx.cohomology(m)
    _, Hn = cx.cohomology(n)
    tab = np.zeros((Hm.shape[1], Hn.shape[1], len(cx.cohomology(m + n)[1].T)), dtype=np.int64)
    for i in range(Hm.shape[1]):
        for j in range(Hn.shape[1]):
            tab[i, j] = cx.class_of(m + n, cup(cx, m, Hm[:, i], n, Hn[:, j]))
    return tab


# ---------------------------------------------------------------------------
# 5.  Baues–Wirsching cohomology of a finite category
# ---------------------------------------------------------------------------
#
# The factorization category F(C): objects = morphisms f of C, a morphism
# f -> f' is a pair (u, v) with f' = u∘f∘v.  A natural system is a functor
# D : F(C) -> k-mod; write u_* = D(u,1), v^* = D(1,v).
#
#   F^n(C; D) = Π_{x_0 <- x_1 <- ... <- x_n} D(f_1∘...∘f_n)
#   (δF)(f_1..f_{n+1}) = f_1* F(f_2..f_{n+1})
#                       + Σ_i (-1)^i F(.., f_i∘f_{i+1}, ..)
#                       + (-1)^{n+1} f_{n+1}^* F(f_1..f_n)
#
# H^n_BW(C; D) is its cohomology.  All chains are used, identities included
# (the unnormalized complex — normalized gives the same cohomology).

class NaturalSystem:
    """
    basis(f)          -> list of labels for a basis of D(f)
    push(u, f, b)     -> dict label->coeff : u_* b ∈ D(u∘f),   b ∈ D(f)
    pull(v, f, b)     -> dict label->coeff : v^* b ∈ D(f∘v)
    """

    def __init__(self, name, basis, push, pull):
        self.name, self.basis, self.push, self.pull = name, basis, push, pull


def constant_system(C):
    """D(f) = k with trivial action: BW cohomology = cohomology of the nerve BC."""
    return NaturalSystem("constant k",
                         lambda f: ['1'], lambda u, f, b: {'1': 1}, lambda v, f, b: {'1': 1})


def parallel_system(C):
    """D(f: x->y) = k·Hom_C(x, y),   u_* g = u∘g,  v^* g = g∘v.
    This is the natural system for which H_BW(C; D) ≅ HH(kC)."""
    return NaturalSystem("k Hom(x,y)",
                         lambda f: C.hom(C.dom[f], C.cod[f]),
                         lambda u, f, g: {C.comp[(u, g)]: 1},
                         lambda v, f, g: {C.comp[(g, v)]: 1})


def functor_system(C, F, side="contravariant-left"):
    """
    A k-linear representation F of C (dict object -> dim, dict morphism -> matrix)
    gives the natural system  D(f: x->y) = F(y)   ('left module' coefficients),
    u_* = F(u),  v^* = id.  Then H_BW(C; D) = Ext_{kC}(k, F) = H^*(C; F).
    """
    def basis(f):
        return [i for i in range(F['dim'][C.cod[f]])]

    def push(u, f, b):
        col = F['mat'][u][:, b]
        return {i: int(c) for i, c in enumerate(col) if c}

    def pull(v, f, b):
        return {b: 1}

    return NaturalSystem("functor", basis, push, pull)


def bw_complex(C, D, N, p=32003):
    """Baues–Wirsching cochains F^0..F^{N+1} with δ^0..δ^N."""
    def chains(n):
        if n == 0:
            return [(('@', x), x, x) for x in C.objects]
        out = [((f,), C.dom[f], C.cod[f]) for f in C.morphisms]
        for _ in range(n - 1):
            out = [(ch + (g,), C.dom[g], t) for ch, s, t in out
                   for g in C.morphisms if C.cod[g] == s]
        return out

    def composite(ch):
        if ch[0] == '@':
            return C.ident[ch[1]]
        h = ch[0]
        for g in ch[1:]:
            h = C.comp[(h, g)]
        return h

    basis = []
    for n in range(N + 2):
        basis.append([(ch, b) for ch, s, t in chains(n) for b in D.basis(composite(ch))])
    idx = [{lab: i for i, lab in enumerate(B)} for B in basis]

    d = []
    for n in range(N + 1):
        M = np.zeros((len(basis[n + 1]), len(basis[n])), dtype=np.int64)
        # one row block per (n+1)-chain; each of the n+2 terms of δ adds entries
        for ch, s, t in chains(n + 1):
            f1, fl = ch[0], ch[-1]
            T = ch[1:] if len(ch) > 1 else ('@', C.dom[f1])
            Hd = ch[:-1] if len(ch) > 1 else ('@', C.cod[fl])
            cT, cH = composite(T), composite(Hd)
            target = composite(ch)
            # f1_* F(T)
            for b in D.basis(cT):
                for c, coef in D.push(f1, cT, b).items():
                    M[idx[n + 1][(ch, c)], idx[n][(T, b)]] += coef
            # middle
            for i in range(len(ch) - 1):
                new = ch[:i] + (C.comp[(ch[i], ch[i + 1])],) + ch[i + 2:]
                for b in D.basis(target):
                    M[idx[n + 1][(ch, b)], idx[n][(new, b)]] += (-1) ** (i + 1)
            # (-1)^{n+1} f_{n+1}^* F(H)
            for b in D.basis(cH):
                for c, coef in D.pull(fl, cH, b).items():
                    M[idx[n + 1][(ch, c)], idx[n][(Hd, b)]] += (-1) ** (n + 1) * coef
        d.append(M % p)
    cx = CochainComplex(basis, d, p)
    cx.category = C
    return cx


def BW_dims(C, D, N, p=32003):
    cx = bw_complex(C, D, N, p)
    return [cx.betti(n) for n in range(N + 1)], cx
