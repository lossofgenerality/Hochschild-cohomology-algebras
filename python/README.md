# Python companion: Hochschild and Baues–Wirsching cohomology

A pure-Python counterpart to the Wolfram packages. It computes with exact linear algebra over 𝔽ₚ: p = 32003 stands in for characteristic 0, and small primes show modular behaviour. The only dependency is NumPy.

| File | What it is |
| --- | --- |
| `hochschild.py` | The library. It builds algebras with a multiplicative basis: quiver algebras with monomial relations, and category algebras of finite categories (groups, monoids, posets). It computes the normalized, E-relative Hochschild complex, dimensions and bases of HHⁿ, cup products, the Gerstenhaber bracket, and Baues–Wirsching cohomology of finite categories with any natural system (constant, functor, or `k·Hom` coefficients). |
| `tutorial.py` | A narrated tutorial in eight sections. Each states some theory, computes it, and checks it against the theorem (48 checks). |
| `tutorial_output.txt` | The output of `python3 tutorial.py`. |

## What the tutorial checks

1. The Hochschild complex of k[x]/x², written out by hand: HH⁰ is the centre, HH¹ the Euler derivation, HH² the deformation k[x]/(x² − t).
2. HH of k[x]/xⁿ jumps exactly when char k divides n.
3. The cup product and graded commutativity: HH\*(k[x]/x²) ≅ k[x, u, v]/(x², xu, xv, u²). In characteristic 2, odd classes square to non-zero classes.
4. The Gerstenhaber bracket: HH¹ of the Kronecker algebra is sl₂.
5. Happel's formula for hereditary path algebras.
6. Category algebras: Maschke's theorem, and the centralizer decomposition for 𝔽₂S₃ and 𝔽₃S₃.
7. Baues–Wirsching cohomology. With constant coefficients it gives the nerve (S¹, S², RP^∞). With functor coefficients it gives functor cohomology. The comparison HH\*(kC) ≅ H\*_BW(C; k·Hom) is checked on seven categories.
8. Gerstenhaber–Schack: HH\*(kP) ≅ H\*(Δ(P); k) for posets.

## Usage

```bash
cd python
python3 tutorial.py
```

```python
from hochschild import *
A = quiver_algebra("k[x]/x^2", [0], {'x': (0, 0)}, forbidden=[('x', 'x')])
HH_dims(A, 5)[0]                       # [2, 1, 1, 1, 1, 1]
C2 = group_category("C2", ['1', 'g'], lambda a, b: '1' if a == b else 'g', '1')
BW_dims(C2, parallel_system(C2), 4, p=2)[0]   # [2, 2, 2, 2, 2] = HH*(F2[C2])
```

---

© Computational Classnotes · CC BY-NC-SA 4.0
