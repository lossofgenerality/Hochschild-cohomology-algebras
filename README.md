# Hochschild Cohomology Algebras

An exact Wolfram Language evaluator for the Hochschild cohomology of finite-dimensional algebras. It covers the Gerstenhaber bracket, the cup product, Gerstenhaber and BV algebra structures, and formal deformation theory, with worked applications from the literature. Every identity in the notebooks is checked by computation over ℚ or over 𝔽ₚ, not assumed.

## Contents

| File | What it is |
| --- | --- |
| `HochschildKit.wl` | The core package. It builds algebras from structure constants, with built-in examples. It computes the Hochschild differential, cup product, Gerstenhaber composition and bracket, the dimensions and bases of HHⁿ(A, A), and cup and bracket tables on cohomology. It includes a Maurer–Cartan deformation solver with obstruction classes in HH³, and `IdentityReport`, which verifies 15 cochain- and cohomology-level identities. |
| `HochschildApplications.wl` | Six worked applications: Poisson brackets from HH² (deformation quantization), the Clifford algebra as a quantization of the exterior algebra, group algebras and the centraliser decomposition, Morita invariance, support varieties and growth rates, and Tradler's BV operator (string topology). |
| `HochschildDemo.wl` | Runs every check and example used in the theory notebook. |
| `HochschildCohomology.nb` | The theory notebook: the Hochschild complex, the Gerstenhaber bracket and Jacobi identity, the cup product and homotopy formula, the Gerstenhaber algebra on HH, deformation theory, and dimension tables. |
| `HochschildApplications.nb` | The six applications as runnable examples. |
| `HochschildApplicationsSurvey.nb` | A literature survey of where Hochschild cohomology is applied, with clickable references, typeset theory and the worked examples. |
| `HochschildBookNotation.nb` | Hochschild algebra written the way the books write it, with infix ⊗ and ⊕ (and ∘ for Gerstenhaber composition). For new learners it includes a geometric introduction to cochains with 2D and 3D drawings, a dictionary between the notebook's objects and the books' notation, the Jacobi identity for the Gerstenhaber bracket, and motivating examples for HH⁰, HH¹ and HH². |

## Requirements

Mathematica or Wolfram Engine. The notebooks were built and tested with Mathematica 15.0.

## Usage

Keep the notebooks in the same folder as the `.wl` packages, open a notebook, and choose **Evaluation ▸ Evaluate Notebook**. Each notebook loads the packages from its own folder.

From a kernel:

```mathematica
SetDirectory["path/to/this/folder"];
Get["HochschildKit.wl"];
Get["HochschildApplications.wl"];

A = QuantumPlane[1];                 (* k[x,y]/(x^2, y^2) *)
HHDimensions[A, 3]                   (* {4, 4, 5, 6} *)
IdentityReport[TruncatedPolynomial[3]]
RunAllApplications[]
```

Each notebook carries a locked creation-date cell under the copyright line. Evaluate Notebook does not re-run it.

## Sign conventions

The bracket is [f, g] = f ∘ g − (−1)^{(p−1)(q−1)} g ∘ f, and the coboundary satisfies δf = (−1)^{p−1}[μ, f]. Two signs were fixed by exact computation:

- **Homotopy formula:** δ(f ∘ g) = f ∘ δg + (−1)^{q−1} δf ∘ g + (−1)^q (g ∪ f − (−1)^{pq} f ∪ g).
- **BV identity:** with this bracket, −[g, f] = (−1)^p Δ(f ∪ g) − (−1)^p Δf ∪ g − f ∪ Δg. This bracket is the negative transpose of the one in Benson–Kessar–Linckelmann.

---

© Computational Classnotes · CC BY-NC-SA 4.0
