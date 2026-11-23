(* ::Package:: *)
(* HochschildApplications.wl -- worked applications of Hochschild cohomology,
   built on HochschildKit.wl.
   Usage:  Get["HochschildKit.wl"]; Get["HochschildApplications.wl"]; RunAllApplications[]
   (c) Computational Classnotes . CC BY-NC-SA 4.0 *)

BeginPackage["HochschildApplications`", {"HochschildKit`"}];

SymmetricPart::usage = "SymmetricPart[phi] = (phi(a,b) + phi(b,a))/2 for a 2-cochain.";
AntisymmetricPart::usage = "AntisymmetricPart[phi] = (phi(a,b) - phi(b,a))/2 for a 2-cochain.";
BiderivationQ::usage = "BiderivationQ[A, b] tests b(ac,d) = a b(c,d) + b(a,d) c and b(a,cd) = c b(a,d) + b(a,c) d.";
PoissonQ::usage = "PoissonQ[A, b] tests that an antisymmetric biderivation b satisfies the Jacobi identity.";
PoissonFromDeformation::usage = "PoissonFromDeformation[A] splits each HH^2 basis class of a commutative algebra into its symmetric (Harrison) and antisymmetric (Poisson) parts.";
CliffordFamily::usage = "CliffordFamily[t] = k<x,y>/(x^2 - t, y^2 - t, xy + yx), basis 1, x, y, xy.";
CliffordTerms::usage = "CliffordTerms[] = {mu0, mu1, mu2} with mu_t = mu0 + t mu1 + t^2 mu2.";
GroupAlgebra::usage = "GroupAlgebra[name, G] the group algebra of a finite permutation group G.";
CentralizerPrediction::usage = "CentralizerPrediction[G, h, n] = Sum over conjugacy classes [g] of h[C_G(g), n], where h[H, n] = dim H^n(H, k).";
MatrixRingOver::usage = "MatrixRingOver[A, m] = M_m(A) from the structure constants of A.";
GroupCohomologyDim::usage = "GroupCohomologyDim[p][C, n] = dim H^n(C, F_p) for the centralisers occurring in the examples (cyclic groups, and S3 at p = 3); p = 0 means characteristic 0.";
BVIdentityQ::usage = "BVIdentityQ[A, eps, pmax] verifies -[g,f] = (-1)^p Delta(f u g) - (-1)^p Delta(f) u g - f u Delta(g) on all pairs of basis classes, and that delta Delta = -Delta delta.";
KleinFourAlgebra::usage = "KleinFourAlgebra[] = k[C2 x C2].";
GrowthRate::usage = "GrowthRate[dims] estimates r with dim HH^n ~ c n^(r-1) (the rate of growth, i.e. complexity) from a list of dimensions.";
TraceForm::usage = "TraceForm[A, eps] the bilinear form <a,b> = eps(ab).";
BVOperator::usage = "BVOperator[A, eps, f] Tradler's BV operator Delta : C^n -> C^(n-1) for a symmetric Frobenius algebra with trace eps.";
RunAllApplications::usage = "RunAllApplications[] runs the six worked examples and returns an Association of results.";

Begin["`Private`"];

(* ---------- helpers: evaluate a 2-cochain or the product on vectors ---------- *)
bil[f_][u_, v_] := v . (u . f);
basis[A_] := IdentityMatrix[AlgebraDimension[A]];

(* ---------- 1. Poisson brackets from HH^2 (HKR in degree 2) ---------- *)
SymmetricPart[phi_] := (phi + Transpose[phi, {2, 1, 3}])/2;
AntisymmetricPart[phi_] := (phi - Transpose[phi, {2, 1, 3}])/2;

BiderivationQ[A_, b_] := Module[{m = bil[Mu[A]], B = bil[b], Eb = basis[A]},
  AllTrue[Flatten@Table[{
      B[m[a, c], d] - (m[a, B[c, d]] + m[B[a, d], c]),
      B[a, m[c, d]] - (m[c, B[a, d]] + m[B[a, c], d])},
     {a, Eb}, {c, Eb}, {d, Eb}], # == 0 &]];

PoissonQ[A_, b_] := Module[{B = bil[b], Eb = basis[A]},
  BiderivationQ[A, b] && AllTrue[Flatten@Table[
     B[a, B[c, d]] + B[c, B[d, a]] + B[d, B[a, c]], {a, Eb}, {c, Eb}, {d, Eb}], # == 0 &]];

PoissonFromDeformation[A_] := Module[{H = HHBasis[A, 2]},
  Table[With[{s = SymmetricPart[h], w = AntisymmetricPart[h]},
    <|"Class" -> CochainForm[A, h],
      "Symmetric part is a cocycle" -> CocycleQ[A, s],
      "Antisymmetric part is a cocycle" -> CocycleQ[A, w],
      "Antisymmetric part is a Poisson bracket" -> PoissonQ[A, w],
      "Bracket" -> CochainForm[A, w],
      "Symmetric part is a coboundary" -> CoboundaryQ[A, s]|>], {h, H}]];

(* ---------- 2. Clifford algebra as a quantization of the exterior algebra ---------- *)
(* basis 1, x, y, xy with x^2 = y^2 = t, yx = -xy *)
cliffordProduct[t_][i_, j_] := Module[{w = {{}, {1}, {2}, {1, 2}}, word, sign = 1, coef = 1, k},
  word = Join[w[[i]], w[[j]]];
  (* bubble sort with sign, then square generators to t *)
  Do[Do[If[word[[k]] > word[[k + 1]], word[[{k, k + 1}]] = word[[{k + 1, k}]]; sign = -sign],
     {k, 1, Length[word] - 1}], {Length[word]}];
  While[Length[word] >= 2 && MemberQ[Partition[word, 2, 1], {a_, a_}],
   k = First@FirstPosition[Partition[word, 2, 1], {a_, a_}];
   word = Delete[word, {{k}, {k + 1}}]; coef *= t];
  sign coef UnitVector[4, First@FirstPosition[w, word]]];

CliffordFamily[t_] := AlgebraFromProduct["Cl_t", {"1", "x", "y", "xy"}, cliffordProduct[t]];

CliffordTerms[] := Module[{c = Mu[CliffordFamily[\[FormalT]]]},
  Table[Coefficient[c, \[FormalT], k], {k, 0, 2}] /. \[FormalT] -> 0];

(* ---------- 3. group algebras and the centraliser decomposition ---------- *)
GroupAlgebra[name_, G_] := Module[{els = GroupElements[G], idx},
  idx = AssociationThread[els -> Range[Length[els]]];
  AlgebraFromProduct[name, ToString /@ els,
   Function[{i, j}, UnitVector[Length[els], idx[PermutationProduct[els[[i]], els[[j]]]]]]]];

conjugacyClassReps[G_] := Module[{els = GroupElements[G]},
  First /@ DeleteDuplicates[
    Table[Union@Table[PermutationProduct[InversePermutation[h], g, h], {h, els}], {g, els}]]];
centralizer[G_, g_] := Select[GroupElements[G], PermutationProduct[#, g] === PermutationProduct[g, #] &];

(* h[order, centraliser element list, n] = dim H^n(C, k) supplied by the caller *)
CentralizerPrediction[G_, h_, n_] := Total[h[centralizer[G, #], n] & /@ conjugacyClassReps[G]];

(* ---------- 4. Morita invariance ---------- *)
MatrixRingOver[A_, m_] := Module[{d = AlgebraDimension[A], c = Mu[A], idx, n, prod},
  idx = Flatten[Table[{i, j, a}, {i, m}, {j, m}, {a, d}], 2]; n = Length[idx];
  prod[p_, q_] := Module[{v = ConstantArray[0, n]},
    If[idx[[p, 2]] == idx[[q, 1]],
     Do[v[[Position[idx, {idx[[p, 1]], idx[[q, 2]], b}][[1, 1]]]] += c[[idx[[p, 3]], idx[[q, 3]], b]], {b, d}]];
    v];
  AlgebraFromProduct["M_" <> ToString[m] <> "(" <> A["Name"] <> ")",
   ("E" <> ToString[#[[1]]] <> ToString[#[[2]]] <> "." <> A["Basis"][[#[[3]]]]) & /@ idx, prod]];

(* ---------- 5. complexity from growth ---------- *)
GrowthRate[dims_List] := Module[{d = Rest[dims], k = 0},
  (* smallest k with k-th differences zero: dim HH^n is a polynomial of degree k-1 in n, so the rate of growth is k *)
  If[Length[d] < 2, Return[Missing["TooFewDegrees"]]];
  While[k < Length[d] && ! AllTrue[Differences[d, k], # == 0 &], k++];
  k];

(* dim H^n(C, F_p) for C cyclic, or S3 at p = 3; enough for the S3 and abelian examples *)
GroupCohomologyDim[p_][C_List, n_] := With[{o = Length[C]},
  Which[n == 0, 1, p == 0 || Mod[o, p] != 0, 0, o == 6 && p == 3, If[MemberQ[{0, 3}, Mod[n, 4]], 1, 0], True, 1]];

(* ---------- 6. BV operator on a symmetric Frobenius algebra (Tradler) ---------- *)
TraceForm[A_, eps_] := Table[eps . bil[Mu[A]][a, b], {a, basis[A]}, {b, basis[A]}];

BVOperator[A_, eps_, f_] := Module[{n = CochainDegree[f], Fe, R, G = TraceForm[A, eps], rot},
  If[n <= 0, Return[0]];
  Fe = f . eps;                                   (* Fe(a1..an) = eps(f(a1..an)) *)
  (* term i: Fe(a_i,...,a_n,a_1,...,a_{i-1}); as a function of (a1..an) *)
  rot[i_] := If[n == 1, Fe, Transpose[Fe, RotateLeft[Range[n], i - 1]]];
  R = Sum[(-1)^(i (n - 1)) rot[i], {i, n}];       (* R[a1..a_{n-1}, a_n] = <Delta f(a1..), a_n> *)
  R . Inverse[G]];

BVIdentityQ[A_, eps_, pmax_] := Module[{DD = BVOperator[A, eps, #] &},
  <|"Symmetric trace form" -> (TraceForm[A, eps] === Transpose[TraceForm[A, eps]]),
    "delta Delta = -Delta delta" -> And @@ Table[With[{f = RandomCochain[A, p]},
        AllTrue[Flatten[{HDelta[A, DD[f]] + DD[HDelta[A, f]]}], # == 0 &]], {p, 1, pmax + 1}],
    "BV identity on HH" -> And @@ Flatten@Table[Table[
        CoboundaryQ[A, GBracket[g, f] + ((-1)^p DD[Cup[A, f, g]] - (-1)^p Cup[A, DD[f], g] - Cup[A, f, DD[g]])],
        {f, HHBasis[A, p]}, {g, HHBasis[A, q]}], {p, 0, pmax}, {q, 0, pmax}]|>];

KleinFourAlgebra[] := AlgebraFromProduct["k[C2xC2]", {"e", "a", "b", "ab"},
  Function[{i, j}, UnitVector[4, BitXor[i - 1, j - 1] + 1]]];

(* ---------- runner: the six worked examples ---------- *)
RunAllApplications[] := Module[{res = <||>, S3 = GroupAlgebra["k[S3]", SymmetricGroup[3]], L = QuantumPlane[-1], t},
  (* 1 *)
  res["1 Poisson brackets from HH^2"] = PoissonFromDeformation[QuantumPlane[1]];
  (* 2 *)
  t = CliffordTerms[];
  res["2 Clifford quantization"] = <|"mu0 = exterior product" -> (t[[1]] == Mu[L]),
    "mu1 cocycle" -> CocycleQ[L, t[[2]]], "mu1 class nonzero" -> ! CoboundaryQ[L, t[[2]]],
    "mu1" -> CochainForm[L, t[[2]]], "mu2" -> CochainForm[L, t[[3]]],
    "family associative for all t" -> AssociativeQ[Mu[CliffordFamily[\[FormalT]]]],
    "HH(Lambda)" -> HHDimensions[L, 3], "HH(Cl, t=1)" -> HHDimensions[CliffordFamily[1], 3],
    "HH(Cl, t=-1)" -> HHDimensions[CliffordFamily[-1], 3]|>;
  (* 3 *)
  res["3 group algebra k[S3]"] = Association@Table[
     ("char " <> ToString[p]) -> <|"computed" -> HHDimensions[S3, 3, Modulus -> p],
       "centraliser prediction" -> Table[CentralizerPrediction[SymmetricGroup[3], GroupCohomologyDim[p], n], {n, 0, 3}]|>,
     {p, {0, 2, 3}}];
  (* 4 *)
  res["4 Morita invariance"] = <|"HH(k[x]/x^2)" -> HHDimensions[DualNumbers[], 2],
     "HH(M2(k[x]/x^2))" -> HHDimensions[MatrixRingOver[DualNumbers[], 2], 2],
     "HH(M2(k))" -> HHDimensions[MatrixAlgebra[2], 2]|>;
  (* 5 *)
  res["5 complexity"] = Association@Map[#[[1]] -> <|"dims" -> #[[2]], "rate of growth" -> GrowthRate[#[[2]]]|> &,
     {{"k[C2], char 2", HHDimensions[CyclicGroupAlgebra[2], 5, Modulus -> 2]},
      {"k[C2xC2], char 2", HHDimensions[KleinFourAlgebra[], 4, Modulus -> 2]},
      {"k[C2xC2], char 0", HHDimensions[KleinFourAlgebra[], 3]},
      {"Lambda(k^2), char 0", HHDimensions[L, 4]},
      {"k[x]/x^2, char 0", HHDimensions[DualNumbers[], 5]}}];
  (* 6 *)
  res["6 BV structure"] = <|"k[x]/x^2" -> BVIdentityQ[DualNumbers[], {0, 1}, 3],
     "k[x]/x^3" -> BVIdentityQ[TruncatedPolynomial[3], {0, 0, 1}, 2],
     "k[x,y]/(x^2,y^2)" -> BVIdentityQ[QuantumPlane[1], {0, 0, 0, 1}, 2]|>;
  res];

End[];
EndPackage[];
