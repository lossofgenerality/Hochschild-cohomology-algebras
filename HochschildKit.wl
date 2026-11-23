(* ::Package:: *)
(* HochschildKit.wl -- Hochschild cochains, cup product, Gerstenhaber bracket,
   Hochschild cohomology and formal deformations of finite-dimensional algebras.
   (c) Computational Classnotes . CC BY-NC-SA 4.0 *)

BeginPackage["HochschildKit`"];

HochschildAlgebra::usage = "HochschildAlgebra[name, basisNames, c] builds an algebra from structure constants c[[i,j,k]] (e_i e_j = Sum_k c[[i,j,k]] e_k).";
AlgebraFromProduct::usage = "AlgebraFromProduct[name, basisNames, prod] with prod[i,j] the coefficient vector of e_i e_j.";
AlgebraFromMatrices::usage = "AlgebraFromMatrices[name, basisNames, mats] builds the algebra spanned by a list of matrices closed under product.";
DualNumbers::usage = "DualNumbers[] = k[x]/(x^2).";
TruncatedPolynomial::usage = "TruncatedPolynomial[m] = k[x]/(x^m).";
MatrixAlgebra::usage = "MatrixAlgebra[m] = M_m(k).";
UpperTriangular::usage = "UpperTriangular[m] = upper-triangular m x m matrices.";
CyclicGroupAlgebra::usage = "CyclicGroupAlgebra[m] = k[C_m].";
QuantumPlane::usage = "QuantumPlane[q] = k<x,y>/(x^2, y^2, yx - q xy); q = -1 is the exterior algebra on two generators.";
SquareZeroAlgebra::usage = "SquareZeroAlgebra[r] = k[x_1..x_r]/(x_1..x_r)^2.";

AlgebraDimension::usage = "AlgebraDimension[A].";
Mu::usage = "Mu[A] the multiplication as a 2-cochain.";
UnitCochain::usage = "UnitCochain[A] the unit 1 as a 0-cochain.";
CochainDegree::usage = "CochainDegree[f].";
ZeroCochain::usage = "ZeroCochain[A, p].";
RandomCochain::usage = "RandomCochain[A, p] random integer p-cochain.";
EvaluateCochain::usage = "EvaluateCochain[f, {v1,...,vp}] evaluates f on vectors.";

GInsert::usage = "GInsert[f, g, i] = f o_i g (g inserted into slot i of f).";
GCircle::usage = "GCircle[f, g] = Sum_i (-1)^((i-1)(q-1)) f o_i g (Gerstenhaber composition).";
GBracket::usage = "GBracket[f, g] = f o g - (-1)^((p-1)(q-1)) g o f.";
Cup::usage = "Cup[A, f, g] = mu(f(..), g(..)).";
HDelta::usage = "HDelta[A, f] Hochschild coboundary.";
DeltaMatrix::usage = "DeltaMatrix[A, p] matrix of delta : C^p -> C^(p+1) (flattened coordinates).";

HHDimension::usage = "HHDimension[A, p] = dim HH^p(A,A).";
HHDimensions::usage = "HHDimensions[A, pmax] list of dims HH^0..HH^pmax.";
HHBasis::usage = "HHBasis[A, p] cocycles representing a basis of HH^p.";
CocycleQ::usage = "CocycleQ[A, f].";
CoboundaryQ::usage = "CoboundaryQ[A, f] for a cocycle f.";
HHCoordinates::usage = "HHCoordinates[A, f] coordinates of the class of a cocycle f in HHBasis[A, p].";
HHCupTable::usage = "HHCupTable[A, p, q] class coordinates of h_i u h_j.";
HHBracketTable::usage = "HHBracketTable[A, p, q] class coordinates of [h_i, h_j].";

IdentityReport::usage = "IdentityReport[A, opts] verifies all cochain- and cohomology-level identities; returns an Association of True/False.";

PrimaryObstruction::usage = "PrimaryObstruction[A, mu1] coordinates in HH^3 of [mu1 o mu1].";
DeformAlgebra::usage = "DeformAlgebra[A, mu1, n] solves the Maurer-Cartan equation order by order.";
DeformedStructureConstants::usage = "DeformedStructureConstants[A, terms, t].";
AssociativeQ::usage = "AssociativeQ[c] for structure constants (exact or polynomial in a parameter).";
CochainForm::usage = "CochainForm[A, f] readable list of nonzero values f(e_i1,...,e_ip) -> element.";
ObstructionMap::usage = "ObstructionMap[A] the quadratic map Sq : HH^2 -> HH^3, a -> [ (a.h) o (a.h) ], as a function of symbolic coordinates.";
EquivalentDeformationsQ::usage = "EquivalentDeformationsQ[A, mu1, nu1] first-order equivalence.";

Options[HHDimension] = {Modulus -> 0};
Options[HHDimensions] = {Modulus -> 0};
Options[HHBasis] = {Modulus -> 0};

Begin["`Private`"];

(* ---------- algebras ---------- *)
HochschildAlgebra[name_, names_List, c_] := <|"Name" -> name, "Basis" -> names,
   "Structure" -> Normal[c], "Dim" -> Length[names],
   "Unit" -> unitOf[Normal[c]]|>;

unitOf[c_] := Module[{n = Length[c], u, sol, eqs},
  u = Array[uu, n];
  eqs = Flatten[{Table[Thread[u . c[[All, j, All]] == UnitVector[n, j]], {j, n}],
                 Table[Thread[u . c[[j, All, All]] == UnitVector[n, j]], {j, n}]}];
  sol = Quiet@Solve[eqs, u];
  If[sol === {}, Missing["NoUnit"], u /. First[sol]]];

AlgebraFromProduct[name_, names_, prod_] := With[{n = Length[names]},
  HochschildAlgebra[name, names, Table[prod[i, j][[k]], {i, n}, {j, n}, {k, n}]]];

AlgebraFromMatrices[name_, names_, mats_] := Module[{n = Length[mats], M},
  M = Transpose[Flatten /@ mats];
  HochschildAlgebra[name, names,
   Table[LinearSolve[M, Flatten[mats[[i]] . mats[[j]]]], {i, n}, {j, n}]]];

TruncatedPolynomial[m_] := AlgebraFromProduct["k[x]/(x^" <> ToString[m] <> ")",
  Table[If[i == 0, "1", "x^" <> ToString[i]], {i, 0, m - 1}],
  Function[{i, j}, If[i + j - 2 < m, UnitVector[m, i + j - 1], ConstantArray[0, m]]]];
DualNumbers[] := TruncatedPolynomial[2];

MatrixAlgebra[m_] := AlgebraFromMatrices["M_" <> ToString[m],
  Flatten@Table["E" <> ToString[i] <> ToString[j], {i, m}, {j, m}],
  Flatten[Table[SparseArray[{{i, j} -> 1}, {m, m}] // Normal, {i, m}, {j, m}], 1]];

UpperTriangular[m_] := AlgebraFromMatrices["T_" <> ToString[m],
  Flatten@Table["E" <> ToString[i] <> ToString[j], {i, m}, {j, i, m}],
  Flatten[Table[Normal@SparseArray[{{i, j} -> 1}, {m, m}], {i, m}, {j, i, m}], 1]];

CyclicGroupAlgebra[m_] := AlgebraFromProduct["k[C_" <> ToString[m] <> "]",
  Table["g^" <> ToString[i], {i, 0, m - 1}],
  Function[{i, j}, UnitVector[m, Mod[i + j - 2, m] + 1]]];

(* basis 1, x, y, xy ; yx = q xy *)
QuantumPlane[q_] := AlgebraFromProduct["k_q<x,y>/(x^2,y^2,yx-q xy), q=" <> ToString[q, InputForm],
  {"1", "x", "y", "xy"},
  Function[{i, j}, Which[
    i == 1, UnitVector[4, j], j == 1, UnitVector[4, i],
    i == 2 && j == 3, UnitVector[4, 4],
    i == 3 && j == 2, q UnitVector[4, 4],
    True, ConstantArray[0, 4]]]];

SquareZeroAlgebra[r_] := AlgebraFromProduct["k[x1..x" <> ToString[r] <> "]/m^2",
  Join[{"1"}, Table["x" <> ToString[i], {i, r}]],
  Function[{i, j}, Which[i == 1, UnitVector[r + 1, j], j == 1, UnitVector[r + 1, i],
    True, ConstantArray[0, r + 1]]]];

AlgebraDimension[A_] := A["Dim"];
Mu[A_] := A["Structure"];
UnitCochain[A_] := A["Unit"];

(* ---------- cochains: f of degree p is an array n x ... x n (p+1 indices),
   f[[i1,...,ip,k]] = coefficient of e_k in f(e_i1,...,e_ip) ---------- *)
CochainDegree[f_] := ArrayDepth[f] - 1;
ZeroCochain[A_, p_Integer?NonNegative] := ConstantArray[0, ConstantArray[A["Dim"], p + 1]];
RandomCochain[A_, p_, r_: 3] := RandomInteger[{-r, r}, ConstantArray[A["Dim"], p + 1]];
EvaluateCochain[f_, vs_List] := Module[{g = f}, Do[g = v . g, {v, vs}]; g];

(* the zero object 0 stands for the zero cochain of any degree, including degree -1 *)
GInsert[0, _, _] = 0; GInsert[_, 0, _] = 0; GCircle[0, _] = 0; GCircle[_, 0] = 0;
GBracket[0, _] = 0; GBracket[_, 0] = 0; Cup[_, 0, _] = 0; Cup[_, _, 0] = 0; HDelta[_, 0] = 0;
tp[a_, perm_] := If[Length[perm] <= 1, a, Transpose[a, perm]];

GInsert[f_, g_, i_Integer] := Module[{p = CochainDegree[f], q = CochainDegree[g], F, R},
  F = tp[f, Join[Range[2, i], {1}, Range[i + 1, p + 1]]];
  R = g . F;
  tp[R, Join[Range[i, i + q - 1], Range[1, i - 1], Range[i + q, p + q - 1], {p + q}]]];

zeroLike[f_, d_] := If[d < 0, 0, ConstantArray[0, ConstantArray[Last@Dimensions[f], d + 1]]];

GCircle[f_, g_] := Module[{p = CochainDegree[f], q = CochainDegree[g]},
  If[p == 0, zeroLike[f, q - 1],
   Sum[(-1)^((i - 1) (q - 1)) GInsert[f, g, i], {i, p}]]];

GBracket[f_, g_] := Module[{p = CochainDegree[f], q = CochainDegree[g]},
  GCircle[f, g] - (-1)^((p - 1) (q - 1)) GCircle[g, f]];

Cup[A_, f_, g_] := GInsert[GInsert[Mu[A], g, 2], f, 1];

HDelta[A_, f_] := Module[{mu = Mu[A], p = CochainDegree[f]},
  GInsert[mu, f, 2] + Sum[(-1)^i GInsert[f, mu, i], {i, p}] + (-1)^(p + 1) GInsert[mu, f, 1]];

DeltaMatrix[A_, p_Integer] := DeltaMatrix[A, p] = Module[{n = A["Dim"], NN, dims},
  If[p < 0, Return[SparseArray[{}, {n, 1}]]];
  NN = n^(p + 1); dims = ConstantArray[n, p + 1];
  SparseArray@Transpose@Table[Flatten@HDelta[A, ArrayReshape[UnitVector[NN, j], dims]], {j, NN}]];

(* ---------- cohomology ---------- *)
rk[m_, md_] := If[Dimensions[m][[2]] == 0, 0, MatrixRank[m, Modulus -> md]];

HHDimension[A_, p_Integer, OptionsPattern[]] := Module[{md = OptionValue[Modulus], n = A["Dim"]},
  n^(p + 1) - rk[DeltaMatrix[A, p], md] - If[p == 0, 0, rk[DeltaMatrix[A, p - 1], md]]];
HHDimensions[A_, pmax_, opts : OptionsPattern[]] := Table[HHDimension[A, p, opts], {p, 0, pmax}];

HHBasis[A_, p_Integer, OptionsPattern[]] := hhb[A, p, OptionValue[Modulus]];
hhb[A_, p_, md_] := hhb[A, p, md] = Module[{n = A["Dim"], Z, B, basis = {}, cur, r},
  Z = NullSpace[Normal@DeltaMatrix[A, p], Modulus -> md];
  B = If[p == 0, {}, Transpose[Normal@DeltaMatrix[A, p - 1]]];
  B = If[B === {}, {}, RowReduce[B, Modulus -> md]]; B = Select[B, AnyTrue[#, # != 0 &] &];
  cur = B; r = Length[cur];
  Do[If[rk[Join[cur, {z}], md] > r, AppendTo[cur, z]; AppendTo[basis, z]; r++], {z, Z}];
  ArrayReshape[#, ConstantArray[n, p + 1]] & /@ basis];

CocycleQ[A_, 0] = True;
CocycleQ[A_, f_] := AllTrue[Flatten@HDelta[A, f], # == 0 &];

CoboundaryQ[A_, 0] = True;
CoboundaryQ[A_, f_] := Module[{p = CochainDegree[f], Dm},
  If[p == 0, AllTrue[Flatten[f], # == 0 &],
   Dm = Normal@DeltaMatrix[A, p - 1];
   MatrixRank[Join[Transpose[Dm], {Flatten[f]}]] == MatrixRank[Dm]]];

HHCoordinates[A_, 0] := {};
HHCoordinates[A_, f_] := Module[{p = CochainDegree[f], H, Dm, M, sol, k},
  H = Flatten /@ HHBasis[A, p]; k = Length[H];
  If[k == 0, Return[{}]];
  Dm = If[p == 0, {}, Transpose[Normal@DeltaMatrix[A, p - 1]]];
  M = Transpose[Join[H, Dm]];
  sol = LinearSolve[M, Flatten[f]];
  Take[sol, k]];

HHCupTable[A_, p_, q_] := Outer[HHCoordinates[A, Cup[A, #1, #2]] &, HHBasis[A, p], HHBasis[A, q], 1];
HHBracketTable[A_, p_, q_] := Outer[HHCoordinates[A, GBracket[#1, #2]] &, HHBasis[A, p], HHBasis[A, q], 1];

(* ---------- identity verification ---------- *)
zeroQ[x_] := AllTrue[Flatten[{x}], # == 0 &];

IdentityReport[A_, opt : OptionsPattern[{"Degrees" -> 2, "Trials" -> 2}]] :=
 Module[{mu = Mu[A], ds = OptionValue["Degrees"], tr = OptionValue["Trials"], rc, res = <||>},
  rc[p_] := RandomCochain[A, p];
  res["mu associative (mu o mu = 0)"] = zeroQ[GCircle[mu, mu]];
  res["delta^2 = 0"] = And @@ Flatten@Table[zeroQ[HDelta[A, HDelta[A, rc[p]]]], {p, 0, ds + 1}, {tr}];
  res["delta f = (-1)^(p-1) [mu, f]"] = And @@ Flatten@Table[
      With[{f = rc[p]}, zeroQ[HDelta[A, f] - (-1)^(p - 1) GBracket[mu, f]]], {p, 0, ds + 1}, {tr}];
  res["cup associative"] = And @@ Flatten@Table[With[{f = rc[p], g = rc[q], h = rc[1]},
       zeroQ[Cup[A, Cup[A, f, g], h] - Cup[A, f, Cup[A, g, h]]]], {p, 0, ds}, {q, 0, ds}];
  res["delta(f u g) = df u g + (-1)^p f u dg"] = And @@ Flatten@Table[With[{f = rc[p], g = rc[q]},
       zeroQ[HDelta[A, Cup[A, f, g]] - Cup[A, HDelta[A, f], g] - (-1)^p Cup[A, f, HDelta[A, g]]]], {p, 0, ds}, {q, 0, ds}];
  res["pre-Lie identity"] = And @@ Flatten@Table[With[{f = rc[p], g = rc[q], h = rc[r]},
       zeroQ[GCircle[GCircle[f, g], h] - GCircle[f, GCircle[g, h]] -
         (-1)^((q - 1) (r - 1)) (GCircle[GCircle[f, h], g] - GCircle[f, GCircle[h, g]])]],
      {p, 1, ds}, {q, 0, ds}, {r, 0, ds}];
  res["bracket graded antisymmetric"] = And @@ Flatten@Table[With[{f = rc[p], g = rc[q]},
       zeroQ[GBracket[f, g] + (-1)^((p - 1) (q - 1)) GBracket[g, f]]], {p, 0, ds}, {q, 0, ds}];
  res["graded Jacobi"] = And @@ Flatten@Table[With[{f = rc[p], g = rc[q], h = rc[r]},
       zeroQ[(-1)^((p - 1) (r - 1)) GBracket[f, GBracket[g, h]] + (-1)^((q - 1) (p - 1)) GBracket[g, GBracket[h, f]] +
         (-1)^((r - 1) (q - 1)) GBracket[h, GBracket[f, g]]]], {p, 0, ds}, {q, 0, ds}, {r, 0, ds}];
  res["delta[f,g] = [f,dg] + (-1)^(q-1)[df,g]"] = And @@ Flatten@Table[With[{f = rc[p], g = rc[q]},
       zeroQ[HDelta[A, GBracket[f, g]] - GBracket[f, HDelta[A, g]] - (-1)^(q - 1) GBracket[HDelta[A, f], g]]],
      {p, 0, ds}, {q, 0, ds}];
  res["homotopy: d(f o g) = f o dg + (-1)^(q-1) df o g + (-1)^q (g u f - (-1)^(pq) f u g)"] = And @@ Flatten@Table[With[{f = rc[p], g = rc[q]},
       zeroQ[HDelta[A, GCircle[f, g]] - GCircle[f, HDelta[A, g]] - (-1)^(q - 1) GCircle[HDelta[A, f], g]
         - (-1)^q (Cup[A, g, f] - (-1)^(p q) Cup[A, f, g])]], {p, 1, ds}, {q, 1, ds}];
  (* cohomology-level statements on representatives *)
  res["HH: cup graded commutative"] = And @@ Flatten@Table[Table[
       CoboundaryQ[A, Cup[A, f, g] - (-1)^(p q) Cup[A, g, f]], {f, HHBasis[A, p]}, {g, HHBasis[A, q]}],
      {p, 0, ds}, {q, 0, ds}];
  res["HH: Poisson rule (Gerstenhaber form) [f u g, h] = [f,h] u g + (-1)^(p(r-1)) f u [g,h]"] = And @@ Flatten@Table[Table[
       CoboundaryQ[A, GBracket[Cup[A, f, g], h] - Cup[A, GBracket[f, h], g] - (-1)^(p (r - 1)) Cup[A, f, GBracket[g, h]]],
       {f, HHBasis[A, p]}, {g, HHBasis[A, q]}, {h, HHBasis[A, r]}], {p, 0, ds - 1}, {q, 0, ds - 1}, {r, 0, ds}];
  res["HH: Poisson rule (left form) [f, g u h] = (-1)^((p-1)r) [f,g] u h + g u [f,h]"] = And @@ Flatten@Table[Table[
       CoboundaryQ[A, GBracket[f, Cup[A, g, h]] - (-1)^((p - 1) r) Cup[A, GBracket[f, g], h] - Cup[A, g, GBracket[f, h]]],
       {f, HHBasis[A, p]}, {g, HHBasis[A, q]}, {h, HHBasis[A, r]}], {p, 0, ds}, {q, 0, ds - 1}, {r, 0, ds - 1}];
  res["HH: Jacobi on classes (bracket closes on HH)"] = And @@ Flatten@Table[Table[
       CocycleQ[A, GBracket[f, g]], {f, HHBasis[A, p]}, {g, HHBasis[A, q]}], {p, 0, ds}, {q, 0, ds}];
  res["HH^0 = centre"] = With[{n = A["Dim"]}, HHDimension[A, 0] ==
     Length@NullSpace[Flatten[Table[mu[[i, j, k]] - mu[[j, i, k]], {j, n}, {k, n}, {i, n}], 1]]];
  res];

(* ---------- formal deformations ---------- *)

PrimaryObstruction[A_, mu1_] := HHCoordinates[A, GCircle[mu1, mu1]];

DeformAlgebra[A_, mu1_, nmax_Integer] := Module[{terms = {mu1}, D2 = Normal@DeltaMatrix[A, 2], n = A["Dim"], Ob, sol, status = "Unobstructed"},
  If[! CocycleQ[A, mu1], Return[<|"Status" -> "mu1 is not a 2-cocycle"|>]];
  Do[
   Ob = Sum[GCircle[terms[[i]], terms[[k - i]]], {i, 1, k - 1}];
   If[! CocycleQ[A, Ob], status = "Internal: obstruction not a cocycle"; Break[]];
   sol = Quiet@Check[LinearSolve[D2, Flatten[Ob]], $Failed];
   If[sol === $Failed || Head[sol] === LinearSolve,
    status = <|"ObstructedAtOrder" -> k, "ObstructionClass" -> HHCoordinates[A, Ob]|>; Break[]];
   AppendTo[terms, ArrayReshape[sol, {n, n, n}]],
   {k, 2, nmax}];
  <|"Terms" -> terms, "Status" -> status|>];

DeformedStructureConstants[A_, terms_List, t_] := Mu[A] + Sum[t^i terms[[i]], {i, Length[terms]}];

AssociativeQ[c_] := zeroQ[Expand@GCircle[c, c]];

EquivalentDeformationsQ[A_, mu1_, nu1_] := CoboundaryQ[A, mu1 - nu1];

ObstructionMap[A_, a_Symbol: \[FormalA]] := Module[{H = HHBasis[A, 2], v},
  v = Array[a, Length[H]];
  {v, Expand@HHCoordinates[A, Expand@GCircle[v . H, v . H]]}];

CochainForm[A_, f_] := Module[{nm = A["Basis"], p = CochainDegree[f]},
  If[p == 0, Total[f*nm],
   DeleteCases[Table[With[{v = Extract[f, Append[idx, All]]},
      If[AllTrue[v, # == 0 &], Nothing, nm[[idx]] -> Total[v*nm]]],
     {idx, Tuples[Range[A["Dim"]], p]}], Nothing]]];

End[];
EndPackage[];
