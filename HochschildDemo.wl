(* HochschildDemo.wl -- runs every check and example used in the write-up.
   Usage: Get["HochschildKit.wl"]; Get["HochschildDemo.wl"]
   (c) Computational Classnotes . CC BY-NC-SA 4.0 *)

Needs["HochschildKit`"];
SeedRandom[2026];

algebras = <|
  "k[x]/(x^2)" -> DualNumbers[], "k[x]/(x^3)" -> TruncatedPolynomial[3],
  "M_2(k)" -> MatrixAlgebra[2], "T_2(k)" -> UpperTriangular[2],
  "k[C_2]" -> CyclicGroupAlgebra[2], "k[C_3]" -> CyclicGroupAlgebra[3],
  "k[x,y]/(x,y)^2" -> SquareZeroAlgebra[2],
  "k[x,y]/(x^2,y^2)" -> QuantumPlane[1], "Lambda(k^2)" -> QuantumPlane[-1],
  "quantum plane q=2" -> QuantumPlane[2]|>;

(* 1. Hochschild cohomology dimensions *)
dims = Association@KeyValueMap[#1 -> HHDimensions[#2, If[AlgebraDimension[#2] > 3, 3, 4]] &, algebras];
charP = <|"k[x]/(x^2), char 2" -> HHDimensions[DualNumbers[], 4, Modulus -> 2],
          "k[x]/(x^3), char 3" -> HHDimensions[TruncatedPolynomial[3], 4, Modulus -> 3],
          "k[C_2], char 2" -> HHDimensions[CyclicGroupAlgebra[2], 4, Modulus -> 2]|>;

(* 2. every identity, cochain level and cohomology level *)
reports = IdentityReport[#, "Degrees" -> 2] & /@ KeyTake[algebras, {"k[x]/(x^3)", "T_2(k)", "Lambda(k^2)", "k[x,y]/(x^2,y^2)"}];
allPassed = AllTrue[Flatten[Values /@ Values[reports]], TrueQ];

(* 3. Gerstenhaber algebra of k[x]/(x^2) *)
D1 = DualNumbers[];
dualReps = Table[CochainForm[D1, #] & /@ HHBasis[D1, p], {p, 0, 4}];
dualCup = Table[HHCupTable[D1, p, q], {p, 1, 3}, {q, 1, 3}];
dualBracket = Table[HHBracketTable[D1, p, q], {p, 1, 3}, {q, 1, 3}];

(* 4. deformations *)
dualDeformation = DeformAlgebra[D1, First@HHBasis[D1, 2], 5];     (* x^2 = t . 1 , exact at order 1 *)
Q1 = QuantumPlane[1];
qDirection = Normal@SparseArray[{{3, 2, 4} -> 1}, {4, 4, 4}];   (* yx = (1+t) xy *)
qDeformation = DeformAlgebra[Q1, qDirection, 5];
qAssociative = AssociativeQ[DeformedStructureConstants[Q1, qDeformation["Terms"], t]];
{sqVars, sq} = ObstructionMap[Q1];
H2 = HHBasis[Q1, 2];
obstructed = DeformAlgebra[Q1, H2[[1]] + H2[[4]], 4];
rigidM2 = HHDimension[MatrixAlgebra[2], 2] == 0;

Print[Column[{dims, charP, allPassed, dualCup, dualBracket, dualDeformation["Status"],
  qDeformation["Status"], qAssociative, sq, obstructed["Status"], rigidM2}]];
