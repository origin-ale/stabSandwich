"""Write one unified-diff patch per T10 mutant (apply with `patch -p1 -d <workdir>`).

usage: python3 -I gen_patches.py <pristine CampsPP dir> <patch out dir>
Each mutation is an exact string replacement that must match exactly once.
"""
import difflib
import os
import sys

src_root, out = sys.argv[1], sys.argv[2]
os.makedirs(out, exist_ok=True)

E = "src/EvolveQubitCAMPS.jl"
D = "src/Disentangler.jl"
P = "src/MagicTensorsPauliPropagation.jl"
Q = "src/DisentangleQubitCAMPS.jl"
T = "src/TimeEvolution.jl"

MUTANTS = [
    ("E1", E, "  return mod(m, 4)\n", "  return mod(m, 2)\n"),
    ("E3", E, "if comm(P, A) == 0x1", "if comm(P, A) == 0x0"),
    ("E4", E, "G[i] = n == 1 ? (im * P * A) : n == 2 ? -A : ((-im) * P * A)",
              "G[i] = n == 1 ? ((-im) * P * A) : n == 2 ? -A : (im * P * A)"),
    ("E5", E, "Gdag = pauli_rotation_tableau(P, mod(-n, 4))", "Gdag = pauli_rotation_tableau(P, n)"),
    ("E6", E, "Cdag.clifford = Cdag.clifford * Gdag", "Cdag.clifford = Gdag * Cdag.clifford"),
    ("E7", E, "mt.conjugate(mt.QubitPauli(P), mt.get_clifford(ψ; daggered = true))",
              "mt.conjugate(mt.QubitPauli(P), mt.get_clifford(ψ; daggered = false))"),
    ("E8", E, "mt.QubitCliffordUnitary(inv(D))", "mt.QubitCliffordUnitary(D)"),
    ("E9", E, "addmagicstate!(ψ, k, sign * ϕ)", "addmagicstate!(ψ, k, ϕ)"),
    ("E10", E, "    k += 1\n  elseif nature == :logical", "  elseif nature == :logical"),
    ("E11", E, "mt.apply_to_mps!(ψ, rotation_gate(Q, ϕ))", "mt.apply_to_mps!(ψ, rotation_gate(P, ϕ))"),
    ("E12", E, "rotation_gate(heisenberg_image(ψ, P), ϕ); cutoff = thl)",
               "rotation_gate(heisenberg_image(ψ, P), ϕ))"),
    ("E13", E, "ComplexF64[c im*s; im*s c]", "ComplexF64[c im*s; -im*s c]"),
    ("E14", E, "push!(R, mt.QubitPauli(copy(P)), im * sin(ϕ))",
               "push!(R, mt.QubitPauli(copy(P)), -im * sin(ϕ))"),
    ("D1", D, "if XQ[k+1:n] != zeros(n-k)", "if any(XQ[k+2:n])"),
    ("D3", D, "sign = (Q.phase[1] == 0x0) ? +1 : -1", "sign = +1"),
    ("D5", D, "_apply_on_sites!(D, tCX, [i, j])", "_apply_on_sites!(D, tCX, [j, i])"),
    ("P1", P, "mt.QubitPauli(pauliop_from_ppterm(term, N)), ComplexF64(coeff))",
              "mt.QubitPauli(pauliop_from_ppterm(term, N)), ComplexF64(conj(coeff)))"),
    ("Q2", Q, "    p = σ^2\n", "    p = σ\n"),
    ("Q4", Q, "  orthogonalize!(trial, n)\n", ""),
    ("T1", T, "while DisentangleCAMPS.bonddim(ψ) < χ && i < M",
              "while DisentangleCAMPS.bonddim(ψ) <= χ && i < M"),
    ("T2", T, "layer = getlayer(start_gate+1, layer_ends)", "layer = getlayer(start_gate, layer_ends)"),
    ("T6", T, "    normalize!(ψ)\n", ""),
]

for mid, rel, old, new in MUTANTS:
    path = os.path.join(src_root, rel)
    text = open(path, encoding="utf-8").read()
    n = text.count(old)
    if n != 1:
        sys.exit(f"{mid}: pattern found {n} times in {rel}")
    mutated = text.replace(old, new, 1)
    diff = difflib.unified_diff(
        text.splitlines(keepends=True), mutated.splitlines(keepends=True),
        fromfile=f"a/CampsPP/{rel}", tofile=f"b/CampsPP/{rel}", n=3)
    with open(os.path.join(out, f"{mid}.patch"), "w", encoding="utf-8") as f:
        f.writelines(diff)
    print(f"{mid}: {rel}")
print(f"{len(MUTANTS)} patches written to {out}")
