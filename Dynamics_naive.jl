using Revise

using CampsPP
import PauliPropagation as pp
import CliffordMPS as cmps

N = 13
χ = 64
thl = 1e-10
Nmax = 100

ψ = cmps.CAMPS(N)
Zs = [pp.PauliString(N, [:Z], [i], 1/N) for i = 1:N]
obs = pp.PauliSum(Zs)

output = "output/NaiveDynamics.txt"
open(output, "w") do f
    println(f, "# N=$N χ=$χ obs=magnetization")
end

ψ, k, s = camps_rndrotation_dynamics(ψ, χ, obs, output; showprogress=true)
open(output, "a") do f
    println(f, "# CAMPS stopped at t = $s (χ ≥ $χ)")
end

s = pauliprop_rndrotation_dynamics(ψ, s, thl, Nmax, obs, output; showprogress=true)
open(output, "a") do f
    println(f, "# Pauli prop. stopped at t = $s (N_pauli ≥ $Nmax)")
end
