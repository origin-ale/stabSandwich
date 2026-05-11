using Revise

using CampsPP
import PauliPropagation as pp
import CliffordMPS as cmps

using Printf

N = 13
χ = 64
thl = 1e-10
Nmax = 200

ψ = cmps.CAMPS(N)
Zs = [pp.PauliString(N, [:Z], [i], 1/N) for i = 1:N]
obs = pp.PauliSum(Zs)

output = "output/NaiveDynamics.txt"
open(output, "w") do f
    println(f, "# N=$N χ=$χ obs=magnetization")
end

ψ, k, s, evs_camps = camps_rndrotation_dynamics(ψ, χ, obs, output; showprogress=true)
open(output, "a") do f
    println(f, "# CAMPS stopped at t = $s (χ ≥ $χ)")
end

s, evs_pp = pauliprop_rndrotation_dynamics(ψ, s, thl, Nmax, obs, output; showprogress=true)
open(output, "a") do f
    println(f, "# Pauli prop. stopped at t = $s (N_pauli ≥ $Nmax)")
end

evs_tot = []
append!(evs_tot, evs_camps)
append!(evs_tot, evs_pp)
println("Expectation values:")
for ev in evs_tot
    @printf "%.4e\n" ev
end