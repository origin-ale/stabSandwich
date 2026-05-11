using CampsPP

import PauliPropagation as pp
import CliffordMPS as cmps

using Printf

N = 13
t = 30
χ = 64
thl = 1e-10
Nmax = 200

ψ = cmps.CAMPS(N)
Zs = [pp.PauliString(N, [:Z], [i], 1/N) for i = 1:N]
obs = pp.PauliSum(Zs)

gates, phases = rotation_circuit(t, N)

output = "output/CircuitDynamics.txt"
open(output, "w") do f
    println(f, "# N=$N χ=$χ obs=magnetization")
end

ψ_evo, _, s, evs_camps = camps_circuit_dynamics(ψ, χ, gates, phases, obs, output; showprogress = true)
leftover_gates, leftover_angles = leftover_rotgates(s, gates, phases)
s, evs_pp = pauliprop_circuit_dynamics(ψ_evo, s, thl, leftover_gates, leftover_angles, Nmax, obs, output; showprogress=true)
open(output, "a") do f
    println(f, "# Pauli prop. stopped at t = $s (N_pauli ≥ $Nmax)")
end
