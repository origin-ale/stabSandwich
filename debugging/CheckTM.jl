using Revise
using CampsPP
using DisentangleCAMPS
using ITensors, ITensorMPS
using Combinatorics

import PauliPropagation as pp
import CliffordMPS as cmps

using Printf

N = 6
t = 1
ϕ = π/4
θ = π/4

gates, phases = xxz_circuit(ϕ, θ, t, N)

start_bits = [2,4,5]
println("Starting from $start_bits")

obs = transferredmagnetization(N, start_bits)
for onebitinds in powerset(1:N, 3, 3)
  ψ = computationalcamps(N, onebitinds)
  ev = cmps.expectation(ψ, obs)
  if true
    print(onebitinds, " -> ")
    @printf "%.1f \n" real(ev)
  end
end