using CampsPP
import PauliPropagation as pp
using Plots
using Random

N = 14
magic_prob = 1
magic_mode = :xy # Dope on XX-YY with 3π/16 or on ZZ with π/3

pp_start_layer = 3
thl = 1e-5

μ = .6
ϕ = π/4
θ = π/4

layer_ends = layerends(N, N/2, xxz_circuit)
pp_start_gate = Int(layer_ends[pp_start_layer])

if magic_mode == :xy
  magic_phase = 3π/16
  magic_doping = xy_magic
  magic_txt = "3π/16 on XX-YY"
elseif magic_mode == :z
  magic_phase = 3π/16
  magic_doping = z_magic
  magic_txt = "π/3 on ZZ"
else
  error("Unrecognised magic_mode: $magic_mode")
end

rng = MersenneTwister(1)
_, onebitinds = domainwallstate(rng, N, μ)
tm = transferredmagnetization(N, onebitinds)

gates, phases = xxz_circuit(ϕ, θ, N/2, N)
phases = magic_doping(phases, magic_prob; magicphase = magic_phase)

pp_gates = gates[pp_start_gate:end]
pp_angles = -2 .* phases[pp_start_gate:end]

@time psum = pp.propagate(pp_gates, tm, pp_angles; min_abs_coeff = thl)
ev = pp.overlapwithcomputational(psum, onebitinds)
spectrum = reverse(sort(abs.(collect(pp.coefficients(psum)))))
println("Thl $thl → $(length(spectrum)) terms → $ev")
p = bar(spectrum[1:100:end],
	xlabel = "Pauli index (sorted)", ylabel = "|coefficient|",
	title = "Pauli spectrum ($magic_txt)", yscale = :log10, legend = false)
display(p)