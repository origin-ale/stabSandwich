import CliffordMPS as cmps
import PauliPropagation as pp
using QuantumClifford
using DisentangleCAMPS
using ProgressMeter

using Plots
using LaTeXStrings

include("./CiffordMPSPauliPropagation.jl")

Ns = length(ARGS) >= 1 ? parse.(Int, split(ARGS[1], ",")) : [12, 16, 24]
n_samples = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 50
fractions = 0.1:0.1:1.1

results_ee = Dict{Int, Vector{Float64}}()
results_sre = Dict{Int, Vector{Float64}}()

generate_showvalues(toN, sample) = () -> [("t/N",toN), ("sample", sample)]

for N in Ns
	ee_vals = Float64[]
	sre_vals = Float64[]

	progressbar = Progress(length(fractions)*n_samples; desc = "N=$N")
	for f in fractions
		t = round(Int, f * N)
		ee_sum = 0.0
		sre_sum = 0.0

		for sample in 1:n_samples
			ψ = cmps.CAMPS(N)
      paulistrings = [pp.PauliString(N, rand(0:4^N-1), rand([1,im,-1,-im])) for _ in 1:t]
      phases = 2π*rand(Float64, (t,))

			ψ_evo, k = evolve(ψ, t, paulistrings, phases)

			ee_sum += maximum(cmps.eEntropys!(ψ_evo.mps)) / N
			sre_sum += cmps.sEntropy(ψ_evo.mps, N^2; α =2) / N
			next!(progressbar, showvalues=generate_showvalues(f, sample))
		end
		avg_ee = ee_sum / n_samples
		avg_sre = sre_sum / n_samples

		push!(ee_vals, avg_ee)
		push!(sre_vals, avg_sre)
	end

	results_ee[N] = ee_vals
	results_sre[N] = sre_vals
end

open("output/gen_results.dat", "w") do io
	for (i, N) in enumerate(Ns)
		println(io, "# N = $N")
		println(io, "# t/N avg_EE avg_SRE")
		for (j, f) in enumerate(fractions)
			println(io, "$f $(results_ee[N][j]) $(results_sre[N][j])")
		end
		i < length(Ns) && print(io, "\n\n")
	end
end

xs = collect(fractions)
colors = reverse(cgrad(:viridis, length(Ns), categorical = true))

p1 = plot(xlabel = L"t/N", ylabel = L"S_E / N",
          title = L"Entanglement vs n. of $R_P (\theta)$", legend = :topleft, dpi = 300)
for (i, N) in enumerate(Ns)
	plot!(p1, xs, results_ee[N], label = "N=$N", marker = :circle, color = colors[i])
end
savefig(p1, "output/gen_avg_ee.png")

p2 = plot(xlabel = L"t/N", ylabel = L"\mathcal{M}_2 / N",
          title = L"Magic vs n. of $R_P (\theta)$", legend = :topleft, dpi = 300)
for (i, N) in enumerate(Ns)
	plot!(p2, xs, results_sre[N], label = "N=$N", marker = :circle, color = colors[i])
end
savefig(p2, "output/gen_avg_sre.png")

println("Saved plots to gen_avg_ee.png and gen_avg_sre.png, data to gen_results.dat")
