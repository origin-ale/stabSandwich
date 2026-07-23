using DisentangleCAMPS
using CliffordMPS

using QuantumClifford
using ITensors, ITensorMPS
using ProgressMeter
using Plots
using LaTeXStrings
using Statistics

Ns = length(ARGS) >= 1 ? parse.(Int, split(ARGS[1], ",")) : [12, 24, 48]
n_samples = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 100
fractions = 0.1:0.1:1.1

results_ee = Dict{Int, Vector{Float64}}()
results_ee_err = Dict{Int, Vector{Float64}}()
results_sre = Dict{Int, Vector{Float64}}()
results_sre_err = Dict{Int, Vector{Float64}}()

generate_showvalues(toN, sample) = () -> [("t/N",toN), ("sample", sample)]

function save_results(path, Ns, results_ee, results_ee_err, results_sre, results_sre_err, fractions)
	done = filter(N -> haskey(results_ee, N), Ns)
	open(path, "w") do io
		for (i, N) in enumerate(done)
			println(io, "# N = $N")
			println(io, "# t/N avg_EE err_EE avg_SRE err_SRE")
			for (j, f) in enumerate(fractions)
				println(io, "$f $(results_ee[N][j]) $(results_ee_err[N][j]) $(results_sre[N][j]) $(results_sre_err[N][j])")
			end
			i < length(done) && print(io, "\n\n")
		end
	end
end

for N in Ns
	ee_vals = Float64[]
	ee_err_vals = Float64[]
	sre_vals = Float64[]
	sre_err_vals = Float64[]

	progressbar = Progress(length(fractions)*n_samples; desc = "N=$N")
	for f in fractions
		t = round(Int, f * N)
		ee_samples = Float64[]
		sre_samples = Float64[]

		for sample in 1:n_samples
			ψ = CAMPS(N)
      xbits = fill(false, N)
      ybits = fill(false, N)
      ybits[1] = true
      
      paulistrings = fill(PauliOperator(0x0, xbits, ybits),t)
      phases = fill(π/8, t)

			ψ_evo, k = evolve(ψ, t, paulistrings, phases)

			push!(ee_samples, maximum(eEntropys!(ψ_evo.mps)) / N)
			push!(sre_samples, sEntropy(ψ_evo.mps, N^2; α =2) / N)
			next!(progressbar, showvalues=generate_showvalues(f, sample))
		end

		push!(ee_vals, mean(ee_samples))
		push!(ee_err_vals, std(ee_samples) / sqrt(n_samples))
		push!(sre_vals, mean(sre_samples))
		push!(sre_err_vals, std(sre_samples) / sqrt(n_samples))
	end

	results_ee[N] = ee_vals
	results_ee_err[N] = ee_err_vals
	results_sre[N] = sre_vals
	results_sre_err[N] = sre_err_vals

	save_results("output/gen_results.dat", Ns, results_ee, results_ee_err, results_sre, results_sre_err, fractions)
end

xs = collect(fractions)
colors = reverse(cgrad(:viridis, length(Ns), categorical = true))

p1 = plot(xlabel = L"t/N", ylabel = L"S_E / N",
          title = L"Entanglement vs n. of $R_P (\theta)$", legend = :topleft, dpi = 300)
for (i, N) in enumerate(Ns)
	plot!(p1, xs, results_ee[N], yerror = results_ee_err[N], label = "N=$N", marker = :circle, color = colors[i])
end
savefig(p1, "output/gen_avg_ee.png")

p2 = plot(xlabel = L"t/N", ylabel = L"\mathcal{M}_2 / N",
          title = L"Magic vs n. of $R_P (\theta)$", legend = :topleft, dpi = 300)
for (i, N) in enumerate(Ns)
	plot!(p2, xs, results_sre[N], yerror = results_sre_err[N], label = "N=$N", marker = :circle, color = colors[i])
end
savefig(p2, "output/gen_avg_sre.png")

println("Saved plots to gen_avg_ee.png and gen_avg_sre.png, data to gen_results.dat")