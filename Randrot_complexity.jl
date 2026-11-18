using CampsPP
using CampsPP.DisentangleCAMPS
import MagicTensors as mt

using QuantumClifford
using ITensors, ITensorMPS
using Plots
using LaTeXStrings
using Statistics
using Random

@isdefined(ScriptSetup) || include(joinpath(@__DIR__, "scripts", "ScriptSetup.jl"))
using .ScriptSetup

# For t/N = 0.1:0.1:1.1, `n_samples` random states built from t layers of (uniformly random
# N-qubit Clifford, then a rotation about a random Pauli string at angle 4π·rand), recording
# max_b S_E(b)/N of the CAMPS's MPS part and the N²-sample estimate of M₂/N.
#
#   julia --project=. Randrot_complexity.jl [Ns=…] [n_samples=…] [rng_seed=…] [fractions=…]
#   e.g. julia --project=. Randrot_complexity.jl Ns=12,24,48 n_samples=100 rng_seed=1100
#
# Reproducible: for each N, `Xoshiro(rng_seed + N)` drives the Pauli strings, the angles and
# the random Cliffords, so a run with one N reproduces that N's part of a run with several.
# The M₂ sampling (`mt.stabilizer_entropy`) draws from the global RNG and takes no RNG
# argument, so the global RNG is reseeded at the start of each N from a draw of that N's
# generator (not with the same seed: two Xoshiro streams seeded alike would coincide).

params = @parameters begin
	Ns = [12, 24, 48]
	n_samples = 100
	rng_seed = 1100
	fractions = 0.1:0.1:1.1
	outdir = "output/Randrot_complexity"
end
mkpath(outdir)

results_ee = Dict{Int, Vector{Float64}}()
results_ee_err = Dict{Int, Vector{Float64}}()
results_sre = Dict{Int, Vector{Float64}}()
results_sre_err = Dict{Int, Vector{Float64}}()

function save_results(path, Ns, results_ee, results_ee_err, results_sre, results_sre_err, fractions)
	done = filter(N -> haskey(results_ee, N), Ns)
	open(path, "w") do io
		println(io, "# rng_seed = $rng_seed, n_samples = $n_samples")
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

""" One sample: (max bond entropy / N, sampled M₂ / N) after t layers. """
function sample_complexity(rng, N, t)
	ψ = computationalcamps(N, Int[])
	paulistrings = [PauliOperator(0x0, rand(rng, Bool, N), rand(rng, Bool, N)) for _ in 1:t]
	angles = 4π * rand(rng, Float64, (t,))

	ψ_evo, _ = evolve_deepcliffords(rng, ψ, t, paulistrings, angles)

	mps = mt.get_mps(ψ_evo)
	return maximum(mt.entanglement_entropy(mps)) / N, mt.stabilizer_entropy(mps, N^2; α = 2) / N
end

for N in Ns
	rng = Xoshiro(rng_seed + N)
	Random.seed!(rand(rng, UInt64)) # global RNG, used by the M₂ sampling

	ee_vals = Float64[]
	ee_err_vals = Float64[]
	sre_vals = Float64[]
	sre_err_vals = Float64[]

	progress = ProgressLog("N=$N sample", length(fractions)*n_samples; frac = 0)
	for (i_f, f) in enumerate(fractions)
		t = round(Int, f * N)
		ee_samples = Float64[]
		sre_samples = Float64[]

		for sample in 1:n_samples
			ee, sre = sample_complexity(rng, N, t)
			push!(ee_samples, ee)
			push!(sre_samples, sre)
			log_progress!(progress, (i_f - 1) * n_samples + sample, "t/N" => f)
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

	save_results(joinpath(outdir, "gen_results.dat"), Ns, results_ee, results_ee_err, results_sre, results_sre_err, fractions)
end

xs = collect(fractions)
colors = reverse(cgrad(:viridis, length(Ns), categorical = true))

p1 = plot(xlabel = L"t/N", ylabel = L"S_E / N",
          title = L"Entanglement vs n. of $R_P (\theta)$", legend = :topleft, dpi = 300)
for (i, N) in enumerate(Ns)
	plot!(p1, xs, results_ee[N], yerror = results_ee_err[N], label = "N=$N", marker = :circle, color = colors[i])
end
savefig(p1, joinpath(outdir, "gen_avg_ee.png"))

p2 = plot(xlabel = L"t/N", ylabel = L"\mathcal{M}_2 / N",
          title = L"Magic vs n. of $R_P (\theta)$", legend = :topleft, dpi = 300)
for (i, N) in enumerate(Ns)
	plot!(p2, xs, results_sre[N], yerror = results_sre_err[N], label = "N=$N", marker = :circle, color = colors[i])
end
savefig(p2, joinpath(outdir, "gen_avg_sre.png"))

log_message("Saved plots to gen_avg_ee.png and gen_avg_sre.png, data to gen_results.dat, in $outdir")
