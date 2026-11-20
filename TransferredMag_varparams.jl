using CampsPP
using CampsPP.DisentangleCAMPS
using ITensors, ITensorMPS

using Random: MersenneTwister

@isdefined(ScriptSetup) || include(joinpath(@__DIR__, "scripts", "ScriptSetup.jl"))
using .ScriptSetup

# == Parameters ===============================================================
# Override from the command line as key=value (`help` lists them), e.g.
#   julia --project=. TransferredMag_varparams.jl N=12 μs=0.3,0.6 magic_probs=0,1
params = @parameters begin
  N = 12
  t = N ÷ 2                  # number of XXZ layers
  ϕ = π/4
  θ = π/4
  μs = [0.3, 0.6, 1., 10.]
  magic_probs = [0., 1]
  Nsamples = 100
  dope_phase = 3π/16         # doping angle, on XX-YY (`xy_magic`)
  χ = 64                     # CAMPS bond dimension that triggers PP
  thl = 1e-10                # PP truncation threshold
  Nmax_pauli = 1_000_000     # maximum number of Paulis in PP
  warn_on_prestop = true
  # BLAS and Strided threads per sample; the samples themselves run on the Julia
  # threads (-t), so more than 1 here oversubscribes the cores
  blas_threads = 1
  outdir = "output/TransferredMag_varparams"
end
nthr = Threads.nthreads()
setup_threads(blas_threads)

dope_method = "on_XY"   # no spaces: header entries are single k=v words
param_pairs = vec([(magic_prob, μ) for magic_prob in magic_probs, μ in μs])

prefix = "TMDxy"
mkpath(outdir)
output = joinpath(outdir, "$(prefix)_$(N)_$(Nsamples).txt")
output_log = joinpath(outdir, "$(prefix)_$(N)_$(Nsamples)_log.txt")
output_full = joinpath(outdir, "$(prefix)_$(N)_$(Nsamples)_full.txt")
param_info = [param_header(params; exclude = (:outdir,));
  "Δ" => string(θ/ϕ);
  "doping" => "$(round(dope_phase; digits = 3))_$dope_method";
  "n.gates" => string(length(xxz_circuit(1,1, t, N)[1]))]
obsname = "Transferred magnetization"
initialize_output(output, "$obsname (avg over $Nsamples samples)", param_info)
initialize_output(output_log, "$obsname", param_info)
initialize_output(output_full, "$obsname", param_info)

log_message("Running CAMPS-PP XXZ circuit dynamics over $(length(param_pairs)) (magic_prob, μ) \
pairs until t = $t for N=$N. χ=$χ, thl = $thl, Nmax_pauli = $Nmax_pauli, \
Nsamples = $Nsamples. Using $nthr threads. BLAS/Strided threads per sample: $blas_threads.")

# Samples run on several threads: the shared counter and the log are updated under a lock
prog = ProgressLog("sample", length(param_pairs) * Nsamples; frac = 0)
prog_lock = ReentrantLock()
ndone = Ref(0)

for pair_idx in eachindex(param_pairs)
  magic_prob, μ = param_pairs[pair_idx]
  sample_evs = Vector{Any}(undef, Nsamples)

  Threads.@threads for it in 1:Nsamples
    rng = MersenneTwister(100_000 * pair_idx + it)

    layer_ends = layerends(N, t, xxz_circuit)
    gates, phases = xxz_circuit(ϕ, θ, t, N)
    phases = xy_magic(rng, phases, magic_prob; magicphase=dope_phase)

    ψ, onebitinds = domainwallstate(rng, N, μ)
    obs = transferredmagnetization(N, onebitinds)

    evs_it, t_stop = campspp_circuit_dynamics(
      ψ, χ, thl, Nmax_pauli, gates, phases, obs;
      layer_ends = layer_ends)
    sample_evs[it] = evs_it
    lock(prog_lock) do
      if warn_on_prestop && (t_stop < length(phases))
        log_message("WARNING: sample $it for (magic_prob=$magic_prob, μ=$μ) stopped at \
gate $t_stop/$(length(phases))")
      end
      ndone[] += 1
      log_progress!(prog, ndone[], "magic_prob" => magic_prob, "μ" => μ)
    end
  end

  evs = stack_samples(sample_evs)
  save_full_samples(output_full, μ, magic_prob, evs)
  save_stats(output, evs, μ, magic_prob)
end
