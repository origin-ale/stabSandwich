using CampsPP
using CampsPP.DisentangleCAMPS
using ITensors, ITensorMPS

using Random: MersenneTwister

@isdefined(ScriptSetup) || include(joinpath(@__DIR__, "scripts", "ScriptSetup.jl"))
using .ScriptSetup

# == Parameters ===============================================================
# Override from the command line as key=value (`help` lists them), e.g.
#   julia --project=. TransferredMag_varparams.jl N=12 μs=0.3,0.6 magic_counts=0,10
params = @parameters begin
  N = 12
  t = N ÷ 2                  # number of XXZ layers
  ϕ = π/4
  θ = π/4
  μs = [0.3, 0.6, 1., 10.]
  # Numbers of doped XXZ triples per sample, out of (N-1)·t; each sample dopes exactly
  # n triples, drawn uniformly without replacement (default: none, and all of them)
  magic_counts = [0, (N - 1) * t]
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
ntriples = length(last(xxz_circuit(ϕ, θ, t, N))) ÷ 3
isempty(magic_counts) && error("No doping counts given")
all(n -> n isa Integer && 0 <= n <= ntriples, magic_counts) ||
  error("Doping counts must be integers in 0:$ntriples (the number of triples): $magic_counts")
# (n, μ index) pairs, n fastest
param_pairs = vec([(magic_count, j) for magic_count in magic_counts, j in eachindex(μs)])

prefix = "TMDxy"
mkpath(outdir)
output = joinpath(outdir, "$(prefix)_$(N)_$(Nsamples).txt")
output_log = joinpath(outdir, "$(prefix)_$(N)_$(Nsamples)_log.txt")
output_full = joinpath(outdir, "$(prefix)_$(N)_$(Nsamples)_full.txt")
# What each sample actually ran: its initial state and its doped triples
output_states = joinpath(outdir, "$(prefix)_$(N)_$(Nsamples)_states.txt")
output_doped = joinpath(outdir, "$(prefix)_$(N)_$(Nsamples)_doped.txt")
param_info = [param_header(params; exclude = (:outdir,));
  "Δ" => string(θ/ϕ);
  "doping" => "$(round(dope_phase; digits = 3))_$dope_method";
  "n.gates" => string(length(xxz_circuit(1,1, t, N)[1]))]
obsname = "Transferred magnetization"
initialize_output(output, "$obsname (avg over $Nsamples samples)", param_info)
initialize_output(output_log, "$obsname", param_info)
initialize_output(output_full, "$obsname", param_info)
initialize_output(output_states,
  "[sample, initial bitstring (character j is qubit j; 1 = excited), indices of excited qubits; one block per (n, μ)]",
  param_info)
initialize_output(output_doped,
  "[sample, indices of the doped XXZ triples (triple m = gates 3m-2…3m); one block per (n, μ)]",
  param_info)

log_message("Running CAMPS-PP XXZ circuit dynamics over $(length(param_pairs)) (n, μ) \
pairs until t = $t for N=$N. χ=$χ, thl = $thl, Nmax_pauli = $Nmax_pauli, \
Nsamples = $Nsamples. Using $nthr threads. BLAS/Strided threads per sample: $blas_threads.")

# Samples run on several threads: the shared counter and the log are updated under a lock
prog = ProgressLog("sample", length(param_pairs) * Nsamples; frac = 0)
prog_lock = ReentrantLock()
ndone = Ref(0)

for (magic_count, j) in param_pairs
  μ = μs[j]
  sample_evs = Vector{Any}(undef, Nsamples)
  sample_inds = Vector{Vector{Int}}(undef, Nsamples)
  sample_triples = Vector{Vector{Int}}(undef, Nsamples)

  Threads.@threads for it in 1:Nsamples
    # Separate streams: the initial state depends on (μ index, sample) only, so it is
    # the same for every n; the mask on (n, μ index, sample), so each μ gets its own
    rng_state = MersenneTwister(100_000 * j + it)
    rng_dope = MersenneTwister(hash((magic_count, j, it)))

    layer_ends = layerends(N, t, xxz_circuit)
    gates, phases = xxz_circuit(ϕ, θ, t, N)
    mask = doping_mask_count(rng_dope, phases, magic_count)
    phases = xy_magic(phases, mask; magicphase=dope_phase)

    ψ, onebitinds = domainwallstate(rng_state, N, μ)
    sample_inds[it] = onebitinds
    sample_triples[it] = findall(mask)
    obs = transferredmagnetization(N, onebitinds)

    evs_it, t_stop = campspp_circuit_dynamics(
      ψ, χ, thl, Nmax_pauli, gates, phases, obs;
      layer_ends = layer_ends)
    sample_evs[it] = evs_it
    lock(prog_lock) do
      if warn_on_prestop && (t_stop < length(phases))
        log_message("WARNING: sample $it for (n=$magic_count, μ=$μ) stopped at \
gate $t_stop/$(length(phases))")
      end
      ndone[] += 1
      log_progress!(prog, ndone[], "n" => magic_count, "μ" => μ)
    end
  end

  evs = stack_samples(sample_evs)
  save_full_samples(output_full, μ, "n" => magic_count, evs)
  save_stats(output, evs, μ, "n" => magic_count)
  open(output_states, "a") do io
    println(io, CampsPP.label_header("n" => magic_count, μ))
    for (it, inds) in enumerate(sample_inds)
      println(io, it, " ", join(q in inds ? '1' : '0' for q in 1:N), " ", join(inds, " "))
    end
    print(io, "\n\n")  # two blank lines: gnuplot index (block) separator
  end
  open(output_doped, "a") do io
    println(io, CampsPP.label_header("n" => magic_count, μ))
    for (it, ms) in enumerate(sample_triples)
      println(io, join([it; ms], " "))
    end
    print(io, "\n\n")
  end
end
