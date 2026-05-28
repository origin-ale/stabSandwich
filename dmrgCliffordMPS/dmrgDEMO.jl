using Pkg; Pkg.activate(".");

using ITensors
using ITensorMPS
using QuantumClifford
using Random
using CliffordMPS

# -- Model --
N = 12    # 1D system size
J = 1.0   # nearest neighbor XX coupling
hz = 0.5  # onsite Z term
hx = 0.3  # onsite X term

# -- Computational parameters --
T = 10                                # number of overall iterations
cutoff = 1E-6                         # relative SVD cutoff
seSamples = 64                        # stabilizer entropy sampling
seed = 1
mpsLinkdims = 32

# -- iTensor DMRG kwargs --
dmrg_kwargs = Dict(:nsweeps=>1, :cutoff=>cutoff)     

# -- Disentangling parameters --
deStrategy = DisentangleStrategy(:snake)
deCriterion = DisentangleCriterion(:entangle)
deMaxRepetition = 10
deMinDiff = 1E-4
deCutoff = 1E-4

# -----------------------------------------------------------------------------

Random.seed!(seed)
mps = MPS(Complex, siteinds("Qubit",N),"Y+")
#mps = random_mps(siteinds("Qubit",N);linkdims=10)
camps = CAMPS(mps)

h = PauliSum(N)
for n in 1:N
    if n < N
        push!(h, (J, PauliOperator(P"XX", N, [n,n+1])))
    end
    push!(h, (hz, PauliOperator(P"Z", N, [n])))
    push!(h, (hx, PauliOperator(P"X", N, [n])))
end

println("Hamiltonian:\n",h)

# -----------------------------------------------------------------------------

ses = fill(NaN, T)
ees = fill(NaN, T, N-1)
for t in 1:T

    dmrg!(camps, h; dmrg_kwargs...)

    diff, schedules, cost = disentangle!(
        camps, deStrategy, deMaxRepetition;
        min_diff=deMinDiff, criterion=deCriterion, cutoff=cutoff)

    ee = eEntropys!(copy(camps.mps); cutoff=cutoff)
    se = sEntropy(camps, seSamples; α=2.0)

    ees[t,:] = ee
    ses[t] = se
end

# -----------------------------------------------------------------------------

# applies local cliffords that brings the density local density matrix as close as possible to |0> 
local_rotation!(camps)

println(camps)
println("entanglement:",  ees[end,:])
println("magic:", ses[end])


# ---- expected output ----------------------------------------------------------------------------

# Activating project at `/data/git/dmrgCliffordMPS`
# [ Info: Precompiling ITensors [9136182c-28ba-11e9-034c-db9fb085ebd5]
# [ Info: Precompiling ITensorMPS [0d1a4710-d33b-49a5-8f18-73bdf49b47e2]
# [ Info: Precompiling ITensorMPSChainRulesCoreExt [0191ce28-3453-5ba4-b23a-2f3410705429]
# [ Info: Precompiling QuantumClifford [0525e862-1e90-11e9-3e4d-1b39d7109de1]
# Precompiling CliffordMPS
#   5 dependencies successfully precompiled in 115 seconds. 275 already precompiled.
# [ Info: Precompiling CliffordMPS [0b786a47-0bfe-4075-b9f4-31a5037e06a5]
# Hamiltonian:
# +0.000000e+00 + ____________
# +1.000000e+00 + XX__________
# +5.000000e-01 + Z___________
# +3.000000e-01 + X___________
# +1.000000e+00 + _XX_________
# +5.000000e-01 + _Z__________
# +3.000000e-01 + _X__________
# +1.000000e+00 + __XX________
# +5.000000e-01 + __Z_________
# +3.000000e-01 + __X_________
# +1.000000e+00 + ___XX_______
# +5.000000e-01 + ___Z________
# +3.000000e-01 + ___X________
# +1.000000e+00 + ____XX______
# +5.000000e-01 + ____Z_______
# +3.000000e-01 + ____X_______
# +1.000000e+00 + _____XX_____
# +5.000000e-01 + _____Z______
# +3.000000e-01 + _____X______
# +1.000000e+00 + ______XX____
# +5.000000e-01 + ______Z_____
# +3.000000e-01 + ______X_____
# +1.000000e+00 + _______XX___
# +5.000000e-01 + _______Z____
# +3.000000e-01 + _______X____
# +1.000000e+00 + ________XX__
# +5.000000e-01 + ________Z___
# +3.000000e-01 + ________X___
# +1.000000e+00 + _________XX_
# +5.000000e-01 + _________Z__
# +3.000000e-01 + _________X__
# +1.000000e+00 + __________XX
# +5.000000e-01 + __________Z_
# +3.000000e-01 + __________X_
# +5.000000e-01 + ___________Z
# +3.000000e-01 + ___________X

# After sweep 1 energy=-11.925930848887202  maxlinkdim=3 maxerr=2.19E-07 time=29.025
# After sweep 1 energy=-11.926081464604533  maxlinkdim=3 maxerr=7.58E-07 time=0.023
# After sweep 1 energy=-11.926081575701796  maxlinkdim=3 maxerr=7.04E-07 time=0.020
# After sweep 1 energy=-11.926081580456863  maxlinkdim=3 maxerr=7.02E-07 time=0.025
# After sweep 1 energy=-11.926081580738884  maxlinkdim=3 maxerr=7.02E-07 time=0.029
# After sweep 1 energy=-11.92608158078027  maxlinkdim=3 maxerr=7.02E-07 time=0.020
# After sweep 1 energy=-11.9260815807943  maxlinkdim=3 maxerr=7.02E-07 time=0.021
# After sweep 1 energy=-11.926081580799671  maxlinkdim=3 maxerr=7.02E-07 time=0.019
# After sweep 1 energy=-11.926081580801752  maxlinkdim=3 maxerr=7.02E-07 time=0.019
# After sweep 1 energy=-11.926081580802546  maxlinkdim=3 maxerr=7.02E-07 time=0.020
# CAMPS for 12 Qubits.
# ----- Cdag -----
# X₀₁ ⟼ - Z___________
# X₀₂ ⟼ + _Z__________
# X₀₃ ⟼ - __Z_________
# X₀₄ ⟼ + ___Z________
# X₀₅ ⟼ - ____Z_______
# X₀₆ ⟼ + _____Z______
# X₀₇ ⟼ - ______Z_____
# X₀₈ ⟼ + _______Z____
# X₀₉ ⟼ - ________Z___
# X₁₀ ⟼ + _________Z__
# X₁₁ ⟼ - __________Z_
# X₁₂ ⟼ + __________ZZ
# Z₀₁ ⟼ - X___________
# Z₀₂ ⟼ + _X__________
# Z₀₃ ⟼ - __X_________
# Z₀₄ ⟼ + ___X________
# Z₀₅ ⟼ - ____X_______
# Z₀₆ ⟼ + _____X______
# Z₀₇ ⟼ - ______X_____
# Z₀₈ ⟼ + _______X____
# Z₀₉ ⟼ - ________X___
# Z₁₀ ⟼ + _________X__
# Z₁₁ ⟼ + __________XX
# Z₁₂ ⟼ + ___________X
# ----- MPS  -----
# [2, 3, 2, 2, 2, 2, 3, 3, 3, 3, 2]

# entanglement:[0.029874764526721005, 0.009343106045279288, 0.006052461065675604, 0.0047275787958646846, 0.004583257182464007, 0.0046612764756923586, 0.004924218254430998, 0.006620662889742411, 0.010043738823905175, 0.029690519528959545, 0.043609890692081044]
# magic:1.4660094352371953
