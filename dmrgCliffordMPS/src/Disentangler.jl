using CliffordMPS
using ITensors
using Base.Threads

# -- disengangle criterion ------------------------------------------------------------------------

"""Tag type selecting the cost function used to rank disentangling candidates."""
struct DisentangleCriterion{x}
end
DisentangleCriterion(x) = DisentangleCriterion{x}()

function improvement_value(::DisentangleCriterion{:entangle}, specs::Array)
    return float(sum(eEntropy.([spec.^2 for spec in specs])))
end

function improvement_value(::DisentangleCriterion{:chi3}, specs::Array)
    return float(sum([length(spec)^3 for spec in specs]))
end

# -- disengangle stategy --------------------------------------------------------------------------

abstract type DisentangleStrategys end

"""Tag type selecting the scan pattern used when searching for disentangling gates."""
struct DisentangleStrategy{x} <: DisentangleStrategys
end

DisentangleStrategy(x) = DisentangleStrategy{x}()

# -- disentangler ---------------------------------------------------------------------------------

"""
Search for the best two-site Clifford gate that reduces entanglement between `nL` and `nR`.

Returns the selected Clifford index number, the improvement in the chosen cost
function, and the rough SVD cost of the search.
"""
function disentangler(mps::MPS, nL::Int, nR::Int;
    index_set=Clifford2IndexSet(:entangle), criterion=DisentangleCriterion(:entangle), svd_kwargs...)
    @assert nL < nR
    mpsPrep = copy(mps)

    orthogonalize!(mpsPrep, nL)
    applySwaps!(mpsPrep, nL, nR-1; svd_kwargs...)

    id_value = nothing
    values = fill(Inf,length(index_set))

    costs = nR-nL

    # Threads.@threads for i in 0:length(index_set) # This does not make things faster :/ 
    for i in 0:length(index_set)
            tmps = copy(mpsPrep)
        if i == 0
            cliffordIndex = Clifford2Index(0)
        else
            cliffordIndex = Clifford2Index(index_set, i-1)
        end
        specs = []
        push!(specs, applyGate!(tmps, cliffordIndex, nR-1, nR; leftOC=true))
        append!(specs, applySwaps!(tmps, nR-1, nL; svd_kwargs...))
        costs += nR-nL

        value = improvement_value(criterion, specs)
        if i == 0
            id_value = value
        else
            values[i] = value
        end
    end

    min_value, min_num = findmin(values)
    if min_value < id_value
        return clifford_number(Clifford2Index(index_set,min_num-1)), id_value - min_value, costs
    else
        return 0, 0.0, costs
    end
end

"""Search all masked pairs of sites and return the best Clifford move at each bond."""
function disentangler(mps::MPS, mask::AbstractArray{Bool,2};
    index_set=Clifford2IndexSet(:entangle), criterion=DisentangleCriterion(:entangle), svd_kwargs...)
    @assert size(mask,1) == size(mask,2)
    N = size(mask,1)

    best_clif_nums=fill(0::Int,N,N)
    differences=fill(-Inf64,N,N)

    costs = 0

    for ind in CartesianIndices(mask)
        if mask[ind] && ind[1] < ind[2]
            best_clif_nums[ind], differences[ind], cost = disentangler(mps, ind[1], ind[2];
                index_set=index_set, criterion=criterion, svd_kwargs...)
            costs += cost
        end
    end
    return best_clif_nums, differences, costs
end

"""Convenience wrapper that searches all pairs of sites in the MPS."""
function disentangler(mps::MPS; kwargs...)
    N = length(mps)
    mask = fill(true,N,N)
    return disentangler(mps, mask; kwargs...)
end


# -- Disentangle CAMPS ------------------------------------------------------------------------------

function _make_full_schedule(best_clif_nums::AbstractArray{Int, 2}, differences::AbstractArray{Float64, 2})
    total_difference = 0.0
    schedule = []
    while true
        diff, ind = findmax(differences)
        if !(diff>0.0)
            break
        end
        differences[1:ind[2],ind[1]:end] .= -Inf64
        total_difference += diff
        push!(schedule, (best_clif_nums[ind], ind[1], ind[2]))
    end
    return total_difference, schedule, -1
end

"""Execute the full disentangling schedule on a `CAMPS` state."""
function disentangle!(camps::CAMPS, ::DisentangleStrategy{:full}; mask=nothing, 
    index_set=Clifford2IndexSet(:entangle), criterion=DisentangleCriterion(:entangle), svd_kwargs...)
    N = length(camps)

    if mask===nothing
        mask = trues(N,N)
    end

    costs = 0

    best_clif_nums, differences, cost = disentangler(camps.mps, mask;
        index_set=index_set, criterion=criterion, svd_kwargs...)
    diff, schedule = _make_full_schedule(best_clif_nums, differences)
    costs += cost

    for (clif_num, nL, nR) in schedule
        transform!(camps, Clifford2Index(clif_num), nL, nR; svd_kwargs...)
        costs += 2*(nR-nL) - 1
    end
    return diff, schedule, costs
end

"""Disentangle only pairs within a fixed site radius."""
function disentangle!(camps::CAMPS, ::DisentangleStrategy{:radius}; radius=1,
    index_set=Clifford2IndexSet(:entangle), criterion=DisentangleCriterion(:entangle), svd_kwargs...)
    N = length(camps)
    mask = falses(N,N)
    for k in 1:radius
        for i in 1:(N-k)
            mask[i, i+k] = true
        end
    end
    return disentangle!(camps, DisentangleStrategy(:full); mask=mask,  
    index_set=index_set, criterion=criterion, svd_kwargs...)
end

"""Apply a brickwork disentangling pass over alternating nearest-neighbor bonds."""
function disentangle!(camps::CAMPS, ::DisentangleStrategy{:brickwork}; 
    index_set=Clifford2IndexSet(:entangle), criterion=DisentangleCriterion(:entangle), svd_kwargs...)
    N = length(camps)

    diff = 0.0
    schedule = []
    costs = 0
    for i in  [collect(1:2:N-1)...,collect(reverse(2:2:N-1))...]
        best_clif_num, difference, cost =  disentangler(camps.mps, i, i+1;
            index_set=index_set, criterion=criterion, svd_kwargs...)
        costs += cost
        if difference > 0.0
            diff += difference
            transform!(camps, Clifford2Index(best_clif_num), i, i+1; svd_kwargs...)
            costs += 1
            push!(schedule, (best_clif_num, i, i+1))
        end
    end
    return diff, schedule, costs
end

"""Apply a forward-and-backward nearest-neighbor disentangling pass."""
function disentangle!(camps::CAMPS, ::DisentangleStrategy{:snake}; 
    index_set=Clifford2IndexSet(:entangle), criterion=DisentangleCriterion(:entangle), svd_kwargs...)
    N = length(camps)

    diff = 0.0
    schedule = []
    costs = 0

    bond_dims = dim.(linkinds(camps.mps))
    for (bond, dim) in enumerate(bond_dims)
        if dim == 1
            continue
        end

        best_clif_num, difference, cost =  disentangler(camps.mps, bond, bond+1;
            index_set=index_set, criterion=criterion, svd_kwargs...)
        costs += cost
        if difference > 0.0
            diff += difference
            transform!(camps, Clifford2Index(best_clif_num), bond, bond+1; svd_kwargs...)
            costs += 1
            push!(schedule, (best_clif_num, bond, bond+1))
        end
    end

    bond_dims = dim.(linkinds(camps.mps))
    for (bond, dim) in reverse(collect(enumerate(bond_dims)))
        if dim == 1
            continue
        end
        costs += 1
        best_clif_num, difference, cost =  disentangler(camps.mps, bond, bond+1;
            index_set=index_set, criterion=criterion, svd_kwargs...)
        costs += cost
        if difference > 0.0
            diff += difference
            transform!(camps, Clifford2Index(best_clif_num), bond, bond+1; svd_kwargs...)
            costs += 1
            push!(schedule, (best_clif_num, bond, bond+1))
        end
    end

    return diff, schedule, costs
end

"""A no-op disentangling strategy."""
function disentangle!(camps::CAMPS, ::DisentangleStrategy{:none}; kwargs...)
    return 0.0, [], 0
end

"""Repeat a disentangling strategy until the improvement falls below `min_diff`."""
function disentangle!(camps::CAMPS, strategy::DisentangleStrategys, max_iter::Int; min_diff=0.0, kwargs...)
    diffs = Float64[]
    schedules = []
    costss = Int[]

    for _ in 1:max_iter
        diff, schedule, costs = disentangle!(camps, strategy; kwargs...)
        push!(diffs, diff)
        push!(schedules, schedule)
        push!(costss, costs)
        if diff<=min_diff
             break
        end
    end
    return diffs, schedules, costss
end

"""Run the default snake-style disentangling schedule."""
function disentangle!(camps::CAMPS)
    N = length(camps)
    return disentangle!(camps, DisentangleStrategy(:snake), N;
        min_diff=1.0e-6, cutoff=1.0e-8)
end

# -- local rotation -------------------------------------------------------------------------------


_trans = [Clifford1Index(1,0) Clifford1Index(5,0) Clifford1Index(0,0);
          Clifford1Index(1,2) Clifford1Index(5,3) Clifford1Index(0,1)]

"""
Apply local single-qubit Clifford rotations that align each site's Bloch vector.

The returned distances quantify how far the local reduced states are from the
computational basis after the update.
"""
function local_rotation!(camps::CAMPS)
    N = length(camps)
    xyz = real(hcat(ITensorMPS.expect(camps.mps, ["X", "Y", "Z"])...))
    distances = Float64[]
    for n in 1:N
        val, xyz_ind = findmax(xyz[n,:].^2)
        
        indices = [1,2,3]
        deleteat!(indices,xyz_ind)
        push!(distances, sqrt(abs(1-val + xyz[n,indices[1]]^2 + xyz[n,indices[2]]^2)))

        pm_ind = nothing
        if xyz[n,xyz_ind] >= 0
            pm_ind = 1
        else
            pm_ind = 2
        end

        transform!(camps, _trans[pm_ind,xyz_ind], n)
    end
    return distances
end
