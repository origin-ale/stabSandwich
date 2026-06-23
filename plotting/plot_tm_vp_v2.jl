#!/usr/bin/env julia
# Julia port of plot_tm_vp_v2.plt, with a sliding-window power-law fit.
#
# For every curve, a*x^b is fitted over the last 2 points, then the last 3, 4,
# ... and the (weighted) reduced chi-square is recorded. The window size with
# the smallest average reduced chi-square across all curves is then used to
# fit, plot (data with y error bars on log-log axes, plus the fits and the
# reference lines 2t and t^(2/3)) and print the final parameters.
#
# Run from the project root with:  julia --project=. plotting/plot_tm_vp_v2.jl

using Plots
using Printf
using LsqFit

# ----------------------------------------------------------------------------
# Configuration (mirrors the variables at the top of the .plt file)
# ----------------------------------------------------------------------------
const N_plot        = "12"
const Nsamples_plot = "50"

# Doping channel: :z dopes ZZ, :xy dopes XX,YY. Controls the file prefix and
# the second line of the annotation.
const doping_mode = :z
const doping_suffix, doping_label =
    doping_mode == :z  ? ("z",  "Doping with π/3 on ZZ")       :
    doping_mode == :xy ? ("xy", "Doping with 3π/16 on XX,YY")  :
    error("unknown doping_mode $doping_mode")
const prefix_plot = "TMD" * doping_suffix

# Each entry is a list of block indices to plot (0-based: index 0 is the first
# data block). One image is saved per entry.
const blocksets = [
    [0, 3, 6, 9],
    [1, 4, 7, 10],
    [2, 5, 8, 11],
    # [3, 7, 11],
]

const datafile = joinpath("output", "$(prefix)_$(N_plot)_$(Nsamples_plot).txt")

# ----------------------------------------------------------------------------
# Data parsing
# ----------------------------------------------------------------------------
# Returns (blocks, titles) where blocks[i] is a (time, mag, err) NamedTuple of
# vectors for the i-th data block (0-based index i corresponds to blocks[i+1]),
# and titles[i] is the matching "p = ..." header.
function read_blocks(path)
    blocks = Vector{NamedTuple{(:t, :y, :e), NTuple{3, Vector{Float64}}}}()
    titles = String[]
    t = Float64[]; y = Float64[]; e = Float64[]
    flush() = begin
        if !isempty(t)
            push!(blocks, (t = copy(t), y = copy(y), e = copy(e)))
            empty!(t); empty!(y); empty!(e)
        end
    end
    for raw in eachline(path)
        line = strip(raw)
        if isempty(line)
            flush()                       # blank line ends a numeric run
        elseif startswith(line, "#")
            if occursin(r"p\s*=", line)   # a "# ... p = ..." header
                push!(titles, strip(replace(line, r"^#\s*" => "")))
            end
        else
            cols = split(line)
            push!(t, parse(Float64, cols[1]))
            push!(y, parse(Float64, cols[2]))
            push!(e, length(cols) >= 3 ? parse(Float64, cols[3]) : 0.0)
        end
    end
    flush()
    return blocks, titles
end

blocktitle(titles, i) = titles[i + 1]

# ----------------------------------------------------------------------------
# Power-law fit  y = a * x^b  over the last `w` points of a block
# ----------------------------------------------------------------------------
powerlaw_model(x, p) = p[1] .* x .^ p[2]

# Fit y = a*x^b to the last `w` (positive-time) points of `blk`, weighted by
# 1/sigma^2 from the error column. Reverts to an unweighted fit when any stdev
# is <= 0. Returns a NamedTuple with the parameters, their asymptotic standard
# errors, the reduced chi-square chi2/(w-2) (Inf when there are no degrees of
# freedom), whether weights were used, and the time span (tlo, thi) fitted.
function fit_lastw(blk, w; a0 = 1.0, b0 = 0.66)
    pos = blk.t .> 0                       # x^b needs x > 0
    t = blk.t[pos]; y = blk.y[pos]; σ = blk.e[pos]
    n = length(t)
    rng = (n - w + 1):n
    t = t[rng]; y = y[rng]; σ = σ[rng]

    weighted = all(σ .> 0)
    fit = weighted ? curve_fit(powerlaw_model, t, y, 1 ./ σ .^ 2, [a0, b0]) :
                     curve_fit(powerlaw_model, t, y, [a0, b0])
    a, b = fit.param
    a_err, b_err = try
        stderror(fit)
    catch
        (NaN, NaN)
    end

    r = powerlaw_model(t, [a, b]) .- y
    chi2 = weighted ? sum(r .^ 2 ./ σ .^ 2) : sum(r .^ 2)
    dof = w - 2
    redchi2 = dof > 0 ? chi2 / dof : Inf

    return (a = a, b = b, a_err = a_err, b_err = b_err,
            redchi2 = redchi2, weighted = weighted,
            tlo = minimum(t), thi = maximum(t))
end

# ----------------------------------------------------------------------------
# Plotting
# ----------------------------------------------------------------------------
const annotation = "N = $N_plot, ϕ = θ = π/4\n$doping_label\n$Nsamples_plot samples"

# Map an axis fraction f to a coordinate on a log-scaled axis [lo, hi].
logfrac(lo, hi, f) = lo * (hi / lo)^f

# Plain-decimal label for a tick value (no scientific notation, no trailing zeros).
# Needed because the :plain formatter is ignored on log-scaled axes (GR backend).
fmtdec(v) = (s = rstrip(rstrip(@sprintf("%.10f", v), '0'), '.'); isempty(s) ? "0" : s)

# (positions, labels) at every power of ten covering [lo, hi].
function decade_ticks(lo, hi)
    vals = [10.0^p for p in floor(Int, log10(lo)):ceil(Int, log10(hi))]
    return (vals, fmtdec.(vals))
end

# (positions, labels) at every multiple of `step` covering [lo, hi].
function multiple_ticks(lo, hi, step)
    vals = collect((ceil(lo / step) * step):step:(floor(hi / step) * step))
    return (vals, fmtdec.(vals))
end

function main()
    blocks, titles = read_blocks(datafile)

    # All curves across every graph (blocksets cover each index exactly once).
    curve_idx = sort(unique(reduce(vcat, blocksets)))
    npos(i) = count(>(0), blocks[i + 1].t)
    maxw = minimum(npos(i) for i in curve_idx)

    # Window scan (no plotting): for each window size, average the reduced
    # chi-square over all curves and keep the one closest to 1 (the best
    # goodness-of-fit; minimizing it would trivially favor the smallest window).
    println("Window scan over ", length(curve_idx), " curves:")
    avg_redchi2 = fill(Inf, maxw)
    for w in 2:maxw
        avg_redchi2[w] = sum(fit_lastw(blocks[i + 1], w).redchi2 for i in curve_idx) /
                         length(curve_idx)
        @printf("  last %2d points:  avg reduced chi2 = %.5g\n", w, avg_redchi2[w])
    end
    wbest = argmin(abs.(avg_redchi2 .- 1))
    @printf("\nAverage reduced chi2 closest to 1 at last %d points (%.5g).\n\n",
            wbest, avg_redchi2[wbest])

    # Plot each graph using the fits over the last `wbest` points.
    for (k, blocks_idx) in enumerate(blocksets)
        nb = length(blocks_idx)
        A = zeros(nb); B = zeros(nb); Tlo = zeros(nb); Thi = zeros(nb)

        for (j, i) in enumerate(blocks_idx)
            f = fit_lastw(blocks[i + 1], wbest)
            A[j] = f.a; B[j] = f.b; Tlo[j] = f.tlo; Thi[j] = f.thi
            @printf("%s:  a = %.4g +/- %.2g,  b = %.4g +/- %.2g,  reduced chi2 = %.4g%s\n",
                    blocktitle(titles, i), f.a, f.a_err, f.b, f.b_err, f.redchi2,
                    f.weighted ? "" : "  (unweighted)")
        end

        # Determine x-range (positive values only, for log axis).
        xs = Float64[]
        for i in blocks_idx
            append!(xs, filter(>(0), blocks[i + 1].t))
        end
        xmin, xmax = minimum(xs), maximum(xs)
        ymin = 0.1
        ymax = maximum(reduce(vcat, blocks[i + 1].y for i in blocks_idx))

        plt = plot(; size = (850, 600),
                   title = "XXZ Floquet dynamics", titlefontsize = 24,
                   xscale = :log10, yscale = :log10,
                   xlabel = "Time", ylabel = "Transferred magnetization",
                   guidefontsize = 16, legendfontsize = 16,
                   legend = :bottomright, ylims = (ymin, Inf), grid = false,
                   xticks = multiple_ticks(xmin, xmax, 5),
                   yticks = decade_ticks(ymin, ymax))

        for (j, i) in enumerate(blocks_idx)
            blk = blocks[i + 1]
            # Drop non-positive points (Gnuplot silently skips them on log axes).
            m = (blk.t .> 0) .& (blk.y .> 0)
            plot!(plt, blk.t[m], blk.y[m]; yerror = blk.e[m],
                  color = j, marker = :circle, markersize = 3,
                  label = blocktitle(titles, i))
            # Overlay the fit over the span of points it was fitted on.
            xf = exp.(range(log(Tlo[j]), log(Thi[j]); length = 200))
            plot!(plt, xf, A[j] .* xf .^ B[j];
                  color = j, linestyle = :dash, label = "")
        end

        # Reference lines.
        xref = range(xmin, xmax; length = 200)
        plot!(plt, xref, 2 .* xref;     color = :gray,      label = "2t")
        plot!(plt, xref, xref .^ 0.66;  color = :lightgray, label = "t^(2/3)")

        # Annotation at graph fractions (0.05, 0.15).
        ax = logfrac(xmin, xmax, 0.05)
        ay = logfrac(ymin, ymax, 0.15)
        annotate!(plt, ax, ay, text(annotation, 16, :left))

        outfile = joinpath("output", "$(prefix)_$(N_plot)_$(Nsamples_plot)_graph_$(k).png")
        savefig(plt, outfile)
        println("Plotted ", outfile)
    end
end

main()
