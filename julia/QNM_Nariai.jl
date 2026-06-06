using DelimitedFiles
using Printf
using LinearAlgebra
using GenericSchur
using Plots
using Base.Threads

# ============================================================
# Nariai spectral solver
#
# Maple exports matrices for the quadratic pencil
#
#     M0 + i*varpi*M1 + varpi^2*M2 = 0.
#
# This Julia code solves the corresponding generalized
# eigenvalue problem, identifies modes stable under
# N = 90, 95, 100, and then filters them against the exact
# analytic Nariai / Poschl-Teller spectrum.
# ============================================================

# ------------------------------------------------------------
# Field sector and angular momentum.
#
# These must match the Maple assembler:
#
# s = 0 : massless scalar,          ell >= 0
# s = 1 : electromagnetic field,    ell >= 1
# s = 2 : gravitational field,      ell >= 2
#
# Maple filetag is:
#
#     filetag := cat("s", s, "_ell", ell)
#
# so the files are, for example:
#
#     M0-90-s0_ell0.mat
#     M1-90-s0_ell0.mat
#     M2-90-s0_ell0.mat
# ------------------------------------------------------------

const SPIN = 2
const ELL  = 3

function nariai_As(s::Integer, ell::Integer)
    if s == 0
        ell >= 0 || error("For s=0 one needs ell >= 0.")
        return ell * (ell + 1)
    elseif s == 1
        ell >= 1 || error("For s=1 one needs ell >= 1.")
        return ell * (ell + 1)
    elseif s == 2
        ell >= 2 || error("For s=2 one needs ell >= 2.")
        return ell * (ell + 1) - 2
    else
        error("Unknown spin sector s=$s. Use s=0,1,2.")
    end
end

const AS = nariai_As(SPIN, ELL)
const FILETAG = "s$(SPIN)_ell$(ELL)"

# ------------------------------------------------------------
# Absolute location of the Maple output files.
# ------------------------------------------------------------

const DATA_DIR = raw"C:\Users\david\Documents\MAPLEFILESOLDPC\SPECTRAL DENYS\07-Nariai\data"

# Julia output directory.
# Each sector is written in its own subdirectory.
const OUT_DIR = joinpath(DATA_DIR, "Julia", FILETAG)
mkpath(OUT_DIR)

# ------------------------------------------------------------
# Precision
# ------------------------------------------------------------

const FILE_DECIMAL_DIGITS = 50
const WORK_DECIMAL_DIGITS = 70

setprecision(BigFloat, WORK_DECIMAL_DIGITS; base = 10)

bf(s::AbstractString) = parse(BigFloat, s)

# Tiny Re/Im roundoff noise below this threshold is printed as 0
# in readable files. Full-precision files are unchanged.
#
# For the scalar monopole, the numerical Re(varpi) values can be around
# 1e-25 or smaller, although the exact value is zero. Hence 1e-20 is
# a convenient display threshold.
const DISPLAY_ZERO_TOL = bf("1e-20")

# Tiny real parts below this value are chopped to zero in plots only.
const PLOT_ZERO_TOL = bf("1e-12")

# Stable-mode matching tolerance between N = 90,95,100 runs.
const STABLE_RELTOL = bf("1e-4")
const STABLE_ABSTOL = bf("1e-4")

# Filtering tolerance against the exact analytic Nariai spectrum.
# Modes farther than this from the nearest analytic frequency are
# treated as off-ladder numerical roots.
const ANALYTIC_FILTER_TOL = bf("1e-8")

# Number of analytic overtones generated for comparison.
const ANALYTIC_NMAX = 80

# These must match the Maple files.
const N_VALS = (90, 95, 100)

# This code solves:
#
#     (M0 + i*varpi*M1 + varpi^2*M2)b = 0
#

matfile(prefix::AbstractString, n::Integer) =
    joinpath(DATA_DIR, "$(prefix)-$(n)-$(FILETAG).mat")

outfile(name::AbstractString) =
    joinpath(OUT_DIR, name)

function read_big_matrix(path::AbstractString)
    isfile(path) || error("Missing input file: $path")

    A = readdlm(path, BigFloat)
    A = Matrix{BigFloat}(A)

    size(A, 1) == size(A, 2) ||
        error("Matrix is not square: $path has size $(size(A))")

    return A
end

function qep_eigs(M0::AbstractMatrix{T},
                  M1::AbstractMatrix{T},
                  M2::AbstractMatrix{T}) where {T<:Real}

    n = size(M0, 1)

    @assert size(M0, 2) == n "M0 must be square"
    @assert size(M1) == (n, n) "M1 must match size(M0)"
    @assert size(M2) == (n, n) "M2 must match size(M0)"

    Tc = Complex{T}

    Z = zeros(Tc, n, n)
    I_n = Matrix{Tc}(I, n, n)

    M0c = Tc.(M0)
    M1c = Tc.(M1)
    M2c = Tc.(M2)

    # gamma = i gives
    #
    #     M0 + i*varpi*M1 + varpi^2*M2 = 0.
    #
    gamma = complex(zero(T), one(T))

    # Linearization:
    #
    #     A v = varpi B v,
    #
    # with v = [b; varpi*b].
    #
    A = [Z      I_n;
        -M0c   -gamma*M1c]

    B = [I_n   Z;
         Z     M2c]

    vals = try
        eigvals(A, B)
    catch err
        @warn "eigvals(A, B) failed; falling back to eigvals(B \\ A)" exception=(err, catch_backtrace())
        eigvals(B \ A)
    end

    return Tc.(vals)
end

function finite_modes(vals)
    return [ω for ω in vals if isfinite(real(ω)) && isfinite(imag(ω))]
end

function load_and_solve(n::Integer)
    M0 = read_big_matrix(matfile("M0", n))
    M1 = read_big_matrix(matfile("M1", n))
    M2 = read_big_matrix(matfile("M2", n))

    vals = qep_eigs(M0, M1, M2)
    vals_finite = finite_modes(vals)

    if length(vals_finite) != length(vals)
        @warn "Discarded non-finite eigenvalues" n=n total=length(vals) finite=length(vals_finite)
    end

    return vals_finite
end

# ------------------------------------------------------------
# Stable-mode identification across resolutions
# ------------------------------------------------------------

# Find the nearest unused eigenvalue in vals to ω.
function nearest_mode(ω, vals, used::AbstractVector{Bool})
    jmin = 0
    dmin = BigFloat(Inf)

    for j in eachindex(vals)
        used[j] && continue

        d = abs(vals[j] - ω)

        if d < dmin
            dmin = d
            jmin = j
        end
    end

    return jmin, dmin
end

# Find the nearest eigenvalue in vals to ω, without excluding used modes.
function nearest_mode_unrestricted(ω, vals)
    jmin = 0
    dmin = BigFloat(Inf)

    for j in eachindex(vals)
        d = abs(vals[j] - ω)

        if d < dmin
            dmin = d
            jmin = j
        end
    end

    return jmin, dmin
end

# Keep only modes that match across all three runs.
function stable_modes(e1, e2, e3; reltol = STABLE_RELTOL, abstol = STABLE_ABSTOL)
    stable = Complex{BigFloat}[]

    e1m = Complex{BigFloat}[]
    e2m = Complex{BigFloat}[]
    e3m = Complex{BigFloat}[]

    err12 = BigFloat[]
    err23 = BigFloat[]
    err13 = BigFloat[]

    used2 = falses(length(e2))
    used3 = falses(length(e3))

    for ω1 in e1
        j2, d12 = nearest_mode(ω1, e2, used2)
        j2 == 0 && continue

        ω2 = e2[j2]

        j3, d23 = nearest_mode(ω2, e3, used3)
        j3 == 0 && continue

        ω3 = e3[j3]
        d13 = abs(ω1 - ω3)

        scale = max(BigFloat(1), abs(ω1), abs(ω2), abs(ω3))
        tol = max(abstol, reltol * scale)

        # Require all three pairwise distances to be acceptable.
        if max(d12, d23, d13) > tol
            continue
        end

        ωavg = (ω1 + ω2 + ω3) / 3

        push!(stable, ωavg)
        push!(e1m, ω1)
        push!(e2m, ω2)
        push!(e3m, ω3)

        push!(err12, d12)
        push!(err23, d23)
        push!(err13, d13)

        used2[j2] = true
        used3[j3] = true
    end

    # Sort by damping, so modes closest to the real axis appear first.
    p = sortperm(
        1:length(stable),
        by = k -> (abs(imag(stable[k])), real(stable[k]))
    )

    return stable[p], e1m[p], e2m[p], e3m[p], err12[p], err23[p], err13[p]
end

# ------------------------------------------------------------
# Analytic Nariai spectrum
# ------------------------------------------------------------

# Exact spectrum:
#
#     varpi_n = +/- sqrt(As - 1/4) - i*(n + 1/2).
#
# For As = 0, the two branches degenerate into
#
#     varpi = 0, -i, -2i, -3i, ...
#
function analytic_nariai_modes(nmax::Integer)
    T = BigFloat

    rootarg = T(AS) - T(1) / T(4)
    root = sqrt(complex(rootarg, zero(T)))

    modes = Complex{T}[]

    for n in 0:nmax
        damping = complex(zero(T), -(T(n) + T(1) / T(2)))

        push!(modes,  root + damping)
        push!(modes, -root + damping)
    end

    # Remove exact or near-exact duplicates.
    # This is important for the scalar monopole As = 0.
    unique_modes = Complex{T}[]

    for ω in modes
        duplicate = any(abs(ω - η) < bf("1e-40") for η in unique_modes)
        duplicate || push!(unique_modes, ω)
    end

    p = sortperm(
        1:length(unique_modes),
        by = k -> (abs(imag(unique_modes[k])), real(unique_modes[k]))
    )

    return unique_modes[p]
end

# Keep only one stable numerical mode per analytic frequency.
# For each analytic frequency, retain the closest numerical stable mode,
# provided the distance is below ANALYTIC_FILTER_TOL.
function filter_by_analytic_unique_indices(stable, analytic; tol = ANALYTIC_FILTER_TOL)
    best = Dict{Int, Tuple{Int, BigFloat}}()

    for k in eachindex(stable)
        j, d = nearest_mode_unrestricted(stable[k], analytic)
        j == 0 && continue
        d <= tol || continue

        if !haskey(best, j) || d < best[j][2]
            best[j] = (k, d)
        end
    end

    analytic_indices = sort(
        collect(keys(best)),
        by = j -> (abs(imag(analytic[j])), real(analytic[j]))
    )

    stable_indices = Int[]
    analytic_errors = BigFloat[]

    for j in analytic_indices
        push!(stable_indices, best[j][1])
        push!(analytic_errors, best[j][2])
    end

    return stable_indices, analytic_indices, analytic_errors
end

# ------------------------------------------------------------
# File output
# ------------------------------------------------------------

function save_complex_pairs(filename::AbstractString, vals)
    open(filename, "w") do io
        println(io, "# Re(varpi) Im(varpi)")
        for ω in vals
            println(io, real(ω), " ", imag(ω))
        end
    end
end

# Add a suffix before the extension of a filename.
function with_suffix(filename::AbstractString, suffix::AbstractString)
    root, ext = splitext(filename)
    isempty(ext) && (ext = ".txt")
    return root * suffix * ext
end

# Readable formatting for Re/Im values.
# Values that are numerical zero up to DISPLAY_ZERO_TOL are printed as exactly 0.
function readable_value(x::Real; chop::Bool = false)
    y = chop && abs(x) < DISPLAY_ZERO_TOL ? zero(BigFloat) : x
    return @sprintf("% .16f", y)
end

# Readable formatting for convergence errors.
function readable_error(x::Real)
    return @sprintf("%.6e", x)
end

const VALUE_WIDTH = 22
const ERROR_WIDTH = 16

function write_stable_summary_header(io)
    print(io, lpad("k", 5))
    print(io, "  ", lpad("Re(avg)", VALUE_WIDTH))
    print(io, "  ", lpad("Im(avg)", VALUE_WIDTH))
    print(io, "  ", lpad("err12", ERROR_WIDTH))
    print(io, "  ", lpad("err23", ERROR_WIDTH))
    print(io, "  ", lpad("err13", ERROR_WIDTH))
    print(io, "  ", lpad("maxerr", ERROR_WIDTH))
    println(io)
    println(io, repeat("-", 5 + 2 * (VALUE_WIDTH + 2) + 4 * (ERROR_WIDTH + 2)))
end

function write_stable_summary_row(io, k::Integer, ω, d12, d23, d13)
    maxerr = max(d12, d23, d13)

    print(io, lpad(string(k), 5))
    print(io, "  ", lpad(readable_value(real(ω); chop = true), VALUE_WIDTH))
    print(io, "  ", lpad(readable_value(imag(ω); chop = true), VALUE_WIDTH))
    print(io, "  ", lpad(readable_error(d12), ERROR_WIDTH))
    print(io, "  ", lpad(readable_error(d23), ERROR_WIDTH))
    print(io, "  ", lpad(readable_error(d13), ERROR_WIDTH))
    print(io, "  ", lpad(readable_error(maxerr), ERROR_WIDTH))
    println(io)
end

# This writes three files:
#
#   1. filename                         readable summary table
#   2. *_details.txt                    readable per-mode details
#   3. *_full_precision.txt             full-precision data
#
function save_stable_modes(filename::AbstractString,
                           stable, e1m, e2m, e3m,
                           err12, err23, err13;
                           title::AbstractString = "Stable Nariai QNM modes")

    details_filename = with_suffix(filename, "_details")
    full_filename = with_suffix(filename, "_full_precision")

    # File 1: compact readable summary.
    open(filename, "w") do io
        println(io, "# ", title)
        println(io, "# Sector: s = ", SPIN, ", ell = ", ELL, ", As = ", AS)
        println(io, "# Pencil: M0 + i*varpi*M1 + varpi^2*M2 = 0")
        println(io, "# Re/Im are printed with 16 decimal places.")
        println(io, "# Errors are printed in scientific notation with 6 digits.")
        println(io, "# In Re/Im columns only, values with |x| < ", DISPLAY_ZERO_TOL, " are printed as 0.")
        println(io, "# For all digits, see: ", basename(full_filename))
        println(io)

        write_stable_summary_header(io)

        for k in eachindex(stable)
            write_stable_summary_row(io, k, stable[k], err12[k], err23[k], err13[k])
        end
    end

    # File 2: readable detailed report, one block per mode.
    open(details_filename, "w") do io
        println(io, "# ", title, ", detailed readable report")
        println(io, "# Sector: s = ", SPIN, ", ell = ", ELL, ", As = ", AS)
        println(io, "# Pencil: M0 + i*varpi*M1 + varpi^2*M2 = 0")
        println(io, "# Re/Im are printed with 16 decimal places.")
        println(io, "# Errors are printed in scientific notation with 6 digits.")
        println(io, "# In Re/Im columns only, values with |x| < ", DISPLAY_ZERO_TOL, " are printed as 0.")
        println(io)

        for k in eachindex(stable)
            println(io, "mode ", k)
            println(io, "  avg  Re = ", readable_value(real(stable[k]); chop = true),
                        "    Im = ", readable_value(imag(stable[k]); chop = true))
            println(io, "  n1   Re = ", readable_value(real(e1m[k]); chop = true),
                        "    Im = ", readable_value(imag(e1m[k]); chop = true))
            println(io, "  n2   Re = ", readable_value(real(e2m[k]); chop = true),
                        "    Im = ", readable_value(imag(e2m[k]); chop = true))
            println(io, "  n3   Re = ", readable_value(real(e3m[k]); chop = true),
                        "    Im = ", readable_value(imag(e3m[k]); chop = true))
            println(io, "  errors: |n1-n2| = ", readable_error(err12[k]),
                        "    |n2-n3| = ", readable_error(err23[k]),
                        "    |n1-n3| = ", readable_error(err13[k]))
            println(io)
        end
    end

    # File 3: full-precision data.
    open(full_filename, "w") do io
        println(io,
            "# Re(avg) Im(avg) ",
            "Re(n1) Im(n1) ",
            "Re(n2) Im(n2) ",
            "Re(n3) Im(n3) ",
            "|n1-n2| |n2-n3| |n1-n3|"
        )

        for k in eachindex(stable)
            println(io,
                real(stable[k]), " ",
                imag(stable[k]), " ",
                real(e1m[k]), " ",
                imag(e1m[k]), " ",
                real(e2m[k]), " ",
                imag(e2m[k]), " ",
                real(e3m[k]), " ",
                imag(e3m[k]), " ",
                err12[k], " ",
                err23[k], " ",
                err13[k]
            )
        end
    end

    println(title, " written to: ", filename)
    println("Readable details written to: ", details_filename)
    println("Full-precision data written to: ", full_filename)
end

function save_analytic_modes(filename::AbstractString, modes)
    open(filename, "w") do io
        println(io, "# Analytic Nariai QNM modes")
        println(io, "# Sector: s = ", SPIN, ", ell = ", ELL, ", As = ", AS)
        println(io, "# varpi_n = +/- sqrt(As - 1/4) - i*(n + 1/2)")
        println(io, "# Re(varpi) Im(varpi)")
        for ω in modes
            println(io, real(ω), " ", imag(ω))
        end
    end

    println("Analytic Nariai benchmark written to: ", filename)
end

function save_comparison_with_analytic(filename::AbstractString,
                                       stable_physical,
                                       analytic,
                                       analytic_indices,
                                       analytic_errors)
    open(filename, "w") do io
        println(io, "# Filtered numerical modes compared with analytic Nariai modes")
        println(io, "# Sector: s = ", SPIN, ", ell = ", ELL, ", As = ", AS)
        println(io, "# k Re(num) Im(num) Re(analytic) Im(analytic) |difference|")

        for k in eachindex(stable_physical)
            j = analytic_indices[k]
            d = analytic_errors[k]
            println(io,
                k, " ",
                real(stable_physical[k]), " ",
                imag(stable_physical[k]), " ",
                real(analytic[j]), " ",
                imag(analytic[j]), " ",
                d
            )
        end
    end

    println("Numerical/analytic comparison written to: ", filename)
end

function save_off_ladder_modes(filename::AbstractString,
                               stable_raw,
                               selected_raw_indices::AbstractVector{Bool},
                               analytic)
    open(filename, "w") do io
        println(io, "# Raw stable modes not retained by the analytic-ladder filter")
        println(io, "# These are treated as numerical duplicates/artifacts.")
        println(io, "# Sector: s = ", SPIN, ", ell = ", ELL, ", As = ", AS)
        println(io, "# Filter tolerance = ", ANALYTIC_FILTER_TOL)
        println(io, "# k Re(raw) Im(raw) Re(nearest analytic) Im(nearest analytic) |difference|")

        for k in eachindex(stable_raw)
            selected_raw_indices[k] && continue

            j, d = nearest_mode_unrestricted(stable_raw[k], analytic)

            if j == 0
                println(io,
                    k, " ",
                    real(stable_raw[k]), " ",
                    imag(stable_raw[k]), " ",
                    "NaN NaN NaN"
                )
            else
                println(io,
                    k, " ",
                    real(stable_raw[k]), " ",
                    imag(stable_raw[k]), " ",
                    real(analytic[j]), " ",
                    imag(analytic[j]), " ",
                    d
                )
            end
        end
    end

    println("Off-ladder numerical roots written to: ", filename)
end

# ------------------------------------------------------------
# Plot helpers
# ------------------------------------------------------------

function plot_xy(vals; chop_real::Bool = false)
    x = Float64[]
    y = Float64[]

    tol = Float64(PLOT_ZERO_TOL)

    for ω in vals
        xr = Float64(real(ω))
        yi = Float64(imag(ω))

        isfinite(xr) && isfinite(yi) || continue

        if chop_real && abs(xr) < tol
            xr = 0.0
        end

        push!(x, xr)
        push!(y, yi)
    end

    return x, y
end

function analytic_for_plot(analytic, stable_physical)
    isempty(stable_physical) && return analytic

    yphys = Float64[Float64(imag(ω)) for ω in stable_physical]
    ymin = minimum(yphys) - 2.0
    ymax = maximum(yphys) + 1.0

    selected = [ω for ω in analytic if ymin <= Float64(imag(ω)) <= ymax]

    return isempty(selected) ? analytic : selected
end

function recommended_xlims(vals...)
    xs = Float64[]

    for arr in vals
        for ω in arr
            xr = Float64(real(ω))
            isfinite(xr) || continue
            push!(xs, xr)
        end
    end

    isempty(xs) && return (-1.0, 1.0)

    maxabsx = maximum(abs.(xs))

    # For the scalar monopole and other purely imaginary cases,
    # force a visible window around Re(varpi)=0.
    if maxabsx < Float64(PLOT_ZERO_TOL)
        return (-0.1, 0.1)
    end

    pad = max(0.1, 0.15 * maxabsx)
    return (-maxabsx - pad, maxabsx + pad)
end

function recommended_ylims(vals...)
    ys = Float64[]

    for arr in vals
        for ω in arr
            yi = Float64(imag(ω))
            isfinite(yi) || continue
            push!(ys, yi)
        end
    end

    isempty(ys) && return (-1.0, 1.0)

    ymin = minimum(ys)
    ymax = maximum(ys)
    pad = max(0.5, 0.05 * (ymax - ymin))

    return (ymin - pad, ymax + pad)
end

# ------------------------------------------------------------
# Main run
# ------------------------------------------------------------

println("Julia threads: ", Threads.nthreads())
println("Reading Maple data from: ", DATA_DIR)
println("Writing Julia output to: ", OUT_DIR)
println("Spin sector s = ", SPIN)
println("ell = ", ELL)
println("As = ", AS)
println("Using file tag: ", FILETAG)
println("Matrix pencil: M0 + i*varpi*M1 + varpi^2*M2 = 0")

# Check all required files before launching expensive solves.
for n in N_VALS
    for prefix in ("M0", "M1", "M2")
        f = matfile(prefix, n)
        isfile(f) || error("Required input file does not exist: $f")
    end
end

# Solve the three problems.
e1 = load_and_solve(N_VALS[1])
e2 = load_and_solve(N_VALS[2])
e3 = load_and_solve(N_VALS[3])

println("Number of finite modes:")
println("N = ", N_VALS[1], ": ", length(e1))
println("N = ", N_VALS[2], ": ", length(e2))
println("N = ", N_VALS[3], ": ", length(e3))

# Save the raw spectra.
save_complex_pairs(outfile("spectrum-N$(N_VALS[1])-$(FILETAG).dat"), e1)
save_complex_pairs(outfile("spectrum-N$(N_VALS[2])-$(FILETAG).dat"), e2)
save_complex_pairs(outfile("spectrum-N$(N_VALS[3])-$(FILETAG).dat"), e3)

# Identify raw stable modes across resolutions.
stable_raw, e1m_raw, e2m_raw, e3m_raw, err12_raw, err23_raw, err13_raw =
    stable_modes(e1, e2, e3;
        reltol = STABLE_RELTOL,
        abstol = STABLE_ABSTOL
    )

println("Number of raw stable modes found: ", length(stable_raw))

# Save raw stable modes for diagnostics.
save_stable_modes(
    outfile("stable_modes_raw-$(FILETAG).txt"),
    stable_raw, e1m_raw, e2m_raw, e3m_raw,
    err12_raw, err23_raw, err13_raw;
    title = "Raw stable Nariai QNM modes"
)

# Analytic benchmark.
analytic = analytic_nariai_modes(ANALYTIC_NMAX)

save_analytic_modes(
    outfile("analytic_nariai_modes-$(FILETAG).dat"),
    analytic
)

# Filter raw stable modes against the analytic ladder.
physical_indices, analytic_indices, analytic_errors =
    filter_by_analytic_unique_indices(stable_raw, analytic;
        tol = ANALYTIC_FILTER_TOL
    )

selected_raw = falses(length(stable_raw))
selected_raw[physical_indices] .= true

stable = stable_raw[physical_indices]
e1m = e1m_raw[physical_indices]
e2m = e2m_raw[physical_indices]
e3m = e3m_raw[physical_indices]
err12 = err12_raw[physical_indices]
err23 = err23_raw[physical_indices]
err13 = err13_raw[physical_indices]

println("Number of filtered physical modes found: ", length(stable))
println("Number of off-ladder raw stable modes discarded: ", length(stable_raw) - length(stable))

# Save filtered physical modes.
save_stable_modes(
    outfile("stable_modes-$(FILETAG).txt"),
    stable, e1m, e2m, e3m,
    err12, err23, err13;
    title = "Filtered physical Nariai QNM modes"
)

save_comparison_with_analytic(
    outfile("stable_vs_analytic-$(FILETAG).dat"),
    stable,
    analytic,
    analytic_indices,
    analytic_errors
)

save_off_ladder_modes(
    outfile("off_ladder_roots-$(FILETAG).dat"),
    stable_raw,
    selected_raw,
    analytic
)

# ------------------------------------------------------------
# Plots
# ------------------------------------------------------------

plotlyjs()

analytic_plot = analytic_for_plot(analytic, stable)

x1, y1 = plot_xy(e1)
x2, y2 = plot_xy(e2)
x3, y3 = plot_xy(e3)
xs, ys = plot_xy(stable; chop_real = true)
xraw, yraw = plot_xy(stable_raw; chop_real = true)
xa, ya = plot_xy(analytic_plot; chop_real = true)

# Main plot: raw spectra at the three resolutions, filtered stable modes,
# and analytic modes. The x-range is chosen from the analytic/filtered data,
# so roundoff-level Re(varpi) does not create misleading zooms.
p = scatter(x1, y1;
    ms = 2.5,
    alpha = 0.25,
    label = "N=$(N_VALS[1])"
)

scatter!(p, x2, y2;
    ms = 2.5,
    alpha = 0.25,
    label = "N=$(N_VALS[2])"
)

scatter!(p, x3, y3;
    ms = 2.5,
    alpha = 0.25,
    label = "N=$(N_VALS[3])"
)

scatter!(p, xs, ys;
    ms = 5,
    alpha = 1.0,
    label = "filtered stable"
)

scatter!(p, xa, ya;
    ms = 4,
    marker = :x,
    alpha = 1.0,
    label = "analytic"
)

xlo, xhi = recommended_xlims(analytic_plot, stable)
ylo, yhi = recommended_ylims(analytic_plot, stable)

xlims!(p, xlo, xhi)
ylims!(p, ylo, yhi)

xlabel!(p, "Re ϖ")
ylabel!(p, "Im ϖ")
title!(p, "Nariai QNMs: s=$(SPIN), ell=$(ELL)")

savefig(p, outfile("qnm_plot-$(FILETAG).png"))
savefig(p, outfile("qnm_plot-$(FILETAG).html"))

# Filtered stable plot: this is the clean plot to use in the paper.
stableplot = scatter(xs, ys;
    ms = 5,
    alpha = 1.0,
    label = "filtered stable"
)

scatter!(stableplot, xa, ya;
    ms = 4,
    marker = :x,
    alpha = 1.0,
    label = "analytic"
)

xlims!(stableplot, xlo, xhi)
ylims!(stableplot, ylo, yhi)

xlabel!(stableplot, "Re ϖ")
ylabel!(stableplot, "Im ϖ")
title!(stableplot, "Nariai QNMs: s=$(SPIN), ell=$(ELL)")

savefig(stableplot, outfile("qnm_plot_stable-$(FILETAG).png"))
savefig(stableplot, outfile("qnm_plot_stable-$(FILETAG).html"))

# Diagnostic plot: raw stable modes before analytic filtering.
# This is useful to inspect numerical duplicates and off-ladder roots.
rawplot = scatter(xraw, yraw;
    ms = 5,
    alpha = 0.7,
    label = "raw stable"
)

scatter!(rawplot, xa, ya;
    ms = 4,
    marker = :x,
    alpha = 1.0,
    label = "analytic"
)

xlims!(rawplot, xlo, xhi)
ylims!(rawplot, ylo, yhi)

xlabel!(rawplot, "Re ϖ")
ylabel!(rawplot, "Im ϖ")
title!(rawplot, "Raw stable Nariai roots: s=$(SPIN), ell=$(ELL)")

savefig(rawplot, outfile("qnm_plot_raw_stable-$(FILETAG).png"))
savefig(rawplot, outfile("qnm_plot_raw_stable-$(FILETAG).html"))

println("Done.")