using DelimitedFiles
using Printf
using LinearAlgebra
using GenericSchur
using Plots
using Base.Threads

# Absolute location of the Maple output files
const DATA_DIR = raw"C:\Users\david\Documents\MAPLEFILESOLDPC\SPECTRAL DENYS\06-Nonextreme SdS\data"

# Julia output directory
const OUT_DIR = joinpath(DATA_DIR, "Julia")
mkpath(OUT_DIR)

# Precision
const FILE_DECIMAL_DIGITS = 50
const WORK_DECIMAL_DIGITS = 70

setprecision(BigFloat, WORK_DECIMAL_DIGITS; base=10)
bf(s::AbstractString) = parse(BigFloat, s)

# Tiny Re/Im roundoff noise below this threshold is printed as 0
# in the readable stable-mode files only. The full-precision file is unchanged.
const DISPLAY_ZERO_TOL = bf("1e-30")

# These must match the Maple files.
const N_VALS = (90, 95, 100)

# Maple currently uses filetag := round(1000*Lambda0)
# with Lambda0 = 1/50, hence FILETAG = 20.
const FILETAG = 110

# This code currently solves:
#
#     (M0 + i*ω*M1 + ω^2*M2)b = 0
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

    # gamma = i gives M0 + i*ω*M1 + ω^2*M2.
    # Change this line if your convention is different.
    gamma = complex(zero(T), one(T))

    A = [Z      I_n;
        -M0c   -gamma*M1c]

    B = [I_n   Z;
         Z     M2c]

    # Prefer the generalized eigenvalue problem A*v = ω*B*v.
    # Fall back to B\A only if GenericSchur/LinearAlgebra fails for this type.
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

# Keep only modes that match across all three runs.
function stable_modes(e1, e2, e3; reltol = bf("1e-4"), abstol = bf("1e-4"))
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

    # Sort by imaginary part, then real part.
    p = sortperm(1:length(stable), by = k -> (imag(stable[k]), real(stable[k])))

    return stable[p], e1m[p], e2m[p], e3m[p], err12[p], err23[p], err13[p]
end

function save_complex_pairs(filename::AbstractString, vals)
    open(filename, "w") do io
        println(io, "# Re Im")
        for ω in vals
            println(io, real(ω), " ", imag(ω))
        end
    end
end

# Add a suffix before the extension of a filename.
# Example: with_suffix("stable_modes-20.txt", "_details")
# returns "stable_modes-20_details.txt".
function with_suffix(filename::AbstractString, suffix::AbstractString)
    root, ext = splitext(filename)
    isempty(ext) && (ext = ".txt")
    return root * suffix * ext
end

# Readable formatting for Re/Im values.
# Values that are numerical zero up to DISPLAY_ZERO_TOL are printed as exactly 0.
# Change %.16f to %.18f, etc., if you want more visible digits.
function readable_value(x::Real; chop::Bool = false)
    y = chop && abs(x) < DISPLAY_ZERO_TOL ? zero(BigFloat) : x
    return @sprintf("% .16f", y)
end

# Readable formatting for convergence errors.
# Change %.6e to %.10e, etc., if you want more visible digits.
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

# This replaces the old one-line, full-precision dump.
# It now writes three files:
#   1. filename                         readable summary table
#   2. *_details.txt                    readable per-mode details
#   3. *_full_precision.txt             original full-precision data
function save_stable_modes(filename::AbstractString,
                           stable, e1m, e2m, e3m,
                           err12, err23, err13)

    details_filename = with_suffix(filename, "_details")
    full_filename = with_suffix(filename, "_full_precision")

    # File 1: compact readable summary.
    open(filename, "w") do io
        println(io, "# Stable QNM modes")
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
        println(io, "# Stable QNM modes, detailed readable report")
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

    # File 3: original full-precision data, kept for reproducibility.
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

    println("Readable stable-mode summary written to: ", filename)
    println("Readable stable-mode details written to: ", details_filename)
    println("Full-precision stable-mode data written to: ", full_filename)
end

println("Julia threads: ", Threads.nthreads())
println("Reading Maple data from: ", DATA_DIR)
println("Writing Julia output to: ", OUT_DIR)
println("Using file tag: ", FILETAG)

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



# Save the raw spectra
save_complex_pairs(outfile("spectrum-N$(N_VALS[1])-$(FILETAG).dat"), e1)
save_complex_pairs(outfile("spectrum-N$(N_VALS[2])-$(FILETAG).dat"), e2)
save_complex_pairs(outfile("spectrum-N$(N_VALS[3])-$(FILETAG).dat"), e3)

# Identify stable modes.
stable, e1m, e2m, e3m, err12, err23, err13 = stable_modes(e1, e2, e3;
    reltol = BigFloat("1e-4"),
    abstol = BigFloat("1e-4")
)

println("Number of stable modes found: ", length(stable))

# Save stable modes.
save_stable_modes(
    outfile("stable_modes-$(FILETAG).txt"),
    stable, e1m, e2m, e3m,
    err12, err23, err13
)

plotlyjs()

# Helper: convert complex BigFloat eigenvalues to Float64 plotting arrays
function plot_xy(vals)
    x = Float64.(real.(vals))
    y = Float64.(imag.(vals))

    good = isfinite.(x) .& isfinite.(y)

    return x[good], y[good]
end

x1, y1 = plot_xy(e1)
x2, y2 = plot_xy(e2)
x3, y3 = plot_xy(e3)
xs, ys = plot_xy(stable)

p = scatter(x1, y1;
    ms = 2.5,
    alpha = 0.35,
    label = string(N_VALS[1])
)

scatter!(p, x2, y2;
    ms = 2.5,
    alpha = 0.35,
    label = string(N_VALS[2])
)

scatter!(p, x3, y3;
    ms = 2.5,
    alpha = 0.35,
    label = string(N_VALS[3])
)

scatter!(p, xs, ys;
    ms = 5,
    alpha = 1.0,
    label = "stable"
)

xlabel!(p, "Re Ω ")
ylabel!(p, "Im Ω ")


savefig(p, outfile("qnm_plot-$(FILETAG).png"))
savefig(p, outfile("qnm_plot-$(FILETAG).html"))

stableplot = scatter(xs, ys;
    ms = 5,
    alpha = 1.0,
    label = "stable"
)

xlabel!(stableplot, "Re Ω")
ylabel!(stableplot, "Im Ω")


savefig(stableplot, outfile("qnm_plot_stable-$(FILETAG).png"))
savefig(stableplot, outfile("qnm_plot_stable-$(FILETAG).html"))