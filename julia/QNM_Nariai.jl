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
# N = 90, 95, 100. For As > 0 it tests the physical two-horizon Jost determinant,
# independently of any analytic frequency list. Analytic frequencies are
# generated only after selection, for a separate accuracy comparison.
# For As = 0 (the exact Nariai scalar monopole), the free-wave
# Jost Wronskian is 2i*varpi: there are no nonzero scattering QNMs.
# Endpoint-regular pencil roots are retained only as diagnostics;
# a numerical candidate for the static zero threshold is reported
# separately. This statement does not apply to near-extremal SdS.
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
# Each sector is written in its own subdirectory. Corrected outputs
# are kept separate from files produced by the original script.
const OUT_DIR = joinpath(DATA_DIR, "Julia_referee2_point5", FILETAG)
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

# Physical Jost matching tolerance. Every constituent resolution root and
# its average must pass independently. Precision agreement is also required.
const JOST_RESIDUAL_TOL = bf("1e-8")
const JOST_PRECISION_TOL = bf("1e-30")
const JOST_DECIMAL_DIGITS = (80, 120)
const ZERO_THRESHOLD_TOL = bf("1e-8")

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

# Scattering QNM spectrum for the physical sectors with As > 0:
#
#     varpi_n = +/- sqrt(As - 1/4) - i*(n + 1/2).
#
# At As = 0 the formal ladder 0, -i, -2i, ... is not the scattering
# QNM spectrum. The nonzero frequencies fail the selected Jost
# boundary conditions; the constant zero threshold is handled separately.
function analytic_nariai_modes(nmax::Integer)
    T = BigFloat
    AS == 0 && return Complex{T}[]

    rootarg = T(AS) - T(1) / T(4)
    root = sqrt(complex(rootarg, zero(T)))

    modes = Complex{T}[]

    for n in 0:nmax
        damping = complex(zero(T), -(T(n) + T(1) / T(2)))

        push!(modes,  root + damping)
        push!(modes, -root + damping)
    end

    # Remove exact or near-exact duplicates.
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

# Match the analytically known constant threshold without inserting a
# synthetic zero eigenvalue or accepting cancellation in the average.
# The returned index is zero when no numerical candidate meets tolerance.
function zero_threshold_index(stable, e1m, e2m, e3m; tol = ZERO_THRESHOLD_TOL)
    best_index = 0
    best_distance = BigFloat(Inf)
    for k in eachindex(stable)
        distance = max(abs(stable[k]), abs(e1m[k]), abs(e2m[k]), abs(e3m[k]))
        if distance <= tol && distance < best_distance
            best_index = k
            best_distance = distance
        end
    end
    return best_index
end

# Physical scattering test from the independently normalized Jost solutions.
# With q=(1+tanh(chi_star))/2, a+b=1, a*b=As, c=1-i*varpi,
# J_- = exp(-i*varpi*chi_star) * 2F1(a,b;c;q), J_+(chi_star)=J_-(-chi_star).
# At chi_star=0, F=2F1(a,b;c;1/2), P=-i*varpi*F+F_q/2,
# W=-2*F*P. Its normalized form is W/(abs2(F)+abs2(P)).
# No expected frequencies, overtone counts, or nearest-ladder errors enter.
# For As=0, simplify the complete expression first: F=1 and W=2i*varpi.
# At c=0,-1,... the unit-amplitude Jost functions may have normalization
# poles. Such cases are flagged for Laurent analysis, never accepted by
# multiplying W by a factor that vanishes at the exceptional frequency.
function nariai_jost_at_precision(omega, As::Integer, digits::Integer)
    return setprecision(BigFloat, digits; base=10) do
        w = Complex{BigFloat}(omega)
        ii = complex(zero(BigFloat), one(BigFloat))
        if As == 0
            W = 2*ii*w
            return (D=W/(1+abs2(w)), status="evaluated", terms=0)
        end
        q = sqrt(complex(BigFloat(1)/4-As, zero(BigFloat)))
        a, b, c = BigFloat(1)/2+q, BigFloat(1)/2-q, 1-ii*w
        if real(c) <= 0 && abs(c-round(real(c))) <= bf("1e-25")
            return (D=complex(BigFloat(NaN),BigFloat(NaN)),
                    status="unresolved normalization pole", terms=0)
        end
        term = one(w)
        F, Fq = term, zero(w)
        small = 0
        for n in 1:10000
            term *= (a+n-1)*(b+n-1)/((c+n-1)*n*2)
            F += term
            derivative_term = 2*n*term
            Fq += derivative_term
            tol = BigFloat(10)^(-digits+15)*max(BigFloat(1),abs(F),abs(Fq))
            small = max(abs(term),abs(derivative_term)) < tol ? small+1 : 0
            if small >= 8
                P = -ii*w*F+Fq/2
                denominator = abs2(F)+abs2(P)
                if denominator == 0 || !isfinite(denominator)
                    return (D=complex(BigFloat(NaN),BigFloat(NaN)),
                            status="unresolved midpoint evaluation", terms=n)
                end
                return (D=-2*F*P/denominator, status="evaluated", terms=n)
            end
        end
        return (D=complex(BigFloat(NaN),BigFloat(NaN)),
                status="unresolved series convergence", terms=10000)
    end
end

function nariai_jost_check(omega; As=AS)
    # The zero threshold is never an ordinary nonzero scattering QNM.
    if As == 0
        v = nariai_jost_at_precision(omega,As,JOST_DECIMAL_DIGITS[2])
        status = abs(omega) <= ZERO_THRESHOLD_TOL ?
                 "static threshold candidate" : "excluded nonzero free-wave root"
        return (pass=false, residual=abs(v.D), precision_difference=BigFloat(0), status=status)
    end
    lo = nariai_jost_at_precision(omega,As,JOST_DECIMAL_DIGITS[1])
    hi = nariai_jost_at_precision(omega,As,JOST_DECIMAL_DIGITS[2])
    if lo.status != "evaluated" || hi.status != "evaluated"
        return (pass=false, residual=BigFloat(NaN), precision_difference=BigFloat(NaN),
                status=lo.status != "evaluated" ? lo.status : hi.status)
    end
    delta = abs(lo.D-hi.D)
    residual = abs(hi.D)
    if delta > JOST_PRECISION_TOL || !isfinite(residual)
        return (pass=false,residual=residual,precision_difference=delta,
                status="unresolved precision agreement")
    end
    passed = residual <= JOST_RESIDUAL_TOL
    return (pass=passed,residual=residual,precision_difference=delta,
            status=passed ? "passes physical Jost check" : "fails Jost numerical tolerance")
end

function save_jost_audit(filename::AbstractString, vals)
    open(filename,"w") do io
        println(io,"# All finite input roots, with multiplicity; no analytic-frequency selection.")
        println(io,"# D=W/(norm(Jminus,Jminus_prime)*norm(Jplus,Jplus_prime)) at chi_star=0.")
        println(io,"# Tolerance ",JOST_RESIDUAL_TOL,"; decimal precisions ",JOST_DECIMAL_DIGITS)
        println(io,"# A failed numerical tolerance does not prove the absence of a nearby resonance.")
        println(io,"# index Re(varpi) Im(varpi) abs(D) precision_difference status")
        for k in eachindex(vals)
            r=nariai_jost_check(vals[k])
            println(io,k," ",real(vals[k])," ",imag(vals[k])," ",r.residual," ",r.precision_difference," ",r.status)
        end
    end
end

function physical_jost_indices(stable,e1m,e2m,e3m)
    keep=Int[]
    for k in eachindex(stable)
        # Resolution matching alone cannot establish the radiation branches.
        checks=[nariai_jost_check(w) for w in (stable[k],e1m[k],e2m[k],e3m[k])]
        all(r.pass for r in checks) && push!(keep,k)
    end
    return keep
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
        println(io, "# Jost-selected numerical modes compared with analytic Nariai modes")
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
for (n,values) in zip(N_VALS,(e1,e2,e3))
    save_jost_audit(outfile("jost_audit_raw-N$(n)-$(FILETAG).dat"),values)
end

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
    title = "Raw stable Nariai regular-pencil roots (diagnostics)"
)

if AS == 0
    # For the exact Nariai scalar monopole, W = 2i*varpi has no
    # nonzero zeros. Regularity of the factored polynomial alone
    # therefore supplies no ordinary scattering QNMs in this sector.
    analytic = analytic_nariai_modes(ANALYTIC_NMAX)
    @assert isempty(analytic)
    threshold_index = zero_threshold_index(
        stable_raw, e1m_raw, e2m_raw, e3m_raw;
        tol = ZERO_THRESHOLD_TOL
    )
    threshold_indices = threshold_index == 0 ? Int[] : [threshold_index]

    save_stable_modes(
        outfile("zero_frequency_threshold-$(FILETAG).txt"),
        stable_raw[threshold_indices],
        e1m_raw[threshold_indices], e2m_raw[threshold_indices], e3m_raw[threshold_indices],
        err12_raw[threshold_indices], err23_raw[threshold_indices], err13_raw[threshold_indices];
        title = "Numerical candidate for the constant zero threshold (not a scattering QNM)"
    )

    status_file = outfile("scattering_qnm_status-$(FILETAG).txt")
    open(status_file, "w") do io
        println(io, "# Exact Nariai scalar monopole: As = 0")
        println(io, "# Scattering boundary conditions use the selected left/right Jost solutions.")
        println(io, "# Jost Wronskian: W = 2i*varpi; no nonzero scattering QNM frequencies.")
        println(io, "# Number of ordinary scattering QNMs: 0")
        println(io, "# The analytic constant solution at varpi = 0 is a separate static threshold.")
        println(io, "# Numerical zero-threshold matching tolerance: ", ZERO_THRESHOLD_TOL)
        if threshold_index == 0
            println(io, "# Threshold status: no stable numerical candidate meets this tolerance.")
        else
            println(io, "# Threshold status: a numerical candidate meets the tolerance in all three runs.")
            println(io, "# Raw stable index: ", threshold_index)
            println(io, "# Candidate Re(varpi) Im(varpi): ",
                real(stable_raw[threshold_index]), " ", imag(stable_raw[threshold_index]))
        end
        println(io, "# Raw spectra and stable regular-pencil roots remain available as diagnostics.")
        println(io, "# In particular, nonzero formal ladder roots are not scattering QNMs.")
        println(io, "# No analytic-ladder QNM tables or plots are generated for this sector.")
        println(io, "# This classification concerns exact Nariai, not near-extremal SdS.")
    end
    println("No nonzero scattering QNMs in the exact Nariai scalar-monopole sector.")
    println("Numerical zero-threshold candidate found: ", threshold_index != 0)
    println("Scattering classification written to: ", status_file)
else
    # Select by physical horizon Jost matching before creating a frequency list.
    physical_indices = physical_jost_indices(stable_raw,e1m_raw,e2m_raw,e3m_raw)
    save_jost_audit(outfile("jost_audit_stable-$(FILETAG).dat"),stable_raw)

    # Analytic frequencies are a post-selection accuracy comparison only.
    analytic = analytic_nariai_modes(ANALYTIC_NMAX)
    save_analytic_modes(outfile("analytic_nariai_modes-$(FILETAG).dat"),analytic)
    analytic_indices = Int[]
    analytic_errors = BigFloat[]
    for k in physical_indices
        j,d = nearest_mode_unrestricted(stable_raw[k],analytic)
        push!(analytic_indices,j)
        push!(analytic_errors,d)
    end

    selected_raw = falses(length(stable_raw))
    selected_raw[physical_indices] .= true

    stable = stable_raw[physical_indices]
    e1m = e1m_raw[physical_indices]
    e2m = e2m_raw[physical_indices]
    e3m = e3m_raw[physical_indices]
    err12 = err12_raw[physical_indices]
    err23 = err23_raw[physical_indices]
    err13 = err13_raw[physical_indices]

    println("Number passing the independent physical Jost check: ", length(stable))
    println("Number of stable roots not passing the Jost check: ", length(stable_raw) - length(stable))

    # Save filtered physical modes.
    save_stable_modes(
        outfile("stable_modes-$(FILETAG).txt"),
        stable, e1m, e2m, e3m,
        err12, err23, err13;
        title = "Stable Nariai roots passing the physical Jost check"
    )

    save_comparison_with_analytic(
        outfile("stable_vs_analytic-$(FILETAG).dat"),
        stable,
        analytic,
        analytic_indices,
        analytic_errors
    )

    # Retain the complement without relabelling numerical failures as proven
    # nonresonances. The physical diagnostic, rather than ladder distance,
    # determines membership in this file.
    save_jost_audit(
        outfile("stable_roots_not_passing_jost-$(FILETAG).dat"),
        stable_raw[.!selected_raw]
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

    # Main plot: raw spectra at the three resolutions, Jost-checked stable modes,
    # and analytic modes. Axes include every raw finite root; none is clipped
    # merely because it does not match the post-selection analytic reference.
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
        label = "Jost-checked stable"
    )

    scatter!(p, xa, ya;
        ms = 4,
        marker = :x,
        alpha = 1.0,
        label = "analytic"
    )

    xlo, xhi = recommended_xlims(e1,e2,e3,stable_raw,stable)
    ylo, yhi = recommended_ylims(e1,e2,e3,stable_raw,stable)

    xlims!(p, xlo, xhi)
    ylims!(p, ylo, yhi)

    xlabel!(p, "Re ϖ")
    ylabel!(p, "Im ϖ")
    title!(p, "Nariai QNMs: s=$(SPIN), ell=$(ELL)")

    savefig(p, outfile("qnm_plot-$(FILETAG).png"))
    savefig(p, outfile("qnm_plot-$(FILETAG).html"))

    # Jost-selected plot; all raw roots are separately saved and plotted.
    stableplot = scatter(xs, ys;
        ms = 5,
        alpha = 1.0,
        label = "Jost-checked stable"
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

    # Diagnostic plot: raw stable modes before physical Jost selection.
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

    raw_xlo,raw_xhi = recommended_xlims(e1,e2,e3,stable_raw)
    raw_ylo,raw_yhi = recommended_ylims(e1,e2,e3,stable_raw)
    xlims!(rawplot, raw_xlo, raw_xhi)
    ylims!(rawplot, raw_ylo, raw_yhi)

    xlabel!(rawplot, "Re ϖ")
    ylabel!(rawplot, "Im ϖ")
    title!(rawplot, "Raw stable Nariai roots: s=$(SPIN), ell=$(ELL)")

    savefig(rawplot, outfile("qnm_plot_raw_stable-$(FILETAG).png"))
    savefig(rawplot, outfile("qnm_plot_raw_stable-$(FILETAG).html"))

end

println("Done.")