"""
Cellular Potts Model (CPM) — Cell Sorting & Cancer Morphology

Hamiltonian:
  H = Σ_{<i,j>} J(τ_i,τ_j)·(1 - δ_{σ_i,σ_j})   contact energy
    + Σ_σ λ_V·(V(σ) - V_target)²                  volume constraint
    + Σ_σ λ_P·(P(σ) - P_target)²                  perimeter constraint

Types: 0=medium, 1=type-A (or cancer), 2=type-B

Two scenarios:
  1. Cell sorting   : mixed A/B cells spontaneously separate (differential adhesion)
  2. Cancer surface : same cell looks smooth (high λ_P, low T)
                      or jagged (low λ_P, high T)
"""

using Random
using Printf
using Plots

# ================================================================
# Parameters
# ================================================================
Base.@kwdef struct CPMParams
    NX::Int = 150
    NY::Int = 150
    N_A::Int = 30
    N_B::Int = 30
    V_target::Int = 50
    λ_V::Float64 = 5.0
    # --- perimeter (surface-tension) constraint ---
    # λ_P > 0  : penalise long perimeters → smooth/compact cell surface
    # λ_P = 0  : no perimeter constraint   → shape driven only by J and T
    λ_P::Float64 = 0.0
    # Contact energies J[τ_i+1, τ_j+1]  (0=medium, 1=A, 2=B)
    J::Matrix{Float64} = [0.0  10.0  10.0;
                          10.0  2.0  16.0;
                          10.0 16.0   2.0]
    T::Float64 = 15.0
end

const NBRS_8 = ((1,0),(-1,0),(0,1),(0,-1),(1,1),(1,-1),(-1,1),(-1,-1))

# ================================================================
# Initialization
# ================================================================

"""Compute P(σ) = number of boundary bonds for each cell."""
function compute_perimeters(lattice, n_cells, NX, NY)
    perim = zeros(Int, n_cells)
    @inbounds for x in 1:NX, y in 1:NY
        s = lattice[x, y]
        s == 0 && continue
        for (dx, dy) in NBRS_8
            if lattice[mod1(x+dx,NX), mod1(y+dy,NY)] != s
                perim[s] += 1
            end
        end
    end
    return perim
end

"""Place n_total circular cells randomly; return (lattice, cell_type, cell_vol, cell_perim, rng)."""
function init_multi_cells(p::CPMParams, rng)
    lattice  = zeros(Int, p.NX, p.NY)
    n_total  = p.N_A + p.N_B
    types    = shuffle(rng, [fill(1, p.N_A); fill(2, p.N_B)])
    cell_vol = zeros(Int, n_total)
    r = ceil(Int, sqrt(p.V_target / π)) + 1

    for cid in 1:n_total
        placed = false
        for _ in 1:200_000
            cx = rand(rng, r+2 : p.NX-r-1)
            cy = rand(rng, r+2 : p.NY-r-1)
            ok = true
            @inbounds for dx in -r:r
                ok || break
                for dy in -r:r
                    if dx^2+dy^2 ≤ r^2 && lattice[cx+dx, cy+dy] != 0
                        ok = false; break
                    end
                end
            end
            if ok
                v = 0
                @inbounds for dx in -r:r, dy in -r:r
                    if dx^2+dy^2 ≤ r^2
                        lattice[cx+dx, cy+dy] = cid
                        v += 1
                    end
                end
                cell_vol[cid] = v
                placed = true; break
            end
        end
        placed || error("Failed to place cell $cid — try fewer cells or a larger grid.")
    end

    cell_perim = compute_perimeters(lattice, n_total, p.NX, p.NY)
    return lattice, types, cell_vol, cell_perim, rng
end

"""Place one large circular cell at the centre; used for cancer morphology demo."""
function init_single_cell(p::CPMParams, rng)
    lattice = zeros(Int, p.NX, p.NY)
    cx, cy  = p.NX ÷ 2, p.NY ÷ 2
    r       = ceil(Int, sqrt(p.V_target / π))
    v = 0
    @inbounds for dx in -r:r, dy in -r:r
        if dx^2+dy^2 ≤ r^2
            xi, yi = cx+dx, cy+dy
            if 1 ≤ xi ≤ p.NX && 1 ≤ yi ≤ p.NY
                lattice[xi, yi] = 1
                v += 1
            end
        end
    end
    cell_perim = compute_perimeters(lattice, 1, p.NX, p.NY)
    return lattice, [1], [v], cell_perim, rng
end

"""
Voronoi tessellation: assign every pixel to the nearest cell centre.
Results in a fully packed lattice with no medium (type 0) pixels.
"""
function init_voronoi(p::CPMParams, rng)
    NX, NY    = p.NX, p.NY
    n_total   = p.N_A + p.N_B
    types     = shuffle(rng, [fill(1, p.N_A); fill(2, p.N_B)])
    cx        = rand(rng, 1:NX, n_total)
    cy        = rand(rng, 1:NY, n_total)

    lattice   = zeros(Int, NX, NY)
    @inbounds for x in 1:NX, y in 1:NY
        best_cid  = 1
        best_dist = typemax(Int)
        for cid in 1:n_total
            d = (x - cx[cid])^2 + (y - cy[cid])^2
            if d < best_dist
                best_dist = d
                best_cid  = cid
            end
        end
        lattice[x, y] = best_cid
    end

    cell_vol  = zeros(Int, n_total)
    @inbounds for x in 1:NX, y in 1:NY
        cell_vol[lattice[x, y]] += 1
    end

    cell_perim = compute_perimeters(lattice, n_total, NX, NY)
    return lattice, types, cell_vol, cell_perim, rng
end

function init_lattice(p::CPMParams; seed=42, single_cell=false, full_pack=false)
    rng = MersenneTwister(seed)
    if single_cell
        init_single_cell(p, rng)
    elseif full_pack
        init_voronoi(p, rng)
    else
        init_multi_cells(p, rng)
    end
end

# ================================================================
# Energy — returns (ΔH, ΔP_old, ΔP_new)
# ================================================================
@inline get_type(cell_type, s) = s == 0 ? 0 : @inbounds cell_type[s]

"""
Compute ΔH when lattice[x,y] changes from σ_old → sigma_new.
Also returns perimeter changes ΔP_old, ΔP_new needed to update cell_perim.
"""
function delta_H(lattice, cell_type, cell_vol, cell_perim, p_target_vec,
                 p::CPMParams, x::Int, y::Int, sigma_new::Int)

    sigma_old = @inbounds lattice[x, y]
    sigma_old == sigma_new && return (Inf, 0, 0)

    type_old  = get_type(cell_type, sigma_old)
    type_new  = get_type(cell_type, sigma_new)
    J         = p.J

    dH_c  = 0.0
    dp_old = 0   # ΔP(sigma_old): how much old cell's perimeter changes
    dp_new = 0   # ΔP(sigma_new): how much new cell's perimeter changes

    @inbounds for (dx, dy) in NBRS_8
        xn = mod1(x+dx, p.NX)
        yn = mod1(y+dy, p.NY)
        sn = lattice[xn, yn]
        type_n = get_type(cell_type, sn)

        # --- contact energy ---
        # remove bond (x,y)[old] — n
        if sn != sigma_old;  dH_c -= J[type_old+1, type_n+1]; end
        # add    bond (x,y)[new] — n
        if sn != sigma_new;  dH_c += J[type_new+1, type_n+1]; end

        # --- perimeter bookkeeping ---
        # σ_old loses site (x,y):
        #   neighbors of (x,y) that were IN σ_old become boundary → +1
        #   neighbors of (x,y) that were NOT in σ_old lose a boundary bond → -1
        if sn == sigma_old; dp_old += 1; else dp_old -= 1; end
        # σ_new gains site (x,y):
        #   neighbors already in σ_new lose a boundary bond → -1
        #   others gain one → +1
        if sn == sigma_new; dp_new -= 1; else dp_new += 1; end
    end

    # --- volume constraint ---
    dH_v = 0.0
    λ_V  = p.λ_V;  V_T = p.V_target
    @inbounds if sigma_old != 0
        v = cell_vol[sigma_old]
        dH_v += λ_V * ((v-1-V_T)^2 - (v-V_T)^2)
    end
    @inbounds if sigma_new != 0
        v = cell_vol[sigma_new]
        dH_v += λ_V * ((v+1-V_T)^2 - (v-V_T)^2)
    end

    # --- perimeter (surface-tension) constraint ---
    dH_p = 0.0
    if p.λ_P != 0.0
        λ_P = p.λ_P
        @inbounds if sigma_old != 0
            P = cell_perim[sigma_old];  P_T = p_target_vec[sigma_old]
            dH_p += λ_P * ((P + dp_old - P_T)^2 - (P - P_T)^2)
        end
        @inbounds if sigma_new != 0
            P = cell_perim[sigma_new];  P_T = p_target_vec[sigma_new]
            dH_p += λ_P * ((P + dp_new - P_T)^2 - (P - P_T)^2)
        end
    end

    return (dH_c + dH_v + dH_p, dp_old, dp_new)
end

# ================================================================
# Metropolis-Hastings step
# ================================================================
function mh_step!(lattice, cell_type, cell_vol, cell_perim, p_target_vec,
                  p::CPMParams, rng::AbstractRNG)
    x = rand(rng, 1:p.NX)
    y = rand(rng, 1:p.NY)

    dx, dy    = @inbounds NBRS_8[rand(rng, 1:8)]
    xn        = mod1(x+dx, p.NX)
    yn        = mod1(y+dy, p.NY)
    sigma_old = @inbounds lattice[x, y]
    sigma_new = @inbounds lattice[xn, yn]
    sigma_old == sigma_new && return

    # Prevent cell from disappearing
    @inbounds if sigma_old != 0 && cell_vol[sigma_old] <= 1; return; end

    dH, dp_old, dp_new = delta_H(lattice, cell_type, cell_vol, cell_perim,
                                  p_target_vec, p, x, y, sigma_new)

    if dH ≤ 0.0 || rand(rng) < exp(-dH / p.T)
        @inbounds lattice[x, y] = sigma_new
        @inbounds if sigma_old != 0
            cell_vol[sigma_old]  -= 1
            cell_perim[sigma_old] += dp_old
        end
        @inbounds if sigma_new != 0
            cell_vol[sigma_new]  += 1
            cell_perim[sigma_new] += dp_new
        end
    end
end

# ================================================================
# Visualization
# ================================================================
function make_snapshot(lattice, cell_type, step, total_steps, title_extra="")
    NX, NY = size(lattice)
    m = zeros(Int, NX, NY)
    @inbounds for x in 1:NX, y in 1:NY
        s = lattice[x, y]
        m[x, y] = s == 0 ? 0 : cell_type[s]
    end
    cmap = cgrad([:white, :royalblue, :tomato], [0.0, 0.5, 1.0])
    heatmap(m';
        color        = cmap,
        clims        = (0, 2),
        aspect_ratio = :equal,
        axis         = nothing,
        border       = :none,
        title        = @sprintf("%s  step=%d/%d", title_extra, step, total_steps),
        colorbar     = false,
        size         = (520, 520))
end

# ================================================================
# Generic simulation runner
# ================================================================
function run_cpm!(p::CPMParams;
                  total_steps  = 5_000_000,
                  save_every   = 250_000,
                  seed         = 42,
                  single_cell  = false,
                  full_pack    = false,
                  title        = "",
                  output_gif   = "",
                  fps          = 5)

    lattice, cell_type, cell_vol, cell_perim, rng =
        init_lattice(p; seed, single_cell, full_pack)

    # P_target = perimeter of the initial compact circle (keep shape target circular)
    p_target_vec = copy(cell_perim)

    anim = Animation()
    function snap(step)
        pl = make_snapshot(lattice, cell_type, step, total_steps, title)
        frame(anim, pl)
        display(pl)
    end

    snap(0)
    t0 = time()
    for step in 1:total_steps
        mh_step!(lattice, cell_type, cell_vol, cell_perim, p_target_vec, p, rng)
        if step % save_every == 0
            @printf("  [%s] step %8d/%d  (%.2f M/s)\n",
                    title, step, total_steps, step/(time()-t0)/1e6)
            snap(step)
        end
    end
    @printf("Finished in %.1f s\n", time()-t0)

    if !isempty(output_gif)
        gif(anim, output_gif, fps=fps)
        println("Saved → $output_gif")
    end

    return lattice, cell_type, cell_vol, cell_perim, anim
end

# ================================================================
# ================================================================
#  SCENARIO 1 — CELL SORTING
#
#  Two cell types start interleaved.
#  Differential adhesion (J_AB >> J_AA ≈ J_BB) drives spontaneous sorting.
#  Tunable parameters:
#    J_AB  : higher → stronger sorting drive
#    J_AA  : lower → same-type cells stick more
#    T     : lower → faster / more complete sorting
# ================================================================
# ================================================================
function scenario_cell_sorting(;
        NX=150, NY=150, N_A=30, N_B=30,
        J_AA = 2.0,   # same-type adhesion (low = sticky)
        J_BB = 2.0,
        J_AB = 16.0,  # cross-type repulsion (high = sorting)
        λ_V  = 5.0,
        λ_P  = 0.5,   # light perimeter constraint for roundness
        T    = 15.0,
        total_steps = 1_000_000,
        save_every  = 10_000,
        fps         = 15,
        seed = 42)

    # V_target = average area per cell so that the packed lattice is balanced
    V_target = (NX * NY) ÷ (N_A + N_B)

    println("\n" * "="^60)
    println("SCENARIO 1: CELL SORTING (fully packed — no medium)")
    @printf("  J(A,A)=%.1f  J(B,B)=%.1f  J(A,B)=%.1f  T=%.1f  V_target=%d\n",
            J_AA, J_BB, J_AB, T, V_target)
    println("="^60)

    # No medium row/column — medium pixels never appear in packed mode
    J = [0.0   0.0   0.0;
         0.0   J_AA  J_AB;
         0.0   J_AB  J_BB]

    p = CPMParams(NX=NX, NY=NY, N_A=N_A, N_B=N_B,
                  V_target=V_target, λ_V=λ_V, λ_P=λ_P,
                  J=J, T=T)

    run_cpm!(p;
        total_steps, save_every, seed,
        single_cell = false,
        full_pack   = true,
        fps,
        title       = "Cell Sorting (packed)",
        output_gif  = "gif/scenario1_sorting.gif")
end

# ================================================================
# ================================================================
#  SCENARIO 2 — CANCER CELL SURFACE MORPHOLOGY
#
#  A single large cancer cell is placed in medium.
#  Surface smoothness is controlled by:
#
#    λ_P (perimeter constraint):
#        high (e.g. 3–5) → smooth, compact membrane
#        low  (e.g. 0)   → no surface-tension → jagged protrusions
#
#    T (temperature):
#        low  (e.g. 5–10) → small thermal fluctuations → smoother
#        high (e.g. 25+)  → large fluctuations → jagged / dendritic
#
#  The function runs BOTH conditions side by side and saves a comparison.
# ================================================================
# ================================================================
function scenario_cancer_morphology(;
        NX=160, NY=160,
        V_target=600,     # large cancer cell (~radius 14 px)
        λ_V=3.0,
        J_cm=12.0,        # cell-medium adhesion (surface tension)
        total_steps=3_000_000,
        save_every=500_000,
        seed=42)

    println("\n" * "="^60)
    println("SCENARIO 2: CANCER CELL SURFACE MORPHOLOGY")
    println("="^60)

    J_single = [0.0   J_cm;
                J_cm   0.0]   # only 1 cell type, no cell-cell contact

    # --- Condition A: SMOOTH (high λ_P, low T) ---
    println("\n[3A] Smooth surface  (λ_P=4.0, T=8)")
    p_smooth = CPMParams(NX=NX, NY=NY, N_A=1, N_B=0,
                         V_target=V_target, λ_V=λ_V, λ_P=4.0,
                         J=J_single, T=8.0)
    lattice_s, ct_s, _, _, _ = run_cpm!(p_smooth;
        total_steps, save_every, seed,
        single_cell = true,
        title       = "Smooth (λ_P=4, T=8)",
        output_gif  = "gif/scenario3a_smooth.gif")

    # --- Condition B: ROUGH / JAGGED (low λ_P, high T) ---
    println("\n[3B] Jagged surface  (λ_P=0, T=28)")
    p_rough = CPMParams(NX=NX, NY=NY, N_A=1, N_B=0,
                        V_target=V_target, λ_V=λ_V, λ_P=0.0,
                        J=J_single, T=28.0)
    lattice_r, ct_r, _, _, _ = run_cpm!(p_rough;
        total_steps, save_every, seed,
        single_cell = true,
        title       = "Jagged (λ_P=0, T=28)",
        output_gif  = "gif/scenario3b_jagged.gif")

    # --- Side-by-side comparison plot ---
    function cell_image(lattice, cell_type)
        NX2, NY2 = size(lattice)
        m = zeros(Int, NX2, NY2)
        for x in 1:NX2, y in 1:NY2
            s = lattice[x, y]
            m[x, y] = s == 0 ? 0 : cell_type[s]
        end
        m
    end

    cmap = cgrad([:white, :tomato], [0.0, 1.0])
    ms = cell_image(lattice_s, ct_s)
    mr = cell_image(lattice_r, ct_r)

    pl_s = heatmap(ms'; color=cmap, clims=(0,1), aspect_ratio=:equal,
                   axis=nothing, border=:none, colorbar=false,
                   title="Smooth\n(λ_P=4.0, T=8)", titlefontsize=11)
    pl_r = heatmap(mr'; color=cmap, clims=(0,1), aspect_ratio=:equal,
                   axis=nothing, border=:none, colorbar=false,
                   title="Jagged\n(λ_P=0.0, T=28)", titlefontsize=11)

    comparison = plot(pl_s, pl_r; layout=(1,2), size=(900, 460),
                      plot_title="Cancer Cell Surface Morphology",
                      top_margin=12Plots.mm)
    display(comparison)
    savefig(comparison, "fig/scenario3_comparison.png")
    println("\nSaved comparison → fig/scenario3_comparison.png")

    return lattice_s, lattice_r
end

# ================================================================
# ================================================================
#  SCENARIO 3 — PARAMETER SWEEP
#
#  Systematically vary J(A,A), J(A,B), J(B,B) to survey sorting regimes:
#
#    J_AB >> J_AA ≈ J_BB  → strong differential adhesion → clean sorting
#    J_AB ≈ J_AA ≈ J_BB   → no driving force → random mixing
#    J_AB < J_AA ≈ J_BB   → cross-type bonds preferred → permanent mixing
#    J_AA << J_BB ≈ J_AB  → A cells strongly self-adhesive → A engulfs B
#
#  Each condition is run to `total_steps` and its final frame is saved as a GIF.
#  A side-by-side comparison PNG is also written to fig/.
# ================================================================
# ================================================================

"""
Each entry in `conditions` is a NamedTuple:
  (label, J_AA, J_BB, J_AB)
"""
function scenario_energy_sweep(;
        NX=150, NY=150, N_A=30, N_B=30,
        λ_V  = 5.0,
        λ_P  = 0.5,
        T    = 15.0,
        total_steps = 1_000_000,
        save_every  = 1_000_000,   # only save final frame by default
        fps         = 5,
        seed        = 42,
        conditions)

    println("\n" * "="^60)
    println("SCENARIO 3: PARAMETER SWEEP")
    @printf("  T=%.1f  λ_V=%.1f  λ_P=%.1f  steps=%d\n", T, λ_V, λ_P, total_steps)
    println("="^60)

    V_target = (NX * NY) ÷ (N_A + N_B)
    final_plots = []

    for (i, c) in enumerate(conditions)
        @printf("\n[3-%d] %s  J_AA=%.1f  J_BB=%.1f  J_AB=%.1f\n",
                i, c.label, c.J_AA, c.J_BB, c.J_AB)

        J = [0.0     0.0      0.0;
             0.0     c.J_AA   c.J_AB;
             0.0     c.J_AB   c.J_BB]

        p = CPMParams(NX=NX, NY=NY, N_A=N_A, N_B=N_B,
                      V_target=V_target, λ_V=λ_V, λ_P=λ_P,
                      J=J, T=T)

        tag = replace(lowercase(c.label), " " => "_")
        gif_path = "gif/sweep_$(i)_$(tag).gif"

        lattice, cell_type, _, _, _ = run_cpm!(p;
            total_steps, save_every, seed,
            full_pack  = true,
            fps,
            title      = c.label,
            output_gif = gif_path)

        # Build final-frame image for comparison
        m = zeros(Int, NX, NY)
        @inbounds for x in 1:NX, y in 1:NY
            s = lattice[x, y]
            m[x, y] = s == 0 ? 0 : cell_type[s]
        end
        cmap = cgrad([:white, :royalblue, :tomato], [0.0, 0.5, 1.0])
        pl = heatmap(m'; color=cmap, clims=(0,2),
                     aspect_ratio=:equal, axis=nothing, border=:none,
                     colorbar=false,
                     title=@sprintf("%s\nJ_AA=%.0f J_BB=%.0f J_AB=%.0f",
                                    c.label, c.J_AA, c.J_BB, c.J_AB),
                     titlefontsize=9)
        push!(final_plots, pl)
    end

    # Side-by-side comparison
    n = length(final_plots)
    ncols = min(n, 4)
    nrows = ceil(Int, n / ncols)
    comparison = plot(final_plots...;
                      layout=(nrows, ncols),
                      size=(ncols*320, nrows*340),
                      plot_title="Parameter Sweep — Cell Sorting",
                      top_margin=12Plots.mm)
    display(comparison)
    savefig(comparison, "fig/sweep_comparison.png")
    println("\nSaved comparison → fig/sweep_comparison.png")

    return comparison
end

# ================================================================
# MAIN — choose which scenario(s) to run
# ================================================================
println("""
CPM Cell Sorting & Cancer Morphology Simulator
================================================
Usage:
  julia --project=. src/cell_sorting_cpm.jl 1      # cell sorting (default params)
  julia --project=. src/cell_sorting_cpm.jl 2      # cell sorting parameter sweep
  julia --project=. src/cell_sorting_cpm.jl 3      # cancer cell morphology

Output files:
  gif/scenario1_sorting.gif        — cell sorting animation
  gif/sweep_<N>_<label>.gif        — per-condition GIFs for parameter sweep
  fig/sweep_comparison.png         — side-by-side final states (parameter sweep)
  gif/scenario3a_smooth.gif        — smooth cancer cell animation
  gif/scenario3b_jagged.gif        — jagged cancer cell animation
  fig/scenario3_comparison.png     — side-by-side final states (cancer morphology)
""")

scenario = length(ARGS) >= 1 ? ARGS[1] : error("Usage: julia src/cell_sorting_cpm.jl <1|2|3>")

# ------ Scenario 1: Cell Sorting (default parameters) ------
# J_AB >> J_AA ≈ J_BB → strong differential adhesion → clean sorting
if scenario == "1"
    scenario_cell_sorting(
        N_A=30, N_B=30,
        J_AA=2.0, J_BB=2.0, J_AB=16.0,
        T=15.0,
    )
end

# ------ Scenario 2: Cell Sorting — Parameter Sweep ------
# Vary J_AA, J_BB, J_AB to survey sorting/mixing/engulfment regimes.
# Customize `conditions` to explore any J combinations you like.
if scenario == "2"
    scenario_energy_sweep(
        T           = 10.0,
        total_steps = 1_000_000,
        save_every  = 10_000,
        fps         = 15,
        conditions  = [
            (label="Strong sorting",  J_AA=1.0,  J_BB=1.0,  J_AB=20.0),
            (label="Weak sorting",    J_AA=5.0,  J_BB=5.0,  J_AB=9.0),
            (label="Mixing",          J_AA=12.0, J_BB=12.0, J_AB=2.0),
            (label="Engulfment_A",    J_AA=1.0,  J_BB=15.0, J_AB=8.0),
        ]
    )
end

# ------ Scenario 3: Cancer Cell Morphology ------
# λ_P high + T low  → smooth surface
# λ_P low  + T high → jagged/dendritic surface
if scenario == "3"
    scenario_cancer_morphology(
        V_target=600,
        J_cm=12.0,
        total_steps=3_000_000,
        save_every=500_000,
    )
end
