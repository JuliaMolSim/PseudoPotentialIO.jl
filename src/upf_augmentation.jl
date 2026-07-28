## Typed accessors for the ultrasoft/PAW augmentation payload of a `UpfFile`.
#
# This file collects, behind a small documented API, the UPF conventions around
# augmentation charges / PAW partial waves that consumer packages otherwise each
# re-derive by hand. The conventions were cross-checked against Quantum ESPRESSO
# (upflib/uspp.f90, upflib/read_upf_new.f90 and pslibrary reference files), so any
# UpfFile consumer gets the verified reshaping without repeating that work.

"""
Augmentation-charge / PAW payload extracted from a `UpfFile` for an ultrasoft (US) or
projector-augmented-wave (PAW) pseudopotential. Returned by [`upf_augmentation`](@ref).

All quantities are stored exactly as read from the file, on the file's own radial grid
(`UpfFile.mesh.r`/`UpfFile.mesh.rab`, i.e. the FULL, untruncated file grid). In
particular the one energy-valued field, `ae_vloc`, keeps the file's Rydberg units,
consistent with the rest of the file layer (e.g. `PP_LOCAL`).
"""
struct UpfAugmentationData{T}
    """
    Radial augmentation functions `Q_ij^l(r)`, keyed by projector-index pair `(i, j)`
    with `i <= j` (mirror via `Q_ij = Q_ji` for `j < i`; only `l` channels satisfying the
    triangle/parity rule `|l_i - l_j| <= l <= l_i + l_j`, `l - l_i - l_j` even, are
    present as keys of the inner `Dict` -- this doubles as the "allowed l range" indexing
    metadata). Source: `PP_AUGMENTATION/PP_QIJL.i.j.l` (requires `q_with_l = true`;
    legacy `q_with_l = false` files -- the pre-2008 `qfcoef`/`rinner` Taylor-pseudization
    scheme -- are not decoded by this accessor, see [`upf_augmentation`](@ref)).

    CONVENTION (empirically verified against QE): despite the UPF specification text
    "4π ∫ Qij(r) r² dr", the *stored* radial array already carries whatever r²/4π/angular
    normalization is needed so that a PLAIN measure integral `Σ_r qfuncl(r)·rab(r)`
    (no extra `r²` or `4π` factor) reproduces the physical multipole moment directly --
    the same "pre-weighted" convention as e.g. `PP_RHOATOM`.
    """
    qfuncl::Dict{Tuple{Int,Int}, Dict{Int, Vector{T}}}
    """
    Multipole matrix exactly as read from the file (`q[i,j] = 4π∫Qij(r) r² dr`, the `l=0`
    moment). Source: `PP_AUGMENTATION/PP_Q`. Provided only for diagnostics/cross-
    validation against `qq_at` (below); do not use this field for the physical overlap --
    some UPF generators have historically written slightly inconsistent `PP_Q` blocks.
    """
    q::Matrix{T}
    """
    `G=0` augmentation multipole matrix, RECOMPUTED from `qfuncl` at `l=0` via plain
    `Σ_r qfuncl(r)·rab(r)` integration (`i,j` pairs with no `l=0` channel, i.e.
    `l_i != l_j`, get `qq_at = 0`, which is physically correct since the `l=0` moment
    vanishes unless `l_i = l_j`). This is what a physically-correct US/PAW overlap
    operator `S = 1 + Σ|β_i⟩ qq_at(i,j) ⟨β_j|` should use, in preference to `q`.
    """
    qq_at::Matrix{T}
    """
    Angular momentum `l_i` of each projector index `i = 1:number_of_proj`, in the file
    order of `PP_BETA` entries (same ordering as `q`/`qq_at`/`PP_DIJ`'s rows-columns).
    """
    proj_l::Vector{Int}
    """
    `true` for PAW, `false` for plain ultrasoft, from `pseudo_type == "PAW"` -- the same
    test as [`PseudoPotentialIO.is_paw`](@ref), carried along here so that consumers
    holding only this struct do not also need the `UpfFile`. Note that PAW files report
    `is_ultrasoft = true` as well (PAW data is a superset of ultrasoft data in the UPF
    schema), so that flag is not the one to branch on.
    """
    is_paw::Bool
    """
    (PAW only, `nothing` for plain ultrasoft) AE partial waves `φ_i(r)`, stored as
    `r·φ(r)` (one power of r, matching `PP_BETA`'s own `r·β(r)` misnomer convention), on
    the FULL file radial grid, ordered like `proj_l` (`i = 1:number_of_proj`). Read from
    `PP_FULL_WFC/PP_AEWFC.i` -- deliberately NOT `PP_PAW/PP_AEWFC`
    (`UpfFile.paw.aewfcs`), which some generators (e.g. pslibrary) leave empty even
    though the header advertises PAW data.

    For fully-relativistic PAW files that also ship `PP_AEWFC_REL`, this field holds the
    ordinary (large-component) `PP_AEWFC` set; the Dirac small components are in
    `ae_rel_wfcs`. Both enter the AE one-center charge density for a fully-relativistic
    treatment.
    """
    aewfcs::Union{Nothing, Vector{Vector{T}}}
    """
    Fully-relativistic PAW Dirac small components `PP_AEWFC_REL`, stored as `r·f(r)` and
    ordered exactly like `aewfcs`. `nothing` for scalar-relativistic PAW and for
    ultrasoft.
    """
    ae_rel_wfcs::Union{Nothing, Vector{Vector{T}}}
    """
    (PAW only, `nothing` for plain ultrasoft) PS partial waves `φ̃_i(r)`, same storage
    convention as `aewfcs`. Read from `PP_FULL_WFC/PP_PSWFC.i`.
    """
    pswfcs::Union{Nothing, Vector{Vector{T}}}
    """
    (PAW only, `nothing` if the file omits it) all-electron (full-Z) core density on the
    FULL file radial grid, BARE convention (not r²-weighted). Source: `PP_PAW/PP_AE_NLCC`.
    """
    ae_core_density::Union{Nothing, Vector{T}}
    """
    (PAW only, `nothing` if the file omits it) all-electron local potential on the FULL
    file radial grid, in **Rydberg** exactly as stored in the file (like `PP_LOCAL`;
    conversion to Hartree is the consumer's job). Source: `PP_PAW/PP_AE_VLOC`.
    """
    ae_vloc::Union{Nothing, Vector{T}}
    """
    (PAW only) the augmentation-sphere cutoff grid index: beyond this radius
    `aewfcs`/`pswfcs` agree to the file's write precision (bit-for-bit for e.g.
    pslibrary files; tail differences of order 1e-9 for some other generators), so
    one-center radial sums must be truncated here to avoid catastrophic cancellation
    from the (potentially large, for unbound reference channels) shared tail.
    `length(rgrid)` (no truncation) for plain ultrasoft.

    PRIMARY SOURCE: the file's own `PP_AUGMENTATION/cutoff_r_index` attribute -- QE's
    `iraug` (`upflib/read_upf_new.f90`), which QE itself uses to zero
    `pfunc`/`ptfunc`/`qfuncl` beyond `iraug`. The `iraug` attribute (same quantity,
    written by some generators instead) is used when `cutoff_r_index` is absent.

    FALLBACK: when the file provides neither attribute (or the value is out of range),
    the cutoff is instead detected empirically as the last grid index where any
    `aewfcs[beta] != pswfcs[beta]` (plus a 2-point margin). When both are available they
    are cross-checked and a disagreement beyond 5% of the grid size warns; the file
    value, when present and in range, is always authoritative.
    """
    i_cut::Int
end

"""
Build the [`UpfAugmentationData`](@ref) payload for an ultrasoft or PAW `UpfFile`
(`file.header.pseudo_type in ("US", "USPP", "PAW")`). Errors if `q_with_l` is false (the
legacy pre-2008 `qfcoef`/`rinner` augmentation Taylor-pseudization scheme is not decoded
by this accessor -- a naive read of `qijs` in that case would silently produce a wrong,
mostly-zero, augmentation).

`identifier` is used only to label warnings/errors and defaults to `file.identifier`.
"""
function upf_augmentation(file::UpfFile; identifier=file.identifier)
    pseudo_type = file.header.pseudo_type
    pseudo_type in ("US", "USPP", "PAW") ||
        error("upf_augmentation: pseudopotential '$identifier' has pseudo_type " *
              "'$pseudo_type', not one of US/USPP/PAW.")
    aug = file.nonlocal.augmentation
    isnothing(aug) && error("upf_augmentation: pseudopotential '$identifier' has no " *
                            "PP_AUGMENTATION block.")
    aug.q_with_l ||
        error("upf_augmentation: pseudopotential '$identifier' uses the legacy " *
              "`q_with_l = false` augmentation charge format (Taylor-pseudized " *
              "`qfcoef`/`rinner`), which is not supported by this accessor.")

    is_paw = pseudo_type == "PAW"
    rgrid  = file.mesh.r
    drgrid = file.mesh.rab
    nproj  = file.header.number_of_proj
    proj_l = [beta.angular_momentum for beta in file.nonlocal.betas]
    @assert length(proj_l) == nproj

    qfuncl = Dict{Tuple{Int,Int}, Dict{Int, Vector{Float64}}}()
    for qijl in aug.qijls
        if isnothing(qijl.first_index) || isnothing(qijl.second_index)
            # PP_QIJL blocks carrying only a composite index: the packing convention is
            # generator-dependent, so refuse rather than risk a silently wrong (i, j).
            error("upf_augmentation: pseudopotential '$identifier' has a PP_QIJL " *
                  "block without explicit first_index/second_index attributes, " *
                  "which is not supported by this accessor.")
        end
        i, j = qijl.first_index, qijl.second_index
        i > j && ((i, j) = (j, i))  # canonical i <= j storage; Q_ij = Q_ji
        channel = get!(qfuncl, (i, j), Dict{Int, Vector{Float64}}())
        channel[qijl.angular_momentum] = qijl.qijl
    end

    q = aug.q
    qq_at = zeros(Float64, nproj, nproj)
    maxdiff = 0.0
    for i = 1:nproj, j = i:nproj
        haskey(qfuncl, (i, j)) || continue
        channels = qfuncl[(i, j)]
        haskey(channels, 0) || continue  # l=0 moment only exists when l_i == l_j
        # Plain measure sum (not Simpson as in QE); the difference is negligible for
        # functions vanishing at both grid ends, and the 1e-5 warning threshold below
        # is calibrated to it.
        integral = sum(channels[0] .* drgrid)
        qq_at[i, j] = qq_at[j, i] = integral
        maxdiff = max(maxdiff, abs(integral - q[i, j]))
    end
    if maxdiff > 1e-5
        @warn("upf_augmentation: recomputed qq_at multipole disagrees with the UPF " *
              "file's `PP_Q` field by more than 1e-5 for pseudopotential " *
              "'$identifier' (max diff $maxdiff) -- possible data or parsing issue.")
    end

    aewfcs = ae_rel_wfcs = pswfcs = nothing
    if is_paw
        full_wfc = file.full_wfc
        n_ae = isnothing(full_wfc) ? 0 : length(full_wfc.aewfcs)
        n_ps = isnothing(full_wfc) ? 0 : length(full_wfc.pswfcs)
        # FR-PAW UPF files concatenate ordinary PP_AEWFC (1:nproj, Dirac large
        # components) with PP_AEWFC_REL (nproj+1:2nproj, Dirac small components).
        valid_ae = n_ae == nproj || (file.header.has_so && n_ae == 2nproj)
        if isnothing(full_wfc) || !valid_ae || n_ps != nproj
            error("upf_augmentation: PAW pseudopotential '$identifier' does not " *
                  "provide `PP_FULL_WFC` AE/PS partial waves in a supported layout for " *
                  "all $nproj projectors (found $n_ae AE and $n_ps PS).")
        end
        ae_records     = n_ae == 2nproj ? full_wfc.aewfcs[1:nproj] : full_wfc.aewfcs
        ae_rel_records = n_ae == 2nproj ? full_wfc.aewfcs[nproj+1:2nproj] : nothing
        by_index(wfcs) = [w.wfc for w in sort(wfcs; by=w -> w.index)]
        aewfcs      = by_index(ae_records)
        ae_rel_wfcs = isnothing(ae_rel_records) ? nothing : by_index(ae_rel_records)
        pswfcs      = by_index(full_wfc.pswfcs)
    end

    ae_core_density = nothing
    if is_paw && !isnothing(file.paw) && !isempty(file.paw.ae_nlcc)
        ae_core_density = file.paw.ae_nlcc
    end

    ae_vloc = nothing
    if is_paw && !isnothing(file.paw) && !isempty(file.paw.ae_vloc)
        length(file.paw.ae_vloc) == length(rgrid) || error(
            "upf_augmentation: PAW pseudopotential '$identifier': `PP_AE_VLOC` length " *
            "$(length(file.paw.ae_vloc)) ≠ mesh length $(length(rgrid)).")
        ae_vloc = file.paw.ae_vloc  # kept in the file's Rydberg units, like PP_LOCAL
    end

    i_cut = length(rgrid)
    if is_paw
        # By PAW construction each PS partial wave matches its AE partner at and beyond
        # the augmentation-sphere radius and differs inside it (where the nodal AE wave
        # is pseudized away):
        #
        #   grid:  |------- AE != PS -------|--------- AE == PS ---------|
        #          1                        ^ sphere edge                end
        #                                     = last differing index
        #
        # so findlast of "differs" marks the sphere edge, not the grid end. Exact `!=`
        # can overshoot slightly for generators that write the shared tail with ~1e-9
        # noise; overshooting is the safe direction (less truncation). The pre-margin
        # start value 1 is only reached in the degenerate all-channels-identical case.
        i_cut_empirical = 1
        for beta in eachindex(proj_l)
            ae, ps = aewfcs[beta], pswfcs[beta]
            idx = findlast(i -> ae[i] != ps[i], eachindex(rgrid))
            isnothing(idx) && continue  # degenerate channel: AE and PS identical everywhere
            i_cut_empirical = max(i_cut_empirical, idx)
        end
        i_cut_empirical = min(length(rgrid), i_cut_empirical + 2)  # 2-point margin

        # `cutoff_r_index` and `iraug` both carry QE's iraug; modern QE writes the
        # former, some generators only the latter. Take the first in-range value.
        i_cut_file = nothing
        for cutoff in (aug.cutoff_r_index, aug.iraug)
            isnothing(cutoff) && continue
            candidate = round(Int, cutoff)
            if 1 ≤ candidate ≤ length(rgrid)
                i_cut_file = candidate
                break
            else
                @warn("upf_augmentation: pseudopotential '$identifier': augmentation " *
                      "cutoff index $cutoff is out of range 1:$(length(rgrid)); " *
                      "falling back to the empirical augmentation-sphere cutoff " *
                      "detector.")
            end
        end

        if !isnothing(i_cut_file)
            i_cut = i_cut_file
            # Tolerance relative to the grid size: where exactly the AE/PS tails start to
            # differ numerically is generator-dependent, so a fixed margin flags healthy
            # files on fine grids. A genuine parsing/data problem is off by far more.
            tolerance = 0.05 * length(rgrid)
            if abs(i_cut_file - i_cut_empirical) > tolerance
                @warn("upf_augmentation: pseudopotential '$identifier': the file's " *
                      "augmentation cutoff index ($i_cut_file) and the empirically- " *
                      "detected augmentation-sphere cutoff ($i_cut_empirical) disagree " *
                      "by more than 5% of the radial grid -- possible data or parsing " *
                      "issue. Using the file value (authoritative).")
            end
        else
            i_cut = i_cut_empirical
        end
    end

    UpfAugmentationData(qfuncl, q, qq_at, proj_l, is_paw, aewfcs, ae_rel_wfcs,
                        pswfcs, ae_core_density, ae_vloc, i_cut)
end
