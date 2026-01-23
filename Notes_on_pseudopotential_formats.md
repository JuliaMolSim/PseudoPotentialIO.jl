
## `upf_to_json` data structure
```
pseudo_potential
    header
        element::String
        pseudo_type::String
        core_correction::Bool
        spin_orbit::Bool
        paw_core_energy::Float64
        z_valence::Float64
        l_max::Int
        cutoff_radius_index::Int
        mesh_size::Int
        number_of_wfc::Int
        number_of_proj::Int
    radial_grid::Vector{Float64}
    core_charge_density::Vector{Float64}
    local_potential::Vector{Float64}
    beta_projectors::Vector{
        label::String
        angular_momentum::Int
        radial_function::Vector{Float64}
    }
    D_ion::Vector{Float64}
    augmentation::Vector{
        i::Int
        j::Int
        angular_momentum::Int
        radial_function::Vector{Float64}
    }
    atomic_wave_functions::Vector{
        label::String,
        angular_momentum::Int,
        occupation::Float64,
        radial_function::Vector{Float64}
    }
    paw_data
        aug_integrals::Vector{Float64}
        aug_multipoles::Vector{Float64}
        ae_wfc::Vector{
            radial_function::Vector{Float64}
            angular_momentum:Int
        }
        ps_wfc::Vector{
            radial_function::Vector{Float64}
            angular_momentum::Int
        }
        occupations::Vector{Float64}
        ae_core_charge_density::Vector{Float64}
        ae_local_potential::Vector{Float64}
    total_charge_density::Vector{Float64}
```

## `libpspio` data structure (norm-conserving only)
- [x] Norm-conserving (local+non-local)
- [ ] Norm-conserving (semilocal)
- [x] Multiple projectors per angular momentum
- [x] Atomic wavefunctions
- [x] Spin-orbit coupling
- [ ] Ultrasoft
- [ ] PAW

To support ultrasoft, add a block for augmentation.
To support paw, add agumentation, paw, and gipaw blocks.
```
struct PSPData
    # General data
    psp_info::PSPInfo,   # Metadata
        author::String,       # Author
        code_name::String,    # Generator name
        code_version::String, # Generator version
        code_input::String,   # Generator input
        description::String,  # Description
        time::Datetime        # Generation date
    format::Int,         # File format
    symbol::String,      # Atomic symbol
    z::Float64,          # Atomic number
    z_valence::Float64,  # Pseudo-ion charge
    ne_valence::Float64, # Number of electrons (usually == zvalence, except for ions)
    l_max::Int,          # Maximum angular momentum
    wave_eq::Int,        # Type of equation which was solved: Dirac | Scalar relativistic | Schroedinger
    total_energy::Float, # Total energy of pseudo-ion

    # Mesh
    mesh::PSPMesh, # Radial mesh on which _all_ functions are discretized
        type::Int,           # Mesh type
        a::Float64,          # Mesh parameter
        b::Float64,          # Mesh parameter
        np::Int,             # Number of points in mesh
        r::Vector{Float64},  # Mesh points
        rab::Vector{Float64} # Factor required for discrete integration: rab[i] = (dr(x)/dx)_{x=i}

    # States
    qn_to_istate::Matrix{Int}, # Lookup table for quantum number -> state index
    n_states::Int,             # Number of states / wavefunctions
    states::Vector{PSPState},  # States / wavefunctions
        qn::QuantumNumber,  # Quantum numbers (n, l, j) for this state / wavefunction
            n::Int,    # Principal quantum number
            l::Int,    # Angular momentum
            j::Float64 # Total angular momentum
        occ::Float64,       # Occupation of the electronic state
        eigenval::Float64,  # Eigenvalue of the electronic state
        label::String,      # Label describing the state (e.g. "2S" or "4D1.5")
        rc::Float64,        # Cutoff radius used for pseudo-state generation
        wf::PSPMeshFunction # Wavefunction
            mesh::PSPMesh      # In `libpspio`: a pointer to the mesh struct
            f::Vector{Float64} # Function values and derivatives on the mesh

    # Pseudo-potentials
    scheme::Int,                      # Pseudo-potential generation scheme
    n_potentials::Int,                # Number of pseudo-potentials
    potentials::Vector{PSPPotential}, # Pseudo-potentials
        qn::QuantumNumber, # Quantum numbers for the potential
        v::PSPMeshFunction # Pseudo-potential on a radial mesh

    # Non-local projectors
    n_projector_radials::Int,                   # Number of non-local projectors
    projector_energies::Vector{Float64}, # Dij terms for interactions between projectors
    projectors::Vector{PSPProjector},    # Non-local projectors
        qn::QuantumNumber,    # Quantum numbers for the non-local projector
        energy::Float64,      # Non-local projector energy
        proj::PSPMeshFunction # Non-local projector on a radial mesh
    projectors_l_max::Int,               # Maximum angular momentum of projectors
    l_local::Int,                        # Angular momentum channel used to generate the local potential
    v_local::Vector{PSPPotential},       # Local potential for non-local pseudo-potential form

    # Exchange-correlation data (incl. non-linear core corrections)
    xc::PSPExchangeCorrelation, # Exchange-correlation
        exchange::Int,                # Exchange functional in libxc convention
        correlation::Int,             # Correlation functional in libxc convention
        nlcc_scheme::Int,             # Scheme used for obtaining non-linear core correction
        nlcc_pf_scale::Float64,       # Prefactor for model scale (multiplies core-valence crossover radius)
        nlcc_pf_value::Float64,       # Prefactor for model amplitude (multiplies core-valence crossover value)
        nlcc_density::PSPMeshFunction # Core charge density and its derivatives on a radial mesh

    # Valence charge density
    rho_valence::PSPMeshFunction  # Valence charge density
end
```

## `quantumESPRESSO` data structure
```
struct PseudoUPF
    # General data
    generated::String  # Generator name
    author::String     # Author
    date::String       # Generation date
    comment::String    # Description
    psd::String        # Atomic symbol
    type::String       # Pseudo type (NC | SL | US / USPP | PAW | 1/r)
    rel::String        # Type of relativistic treatment (no | scalar | full)
    tvanp::Bool        # Is ultrasoft
    tcoulombp::Bool    # Is Coulomb
    nlcc::Bool         # Has non-linear core corrections
    is_gth::Bool       # Is Goedecker-Teter-Hutter
    is_multiproj::Bool # Has multiple non-local projectors per angular momentum (NC only, assumed true with US & PAW)
    dft::String        # Exchange-correlation type
    zp::Float64        # Pseudo-ion charge
    etotps::Float64    # Total energy of pseudo-ion
    ecutwfc::Float64   # Suggested wavefunction energy cutoff
    ecutrho::Float64   # Suggested charge density energy cutoff

    # Non-local projectors and atomic wavefunctions
    nv::String                  # UPF version string (e.g. 2.0.1)
    lmax::Int                   # Maximum angular momentum in non-local projectors
    lmax_rho::Int               # Maximum angular momentum in charge density (should be 2 * lmax)
    vnl::Matrix{Float64,3}      # Semi-local potential for single-projector norm-conserving pseudo-potentials
    nwfc::Int                   # Number of atomic wavefunctions
    nbeta::Int                  # Number of non-local projectors
    kbeta::Vector{Int}          # Cutoff radius index for each non-local projector
    kkbeta::Int                 # Maximum cutoff radius index for all non-local projectors
    lll::Vector{Int}            # Angular momentum of each non-local projector
    beta::Matrix{Float64}       # Non-local projectors
    els::String                 # Label of each atomic wavefunction
    els_beta::String            # Label of each non-local projector
    nchi::Vector{Int}           # Principal quantum number of each atomic wavefunciton
    lchi::Vector{Int}           # Angular momentum of each atomic wavefunction
    oc::Vector{Float64}         # Occupation of each atomic wavefunction
    epseu::Vector{Float64}      # One-particle pseudo-energy for each atomic wavefunction
    rcut_chi::Vector{Float64}   # Inner cutoff radius for each atomic wavefunction
    rcutus_chi::Vector{Float64} # Outer cutoff radius for each atomic wavefunction (ultrasoft)
    chi::Matrix{Float64}        # Atomic wavefunctions (used for initial wavefunction)
    rho_at::Vector{Float64}     # Atomic charge-density (used for initial charge density)
    mesh::Int                   # Number of points in mesh
    xmin::Float64               # Minimum `x` of the linear mesh
    rmax::Float64               # Maximum radius of the mesh
    zmesh::Float64              # Nuclear charge used for the mesh
    dx::Float64                 # Δx used for the linear mesh
    r::Vector{Float64}          # Radial grid
    rab::Vector{Float64}        # Factor required for discrete integration: rab[i] = (dr(x)/dx)_{x=i}
    rho_atc::Vector{Float64}    # Atomic core charge

    # Local potential
    lloc::Int             # Angular momentum channel used to generate the local potential
    rcloc::Float64        # Radius above which the local potential is equal to the all-electron potential
    vloc::Vector{Float64} # Local atomic potential
    dion::Matrix{Float64} # Dij terms for interactions between projectors

    # Augmentation
    q_with_l::Bool           # `qfunc` is pseudized separately for each angular momentum
    nqf::Int                 # Number of Q coefficients
    nqlc::Int                # Number of angular momenta in Q
    qqq_eps::Float64         # `qfunc` is null if its norm is < `qqq_eps`
    rinner::Vector{Float64}  # Inner cutoff radius `r_L`
    qqq::Matrix{Float64}     # Qij
    qfunc::Matrix{Float64}   # `Qij(|r|)` for `|r| > r_L` without angular momentum dependence
    qfuncl::Array{Float64,3} # `Qij(|r|)` for `|r| > r_L` with angular momentum dependence (optional, except for PAW)
    qfcoef::Array{Float64,4} # Coefficients for `Q` for `|r| < r_L`

    # Atomic wavefunctions
    has_wfc::Bool          # Has all-electron and pseudo- wavefunctions (different from chi; one for each projector)
    aewfc::Matrix{Float64} # All-electron wavefunctions for each non-local projector
    pswfc::Matrix{Float64} # Pseudo-wavefunctions for each non-local projector

    # Spin-orbit coupling
    has_so::Bool            # Has spin-orbit coupling
    nn::Vector{Float64}     # Principal quantum number for each atomic wavefunction
    rcut::Vector{Float64}   # Inner? cutoff radius for each non-local projector
    rcutus::Vector{Float64} # Outer? cutoff radius for each non-local projector (ultrasoft)
    jchi::Vector{Float64}   # Total angular momentum for each atomic wavefunction (j=l+1/2 or j=l-1/2)
    jjj::Vector{Float64}    # Total angular momentum for each non-local projector (j=l+1/2 or j=l-1/2)

    # PAW
    paw_data_format::Int # Version of the PAW data format
    tpawp::Bool          # Is PAW
    paw::PAWData
        ae_rho_atc::Vector{Float64} # All-electron core charge
        pfunc::Array{Float64,3}     # Psi_i(r) * Psi_j(r)
        pfunc_rel::Array{Float64,3} # Psi_i(r) * Psi_j(r) small component
        ptfunc::Array{Float64,3}    # "as `pfunc`, but for pseudo
        aewfc_resl::Matrix{Float64} # "as `pfunc_rel`, but for pseudo"
        ae_vloc::Vector{Float64}    # All-electron local potential
        oc::Vector{Float64}         # Starting occupation used in initializing becsum for each non-local projector (US is for each wfc)
        augmom::Array{Float64,3}    # Multipole all-electron "pseudo"
        raug::Float64               # Cutoff radius for "augfunction"
        iraug::Int                  # Index on radial grid close to and greater than `raug`
        lmax_aug::Int               # Maximum angular momentum of augmentation functions
        core_energy::Float64        # Constant shift to recover all-electron energy
        augshape::String            # Shape of augmentation charge
    
    # GIPAW
    has_gipaw::Bool                       # Has GIPAW data
    paw_as_gipaw::Bool                    # "EMINE"
    gipaw_data_format::Int                # Version of the GIPAW data format
    gipaw_ncore_orbitals::Int             # ?
    gipaw_core_orbital_n::Vector{Float64} # ?
    gipaw_core_orbital_l::Vector{Float64} # ?
    gipaw_wfs_nchannels::Int              # ?
    gipaw_wfs_el::Vector{String}          # ?
    gipaw_wfs_ll::Vector{Int}             # ?
    gipaw_wfs_ae::Matrix{Float64}         # ?
    gipaw_wfs_rcut::Vector{Float64}       # ?
    gipaw_wfs_rcutus::Vector{Float64}     # ?
    gipaw_wfs_ps::Matrix{Float64}         # ?
end
```

# Unified Pseudopotential Format, v.2.0.1

The Unified Pseudopotential Format (UPF), currently at v.2.0.1,
is designed to store different kinds of pseudopotentials:

- Norm-conserving (NC) pseudopotentials (PP) in nonlocal form
- As above, in both semilocal (SL) and nonlocal (NL) form
- Ultrasoft (US) PP (aka Vanderbilt)
- Projector Augmented Waves (PAW) datasets
- Additional data for the GIPAW reconstruction of all-electron (AE) charge

into a flexible file format, having a XML-like structure.
All quantities depending upon r are stored in numerical form on a radial grid.

## General structure

The file is formatted, starts with a line
```
    <UPF version="2.0.1">
```
ends with a line
```
    </UPF>
```
The file contains "fields", introduced by starting and ending tags
as in the example below:
```
    <FOO>
      (content of field FOO)
    </FOO>
```
The opening tag FOO may have attributes, as in the following example:
```
    <FOO attr="BAR">
      (content of field FOO)
    </FOO>
```
Fields may contain numeric data, character strings, comments,
or other fields (subfields). Comments are introduced by `<!--`,
terminated by `-->`. The following restrictions apply:
- Tags are case-sensitive and should be capitalized
- Tags can contain only letters and digits
- Spaces are not allowed between `<` and the tag name
- Maximum line length is 80 characters
- Comments must be terminated on the same line it is opened
- Numeric data in fields should be readable using fortran free format

Note that
- Blank spaces preceding tags are ignored
- Trailing characters after a tag in the same line are ignored
- Blank lines inside a field are ignored.
- All attributes must be between quotes (attr="..." or attr ='...')

## Pre-defined fields

First-level fields defined in v.2.0.1 of the UPF format are
([] = optional, or present only in some cases):

- [ PP_INFO ]
- PP_HEADER
- PP_MESH]
- [ PP_NLCC ]
- PP_LOCAL
- [ PP_SEMILOCAL ]
- PP_NONLOCAL
- PP_PSWFC
- [ PP_FULL_WFC ]
- PP_RHOATOM
- [ PP_PAW ]

These fields must appear in the order specified above.
Fields that are not defined are ignored.

## Field specifications

All quantities are in atomic Rydberg units: e^2=2, m=1/2, hbar=1. Lengths
are in Bohr (0.529177 A), energies in Ry (13.6058 eV). Potentials are
multiplied by e so they have the units of energy. 

Arrays are written in the Fortran way (leftmost index runs faster) and can be
written in free format, not exceeding the limit of 80 characters per line.


### PP_INFO

The PP_INFO field should contain any piece of information that is deemed useful
to re-generate the pseudopotential. It is not meant to contain data that is
actually read and used. The recommanded structure is the following:
```
<PP_INFO>
  Generated using XXX code v.N
  Author: Jon Doe
  Generation date: 32Oct1976
  Pseudopotential type: SL|NC|1/r|US|PAW
  Element:  Tc
  Functional:  SLA  PW   PBX  PBC
  Suggested minimum cutoff for wavefunctions:  N Ry
  Suggested minimum cutoff for charge density: M Ry
  Non-/scalar-/fully-relativistic pseudopotential
  Local potential generation info (L, rcloc, pseudization)
  Pseudopotential is spin-orbit/contains GIPAW data
  Valence configuration:
  nl, pn, l, occ, Rcut, Rcut US, E pseu
  els(1),  nns(1),  lchi(1),  oc(1),  rcut(1),  rcutus(1),  epseu(1)
  ...
  els(n),  nns(n),  lchi(n),  oc(n),  rcut(n),  rcutus(n),  epseu(n)
  Generation configuration:
     as above, including all states used in generation
  Pseudization used: Martins-Troullier/RRKJ
  <PP_INPUTFILE>
    Copy of the input file used in generation
  </PP_INPUTFILE>
</PP_INFO>
```

### PP_HEADER

PP_HEADER contains important information, only as attributes. Structure:

```
<PP_HEADER attr1="value1" ... attrN="valueN">
</PP_HEADER>
```
where

| Attribute       | Meaning, allowed values |
| ---             | ---                    |
| generated       | Generation code (character) | 
| author          | Author (character)          | 
| date            | Generation date (character) | 
| comment         | Brief description (character)| 
| element         | A valid chemical symbol: H, He, Li, ... (character)| 
| pseudo_type     | Type of PP (see below): NC, SL, 1/r, US, PAW (character)|
| relativistic    | scalar, full, nonrelativistic (character)| 
| is_ultrasoft    | T if Ultrasoft (logical) |
| is_paw          | T if PAW       (logical) |
| is_coulomb      | T if the PP is just a bare Coulomb potential (logical) |
| has_so          | T if fully relativistic PP with spin-orbit terms (logical) |
| has_wfc         | T if all-electron orbitals in field PP_FULL_WFC (logical) |
| has_gipaw       | T if data for GIPAW reconstructions is present in n PP_GIPAW (logical) |
| paw_as_gipaw    | T if GIPAW reconstruction info is included (?) (logical) |
| core_correction | T if non-linear core correction is included (logical)  |
| functional      | A valid identifier for the exchange-correlation functional (max 25 characters) |
| z_valence       | Valence charge (Z_v real) |
| total_psenergy  | Total pseudo-valence energy of PP (real)|
| wfc_cutoff      | Optional: suggested PW cutoff (ecutwfc) for expansion of KS orbitals (real)|
| rho_cutoff      | Optional: suggested PW cutoff (ecutrho) for expansion of charge density (real)|
| l_max           | Max angular momentum component present in PP (integer, 0 to 3) |
| l_max_rho       | As above, in atomic charge density (PAW only) (integer) |
| l_local         | Angular momentum chosen as local part (-1 if none, 0 to l_max-1) (integer)|
| mesh_size       | Number of points (mesh) in radial grid (integer)|
| number_of_wfc   | Number of atomic pseudo-orbitals (nwfc, see PP_PSWFC, may differ from nbeta below) (integer) |
| number_of_proj  | Number of nonlocal projectors (nbeta, see PP_NONLOCAL) (integer) |


Allowed values for pseudo_type:
| Value  | Meaning  |
| ---    | ---      |
| NC     | Norm-Conserving PP, fully nonlocal form only     |
| SL     | Norm-Conserving PP, nonlocal and semilocal forms available |
| US     | Ultrasoft (Vanderbilt) PP. Implies: is_uspp="T"   |
| 1/r    | Coulomb potential. Implies: is_coulomb="T"        |
| PAW    | Projector-Augmented Wave set. Implies: is_paw="T" |

### PP_MESH

PP_MESH contains important information on the radial grid. Structure:
```
<PP_MESH dx="dx" mesh="mesh" xmin="xmin" rmax="rmax" zmesh="Zmesh">
  <PP_R>
     r(1) r(2) ...  r(mesh)
  </PP_R>
  <PP_RAB>
     rab(1) rab(2) ... rab(mesh)
  </PP_RAB>
</PP_MESH>
```
Meaning of the variables:
 - dx, mesh, xmin, rmax, Zmesh: radial grid parameters (mesh is integer,
   all others are real numbers)
 - r(1:mesh): radial grid points (a.u.). Can be one of the following:
   * $`r_i = exp(xmin) exp((i-1)*dx)/Zmesh`$ or
   * $`r_i = exp(xmin)(exp((i-1)*dx)-1)/Zmesh`$
- rmax=r(mesh)
- rab(1:mesh): factor required for discrete integration, so that
  $`\int f(r)dr = \sum_i f(i) rab(i)`$.

### PP_NLCC

Optional, needed only for PP with core corrections. Structure:
```
<PP_NLCC>
  rho_c(1) rho_c(2) ... rho_c(mesh)
</PP_NLCC>
```
Meaning of the rho_c(mesh) variable:
- core charge for nonlinear core correction. The integral is
  $`Z_c=\int \rho_c(r) r^2 dr d\Omega`$.

### PP_LOCAL

PP_LOCAL contains the local part] of the PP. Structure:
```
<PP_LOCAL>
   vloc(1) vloc(2) ... vloc(mesh)
</PP_LOCAL>
```
Meaning of the vloc(mesh) variable:
- local potential (Ry a.u.), sampled on the radial grid. Contains the
  long-range term $`-Z_v e^2/r`$.

### PP_SEMILOCAL

Optional, for NC PP with semilocal form only. Structure:
```
<PP_SEMILOCAL>
  <PP_VNL.1 L="l1" J="j1" >
       V(1,l1) V(2,l1) ... V(mesh,l1)
  </PP_VNL.1>
  <PP_VNL.2 L="l2" J="j2" >
       V(1,l2) V(2,l2) ... V(mesh,l2)
  </PP_VNL.2>
  ...
</PP_SEMILOCAL>
```
Contains nbeta sub-fields, each one containing potential V(1:mesh,l) for the
specified (integer) value of L. The attribute J (real valued) is present in
the fully relativistic case: J=0.5 if L=0, J=L+0.5 or J=L-0.5 if L>0.

### PP_NONLOCAL

PP_NONLOCAL contains the nonlocal part of the PP, namely:
- the "beta functions" (projectors) $`\beta_n(r)`$
- the "D matrix" $`D_{ij}`$ (diagonal for NC PP)
- augmentation terms (US and PAW only):
  * the "q functions" $`q_{ij}(r)`$ 
  * the multipoles, for PAW only
  * the integrals $`Q_{ij} = 4\pi\int q_{mn}(r) r^2 dr`$.

Note: for US PP, PAW, and for NC PP in the Kleinman-Bylander form,
$`V_{NL} = \sum_{ij} |\beta_i\rangle D_ {ij}\langle\beta_j|`$.

Structure:
```
<PP_NONLOCAL>
  <PP_BETA.1 index="i1" angular_momentum="l" [ tot_ang_mom="j" label="Nl" ]
             [ + optional attributes ] >
       beta(1,i1) beta(2,i1) ... beta(mesh,i1)
  </PP_BETA.1>
  <PP_BETA.2 index="i2" angular_momentum="l" [ tot_ang_mom="j" label="Nl" ]
             [ + optional attributes ] >
       beta(1,i2) beta(2,i2) ... beta(mesh,i2)
  </PP_BETA.2>
  ...

  <PP_DIJ>
     d(1:nbeta,1:nbeta)
  </PP_DIJ>
  
  <PP_AUGMENTATION q_with_l="q_with_l" nqf="nqf" nqlc="nqlc"
                   [ + optional attributes ] > 
     <PP_Q>
        qq(1:nbeta,1:nbeta)
     </PP_Q>
     <PP_MULTIPOLES>
        mom(1:nbeta,1:nbeta,0:2*l_max)
     </PP_MULTIPOLES>
     [ <PP_QFCOEF>
          qfcoef(1:nqf, 1:nqlc, 1:nbeta, 1:nbeta)
       </PP_QFCOEF>
       <PP_RINNER>
          rinner(1) ... rinner(nqlc)
       </PP_RINNER>
     ]
     <PP_QIJL.1.1.0  composite_index="nmb" angular_momentum="l"> 
         qfuncl(1,nmb,l) ... qfuncl(mesh,nmb,l)
    </PP_QIJL.1.1.0>
     <PP_QIJL.1.2.0  composite_index="nmb" angular_momentum="l">
         qfuncl(1,nmb,l) ... qfuncl(mesh,nmb,l)
     </PP_QIJL.1.2.0>
        ...
     <PP_QIJ.1.1  composite_index="nmb" > 
         qfunc(1,nmb) ... qfunc(mesh,nmb)
    </PP_QIJ.1.1>
     <PP_QIJ.1.2  composite_index="nmb" >
         qfunc(1,nmb) ... qfunc(mesh,nmb)
     </PP_QIJ.1.2>
        ...
   </PP_AUGMENTATION>
</PP_NONLOCAL>
```
| Variable  | Physical meaning |
| ---       | ---              |
| beta(i,j) |  $`\beta_j(r_i)`$ (beta functions or projectors)|
| d(i,j)    |  $`D_{ij}`$ (D matrix)     |
| qq(i,j)   |  $`Q_{ij}`$ integrals: $`Q_{ij} = 4\pi\int q_{ij}(r) r^2 dr`$ |
| mom(i,j,l)|  ??? |
| qfcoef(n,l,i,j)| Small-radius expansion coefficients of q functions|
| rinner(l) | q functions are pseudized for $`r < r_{inner}`$|
| qfuncl(n,ij,l) | q functions decomposed into multipoles, if q_with_l="T" (ij is a composite index)|
| qfunc(n,ij) | q functions $`q_{ij}(r_n)`$, if q_with_l="F" (ij is a composite index)|

qfcoef and rinner are to be considered obsolescent: they follow the convention
of the US PP code of David Vanderbilt. If nqf=0 (see below) qfcoef and rinner are not read.

| Tag       | Attribute       | Meaning, allowed values |
| ---       | ---             | ---                     |
| PP_BETA.N | index           | Second beta functions index (integer, 1 to nbeta)|
|           | angular_momentum| Angular momentum L (integer, 0 to l_max)   |
|	    | tot_ang_mom     | Total angular momentum J for fully relativistic case (real, J=0.5 if L=0, J=L+0.5 or J=L-0.5 if L>0) |
|	    | label           | Label for this beta function (character, e.g. 4S) |
|           | cutoff_radius          | Beta function is zero beyond cutoff_radius (real)|
|	    | cutoff_radius_index    | Index of cutoff_radius in the grid (integer)|
|           | ultrasoft_cutoff_radius| Optional: used during PP generation |
| PP_AUGMENTATION | q_with_l | T if q functions are decomposed into angular momentum components (logical)|
|                 | nqf      | Number of expansion coefficients for q functions (integer, may be zero) |
|                 | nqlc     | Number of angular momenta terms in q functions (integer, may be zero) |
|                 | shape    | PAW only: Shape of augmentation function (PSQ or GAUSS or BESSEL, character) |
|                 | cutoff_r | PAW only: q functions zero beyond cutoff_r (real) |
|                 | cutoff_r_index | PAW only: Index of cutoff_r in the grid (integer)|
|                 | augmentation_epsilon | PAW only: obscure parameter used during generation (real) |
|		  | l_max_aug | PAW only: Max value of L appearing in q functions (integer) |

A note on units: the beta functions and D matrix are defined as in Vanderbilt’s
US PP code and should have Bohr^{-1/2} and Ry units, respectively.
Since they enter the calculation only as (beta*D*beta), one might as well
have "dion" in Ry^{-1} units and "beta" in Ry*Bohr^{-1/2} units, in line with
what suggested by the Kleinman-Bylader transformation. Some converters actually
do exacrtly that.

***TO BE UNDERSTOOD***: qfunc is q(r) or r^2 q(r)???

TILL HERE
=========================================================================
### PP_PSWFC

```
<PP_PSWFC>
  els(1) lchi(1) oc(1)  "Wavefunction"
  chi(1,1) chi(2,1) ...  chi(mesh,1)
  ..........
  els(natwfc) lchi(natwfc) oc(natwfc)  "Wavefunction"
  chi(1,natwfc) chi(2,natwfc) ... chi(mesh,natwfc)
</PP_PSWFC>
```

    chi(mesh,i): χi(r), i-th radial atomic (pseudo-)orbital (radial part of the KS equation, multiplied by r)
    els(natwf), lchi(natwf), oc(natwf): as in PP_HEADER

PP_RHOATOM

```
<PP_RHOATOM>
   rho_a(1) rho_a(2) ... rho_a(mesh)
</PP_RHOATOM>
```

    rho_a(mesh): radial atomic (pseudo-)charge .This is 4π r2 times the true charge.

Additional Fields
PAW

If a PAW dataset is contained in the UPF file then the additional structure <PAW> is present; it contains the fields listed in the following sections.

PP_PAW_FORMAT_VERSION

<PP\_PAW\_FORMAT_VERSION>
  version number
</PP\_PAW\_FORMAT_VERSION>

Contains version number, current version is 0.1.

PP_AUGMENTATION

<PP_AUGMENTATION>
  Shape of augmentation charge:
  BESSEL | GAUSS | PSQ | ...
  r\_match\_augfun, irc  "augmentation max radius"
  lmax_aug             "augmentation max angular momentum"
  "Augmentation multipoles:"
  nb  = 1,nbeta
    nb1 = 1,nbeta
      l   = 0,lmax_aug
        augmom(nb,nb1,l)
      enddo
    enddo
  enddo
  "Augmentation functions:"
   do l = 0,lmax_aug
       do nb = 1,nbeta
       do nb1 = 1,nbeta
           if (abs(augmom(nb,nb1,l)) > 0)
              "label of the augmentation function"
              augfun(k), k = 1, mesh
           endif
       enddo
       enddo
   enddo
</PP_AUGMENTATION>

Data required to build augmentation functions:

    BESSEL|GAUSS|PSQ|…: function used to pseudize the augmentation
    r_match_augfun: range beyond which all the augmentation functions are zero (a.u.)
    irc: index of the radial grid closer to, and greater than, r_match_augfun
    augmom: multipole of augmentation channel (nb,nb1); it is computed as:
    mlnb,nb1 =∫ dr r2 rl (χAEnb(r) χAEnb1(r) – χPSnb(r) χPSnb1(r))
    augfun: the augmentation function, it is stored only if the corresponding augmentation multipole is different from zero.

PP_AE_RHO_ATC

<PP\_AE\_RHO_ATC>
    aeccharg(k), k = 1,mesh
</PP\_AE\_RHO_ATC>

All-electron atomic density on the radial grid.

PP_AEWFC

<PP_AEWFC>
   do nb = 1,nbeta
      aewfc(k, nb), k = 1,mesh
   end do
</PP_AEWFC>

All-electron wavefunctions used for the generation of the dataset; there is one wavefunction for each beta projector.

PP_PSWFC_FULL

<PP\_PSWFC\_FULL>
   do nb = 1,nbeta
      pswfc(k, nb), k = 1,mesh
   end do
</PP\_PSWFC\_FULL>

Pseudo wavefunction used for the generation of the dataset; note that in the PP_PSWFC field only the occupied wavefunctions are stored while for PAW you need a wavefunction for each projector.

PP_AEVLOC

<PP_AEVLOC>
   do nb = 1,nbeta
      aewfc(k, nb), k = 1,mesh
   end do
</PP_AEVLOC>

All-electron local potential

PP_KDIFF

<PP_KDIFF>
  nb  = 1,nbeta
    nb1 = 1,nbeta
      kdiff(nb, nb1)
    enddo
  enddo
<PP_KDIFF>

Kinetic energy difference between all-electron and pseudo component of each augmentation channel.

PP_OCCUP

<PP_OCCUP>
   do nb = 1,nbeta
      occ(nb)
   end do
</PP_OCCUP>

Occupations of atomic orbitals.

PP_GRID_RECON

<PP\_GRID\_RECON>
   "Minimal info to reconstruct the radial grid:"
   grid%dx,   "  dx"
   grid%xmin, "  xmin"
   grid%rmax, "  rmax"
   grid%zmesh,"  zmesh"

   <PP\_SQRT\_R>
   grid%sqr(k), k=1,mesh
   </PP\_SQRT\_R>
</PP\_GRID\_RECON>

Addition data necessary to accurately reconstruct the radial grid used for the dataset generation.

GIPAW

GIPAW additional data is necessary to reconstruct all-electron charge density using the gipaw.x program included in QE distribution.

PP_GIPAW_FORMAT_VERSION

<PP\_GIPAW\_FORMAT_VERSION>
  version number
</PP\_GIPAW\_FORMAT_VERSION>

Contains version number, current version is 0.1.

GIPAW_CORE_ORBITALS

<GIPAW\_CORE\_ORBITALS>
n\_core\_orbitals "number of core orbitals"
    <GIPAW\_CORE\_ORBITAL>
       n, l  "orbital n and l quantum numbers"
       core_orbital(k), k = 1,mesh
    </GIPAW\_CORE\_ORBITAL>

Repeated for each core orbital

</GIPAW\_CORE\_ORBITALS>

Core orbitals.

GIPAW_LOCAL_DATA

<GIPAW\_LOCAL\_DATA>
   <GIPAW\_VLOCAL\_AE>
      vlocal_ae(k), k = 1,mesh
   </GIPAW\_VLOCAL\_AE>
   <GIPAW\_VLOCAL\_PS>
      vlocal_ps(k), k = 1,mesh
   </GIPAW\_VLOCAL\_PS>
</GIPAW\_LOCAL\_DATA>

All electron and pseudo local potentials, sampled on the radial grid.

GIPAW_ORBITALS

<GIPAW_ORBITALS>
   <GIPAW\_AE\_ORBITAL>
         el(nb), ll(nb)
         wfs_ae(k,nb), k = 1, mesh
   </GIPAW\_AE\_ORBITAL>
   <GIPAW\_PS\_ORBITAL>
         rcut(nb), rcutus(nb)
         wfs_ae(k,nb), k = 1, mesh
   </GIPAW\_PS\_ORBITAL>
  Repeated for each valence orbital.
</GIPAW_ORBITALS>

    el: principal quantum number (0,1,2..)
    ll: angular momentum quantum number
    rcut: inner cutof radius (a.u.)
    rcutus: outer cutoff radius (a.u.)
    wfc_ae: all-electron wavefunction sample on radial grid
    wfc_ps: pseudo wavefunction sample on radial grid