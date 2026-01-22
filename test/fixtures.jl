_FILE_EXT_LOADERS = Dict(".upf" => UpfFile,
                         ".psp8" => Psp8File,
                         ".hgh" => HghFile)

# Resolve the name of a PsP family artifact into its directory
# or return the input if it's already a directory
function _resolve_family(family_name_or_dir::AbstractString)
    isdir(family_name_or_dir) && return family_name_or_dir
    try
        return @artifact_str(family_name_or_dir)
    catch
        error("PsP family $family does not exist")
    end
end

function _find_pseudos(dir::AbstractString)
    (_, _, filenames) = first(walkdir(dir))
    psp_filenames = filter(f -> !startswith(f, "."), filenames)
    psp_filenames = filter(f -> lowercase(splitext(f)[2]) in keys(_FILE_EXT_LOADERS),
                           psp_filenames)
    return psp_filenames
end


TEST_FAMILIES = (
    "hgh_lda_hgh",
    "hgh_pbe_upf",
    "gbrv_lda_1.5_upf",
    "sg15_2022.02.06_upf",
    "sssp_pbe_precision_1.1.2_upf",
    "pd_nc_sr_pbe_standard_0.4.1_upf",
    "pd_nc_sr_pbe_standard_0.4.1_psp8",
    "pd_nc_fr_pbesol_standard_0.4_psp8",
)
TEST_FILEPATHS = []
for family in TEST_FAMILIES
    family_dir = _resolve_family(family)
    psp_filenames = _find_pseudos(family_dir)
    psp_filepaths = map(filename -> joinpath(family_dir, filename), psp_filenames)
    append!(TEST_FILEPATHS, psp_filepaths)
end

UPF_FAMILIES = (
    "gbrv_pbe_1.5_upf",
    "hgh_lda_upf",
    "sg15_2022.02.06_upf",
    "sssp_pbe_efficiency_1.1.2_upf",
)
UPF_FILEPATHS = []
for family in UPF_FAMILIES
    family_dir = _resolve_family(family)
    psp_filenames = _find_pseudos(family_dir)
    psp_filepaths = map(filename -> joinpath(family_dir, filename), psp_filenames)
    append!(UPF_FILEPATHS, psp_filepaths)
end

HGH_FAMILIES = (
    "hgh_lda_hgh",
    "hgh_pbe_hgh",
)
HGH_FILEPATHS = []
for family in HGH_FAMILIES
    family_dir = _resolve_family(family)
    psp_filenames = _find_pseudos(family_dir)
    psp_filepaths = map(filename -> joinpath(family_dir, filename), psp_filenames)
    append!(HGH_FILEPATHS, psp_filepaths)
end

PSP8_FAMILIES = (
    "pd_nc_sr_pbe_standard_0.4.1_psp8",
    "pd_nc_fr_pbesol_standard_0.4_psp8",
)
PSP8_FILEPATHS = []
for family in PSP8_FAMILIES
    family_dir = _resolve_family(family)
    psp_filenames = _find_pseudos(family_dir)
    psp_filepaths = map(filename -> joinpath(family_dir, filename), psp_filenames)
    append!(PSP8_FILEPATHS, psp_filepaths)
end

UPF1_CASES = [
    ("gbrv_lda_1.5_upf", "ag_lda_v1.4.uspp.F.UPF"),
    ("sssp_pbe_precision_1.1.2_upf", "B_pbe_v1.01.uspp.F.UPF"),
    ("gbrv_pbesol_1.5_upf", "si_pbesol_v1.uspp.F.UPF"),
    ("gbrv_pbe_1.5_upf", "mg_pbe_v1.4.uspp.F.UPF"),
]
UPF1_CASE_FILEPATHS = Dict()
for (family, filename) in UPF1_CASES
    if !haskey(UPF1_CASE_FILEPATHS, filename)
        UPF1_CASE_FILEPATHS[filename] = joinpath(_resolve_family(family), filename)
    else
        error("UPF1_CASE_FILEPATHS already contains a case for $filename")
    end
end

UPF2_CASES = [
    ("pd_nc_fr_pbesol_stringent_0.4_upf", "Mg.upf"),
    ("sssp_pbe_efficiency_1.1.2_upf", "Si.pbe-n-rrkjus_psl.1.0.0.UPF"),
    ("sssp_pbe_efficiency_1.1.2_upf", "Al.pbe-n-kjpaw_psl.1.0.0.UPF"),
    ("sssp_pbe_efficiency_1.1.2_upf", "Dy.GGA-PBE-paw-v1.0.UPF"),
    ("hgh_pbe_upf", "He.pbe-hgh.UPF"),
    ("hgh_pbe_upf", "Al.pbe-hgh.UPF"),
    ("pd_nc_sr_pbe_standard_0.4.1_upf", "H.upf"),
    ("pd_nc_sr_pbe_standard_0.4.1_upf", "Fe.upf"),
]
UPF2_CASE_FILEPATHS = Dict()
for (family, filename) in UPF2_CASES
    if !haskey(UPF2_CASE_FILEPATHS, filename)
        UPF2_CASE_FILEPATHS[filename] = joinpath(_resolve_family(family), filename)
    else
        error("UPF2_CASE_FILEPATHS already contains a case for $filename")
    end
end

PSP8_CASES = [
     ("pd_nc_sr_pbe_standard_0.4.1_psp8", "H.psp8"),
     ("pd_nc_sr_pbe_standard_0.4.1_psp8", "Fe.psp8"),
     ("pd_nc_fr_pbesol_standard_0.4_psp8", "Ti.psp8"),
     ("pd_nc_sr_pbesol_stringent_0.4.1_psp8", "Zn.psp8"),
]
PSP8_CASE_FILEPATHS = Dict()
for (family, filename) in PSP8_CASES
    if !haskey(PSP8_CASE_FILEPATHS, filename)
        PSP8_CASE_FILEPATHS[filename] = joinpath(_resolve_family(family), filename)
    else
        error("PSP8_CASE_FILEPATHS already contains a case for $filename")
    end
end

HGH_CASES = [
    ("hgh_lda_hgh", "c-q4.hgh"),
    ("hgh_lda_hgh", "ni-q18.hgh"),
]
HGH_CASE_FILEPATHS = Dict()
for (family, filename) in HGH_CASES
    if !haskey(HGH_CASE_FILEPATHS, filename)
        HGH_CASE_FILEPATHS[filename] = joinpath(_resolve_family(family), filename)
    else
        error("HGH_CASE_FILEPATHS already contains a case for $filename")
    end
end

NUMERIC_CASE_FILEPATHS = Dict(
    "B_pbe_v1.01.uspp.F.UPF" => UPF1_CASE_FILEPATHS["B_pbe_v1.01.uspp.F.UPF"],
    "si_pbesol_v1.uspp.F.UPF" => UPF1_CASE_FILEPATHS["si_pbesol_v1.uspp.F.UPF"],
    "He.pbe-hgh.UPF" => UPF2_CASE_FILEPATHS["He.pbe-hgh.UPF"],
    "Al.pbe-hgh.UPF" => UPF2_CASE_FILEPATHS["Al.pbe-hgh.UPF"],
    "Si.pbe-n-rrkjus_psl.1.0.0.UPF" => UPF2_CASE_FILEPATHS["Si.pbe-n-rrkjus_psl.1.0.0.UPF"],
    "H.upf" => UPF2_CASE_FILEPATHS["H.upf"],
    "Fe.upf" => UPF2_CASE_FILEPATHS["Fe.upf"],
    "H.psp8" => PSP8_CASE_FILEPATHS["H.psp8"],
    "Fe.psp8" => PSP8_CASE_FILEPATHS["Fe.psp8"]
)

UPF2_HGH_PAIRS = [
    ("H.pbe-hgh.UPF", "h-q1.hgh"),
    ("He.pbe-hgh.UPF", "he-q2.hgh"),
    ("B.pbe-hgh.UPF", "b-q3.hgh"),
    ("C.pbe-hgh.UPF", "c-q4.hgh"),
    ("N.pbe-hgh.UPF", "n-q5.hgh"),
    ("O.pbe-hgh.UPF", "o-q6.hgh"),
    ("F.pbe-hgh.UPF", "f-q7.hgh"),
    ("Ne.pbe-hgh.UPF", "ne-q8.hgh"),
    ("Mg.pbe-hgh.UPF", "mg-q2.hgh"),
    ("Al.pbe-hgh.UPF", "al-q3.hgh"),
    ("Si.pbe-hgh.UPF", "si-q4.hgh"),
    ("P.pbe-hgh.UPF", "p-q5.hgh"),
    ("S.pbe-hgh.UPF", "s-q6.hgh"),
    ("Cl.pbe-hgh.UPF", "cl-q7.hgh"),
    ("Ar.pbe-hgh.UPF", "ar-q8.hgh"),
    ("Ga.pbe-hgh.UPF", "ga-q3.hgh"),
    ("Ge.pbe-hgh.UPF", "ge-q4.hgh"),
    ("As.pbe-hgh.UPF", "as-q5.hgh"),
    ("Se.pbe-hgh.UPF", "se-q6.hgh"),
    ("Br.pbe-hgh.UPF", "br-q7.hgh"),
    ("Kr.pbe-hgh.UPF", "kr-q8.hgh"),
    ("Ru.pbe-hgh.UPF", "ru-q8.hgh"),
    ("Rh.pbe-hgh.UPF", "rh-q9.hgh"),
    ("Pd.pbe-hgh.UPF", "pd-q10.hgh"),
    ("In.pbe-hgh.UPF", "in-q3.hgh"),
    ("Sn.pbe-hgh.UPF", "sn-q4.hgh"),
    ("Sb.pbe-hgh.UPF", "sb-q5.hgh"),
    ("Te.pbe-hgh.UPF", "te-q6.hgh"),
    ("I.pbe-hgh.UPF", "i-q7.hgh"),
    ("Xe.pbe-hgh.UPF", "xe-q8.hgh"),
    ("La.pbe-hgh.UPF", "la-q11.hgh"),
    ("Ta.pbe-hgh.UPF", "ta-q5.hgh"),
    ("W.pbe-hgh.UPF", "w-q6.hgh"),
    ("Re.pbe-hgh.UPF", "re-q7.hgh"),
    ("Os.pbe-hgh.UPF", "os-q8.hgh"),
    ("Ir.pbe-hgh.UPF", "ir-q9.hgh"),
    ("Pt.pbe-hgh.UPF", "pt-q10.hgh"),
    ("Tl.pbe-hgh.UPF", "tl-q3.hgh"),
    ("Pb.pbe-hgh.UPF", "pb-q4.hgh"),
    ("Bi.pbe-hgh.UPF", "bi-q5.hgh"),
    ("Po.pbe-hgh.UPF", "po-q6.hgh"),
    ("At.pbe-hgh.UPF", "at-q7.hgh"),
    ("Rn.pbe-hgh.UPF", "rn-q8.hgh")
]
UPF2_HGH_FILEPATHS = map(UPF2_HGH_PAIRS) do (upf2_filename, hgh_filename)
    (
        joinpath(_resolve_family("hgh_pbe_upf"), upf2_filename),
        joinpath(_resolve_family("hgh_pbe_hgh"), hgh_filename)
    )
end

UPF2_PSP8_ELEMENTS = ("Ag", "Al", "Ar", "As", "Au", "B", "Ba", "Be", "Bi", "Br", "C", "Ca",
                      "Cr", "Cu", "Cd", "Cl", "Co", "Cs", "F", "Fe", "Ga", "Ge", "H", "He",
                      "Hf", "Hg", "I", "In", "Ir", "K", "Kr", "La", "Li", "Lu", "Mg", "Mn",
                      "Mo", "N", "Na", "Nb", "Ne", "Ni", "O", "Os", "P", "Pb", "Pd", "Po",
                      "Pt", "Rb", "Re", "Rh", "Rn", "Ru", "S", "Sb", "Sc", "Se", "Si", "Sn",
                      "Sr", "Ta", "Tc", "Te", "Ti", "Tl", "V", "W", "Xe", "Y", "Zn", "Zr")
UPF2_PSP8_FILEPATHS = map(UPF2_PSP8_ELEMENTS) do element
    (
        element,
        joinpath(_resolve_family("pd_nc_sr_pbe_standard_0.4.1_upf"), "$(element).upf"),
        joinpath(_resolve_family("pd_nc_sr_pbe_standard_0.4.1_psp8"), "$(element).psp8")
    )
end
