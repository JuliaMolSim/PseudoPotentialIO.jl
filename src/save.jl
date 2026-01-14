function save_psp_file(path::AbstractString, psp::PsPFile, file_format::String)
    if file_format == "UPF v2.0.1"
        save_psp_file(path, psp, 2)
    elseif file_format == "UPF v1.old"
        @error "Not support save to UPF v1.old"
    else
        @error "Not support save from $(format(psp)) to $file_format"
    end
end