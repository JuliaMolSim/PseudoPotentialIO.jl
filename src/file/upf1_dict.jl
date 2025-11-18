function extract_block(tag, lines)
    start_tag = "<$tag"
    end_tag = "</$tag"

    istart = findfirst(contains.(lines, start_tag))
    isnothing(istart) && error("$start_tag> not found")

    iend = findfirst(contains.(lines, end_tag))
    isnothing(iend) && error("$end_tag> not found")

    return lines[istart+1:iend-1]
end

function sanitize_header(dct)
    lines = dct["content"]
    splitlines = split.(lines)
    return Dict(
        "element" => splitlines[2][1],
        "is_ultrasoft" => splitlines[3][1] == "US",
        "core_correction" => splitlines[4][1] == "T",
        "functional" => join(splitlines[5][1:end-3], " "),
        "z_valence" => parse(Float64, splitlines[6][1]),
        "total_psenergy" => parse(Float64, splitlines[7][1]),
        "rho_cutoff" => parse(Float64, splitlines[8][2]),
        "l_max" => parse(Int, splitlines[9][1]),
        "mesh_size" => parse(Int, splitlines[10][1]),
        "number_of_wfc" => parse(Int, splitlines[11][1]),
        "number_of_proj" => parse(Int, splitlines[11][2]),
    )
end

function sanitize_beta(dct)
    lines = dct["content"]
    index, l, beta, _ = split(lines[1])
    @assert beta == "beta"
    columns = length(split(lines[2]))
    content = [v for line in lines[begin+1:end] for v in parse.(Float64, split(line))]
    return Dict(
        "type" => "real",
        "size" => length(content),
        "columns" => columns,
        "index" => parse(Int, index),
        "angular_momentum" => parse(Int, l),
        "content" => content,
    )
end

function sanitize_dij(dct)
    lines = dct["content"]
    length = max([parse(Int, split(line)[2]) for line in lines[begin+1:end]])
    content = zeros(Float64, length, length)
    for row in lines[begin+1:end]
        i, j, val = split(row)
        content[parse(Int, i), parse(Int, j)] = parse(Float64, val)
    end
    return Dict(
        "type" => "real",
        "size" => length^2,
        "content" => content,
    )
end

function sanitize(tag, dct)
    if tag == "header"
        return sanitize_header(dct)
    elseif tag == "beta"
        return sanitize_beta(dct)
    elseif tag == "dij"
        return sanitize_dij(dct)
    elseif tag == "pswfc"
        return sanitize_pswfc(dct)
    elseif any(tag .== ["r", "rab", "local", "rhoatom", "nlcc"])
        @assert haskey(dct, "content") "Malformed UPF1 file: missing content in <$tag> block"
        return [value for line in dct["content"] for value in parse.(Float64, split(line))]
    else
        return dct
    end
end

function upf1_block_to_dict(lines::Vector{String})
    dct = Dict()
    # Drop all the lines before the first tag
    inext = findfirst(startswith.(lines, "<"))
    textlines = isnothing(inext) ? lines : lines[inext:end]
    text = filter(!isempty, strip.(textlines))
    !isempty(text) && (dct["content"] = text)
    # Loop over all (nested) children
    while true
        # Find the next opening tag
        inext = findfirst(startswith.(text, "<"))
        # If there are no more tags, we are finished
        isnothing(inext) && break
        match = text[inext]
        @assert !isempty(match) "Malformed UPF1 file"
        tag = split(strip(match)[begin+1:end-1])[1]
        # Extract and parse that block
        block = extract_block(tag, lines)
        # Remove everything up to the end of the new block
        deleteat!(lines, 1:inext+length(block)-1)
        # Parse the contents of the block with a generic parser
        subdct = upf1_block_to_dict(block[start+1:end-1])
        tag = lowercase(replace(tag, "PP_" => ""))
        # Sanitize particular blocks
        subdct = sanitize(tag, subdct)
        # Store the dict
        if haskey(dct, tag)
            isa(dct[tag], Vector) || (dct[tag] = [dct[tag]])
            push!(dct[tag], subdct)
        else
            dct[tag] = subdct
        end
    end
    return dct
end

upfv1contents_to_dict(str::String) = (split(str, "\n"))
