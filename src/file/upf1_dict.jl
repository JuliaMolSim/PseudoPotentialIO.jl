function extract_block(tag, lines)
    start_tag = "<$tag"
    end_tag = "</$tag"

    istart = findfirst(contains.(lines, start_tag))
    isnothing(istart) && error("$start_tag> not found")

    iend = findfirst(contains.(lines, end_tag))
    isnothing(iend) && error("$end_tag> not found")

    return lines[istart:iend]
end

function sanitize_header(dct)
    lines = pop!(dct, "lines")
    splitlines = split.(lines)
    return Dict{String,Any}("element" => splitlines[2][1],
                "is_ultrasoft" => splitlines[3][1] == "US",
                "core_correction" => splitlines[4][1] == "T",
                "functional" => join(splitlines[5][1:(end - 3)], " "),
                "z_valence" => parse(Float64, splitlines[6][1]),
                "total_psenergy" => parse(Float64, splitlines[7][1]),
                "rho_cutoff" => parse(Float64, splitlines[8][2]),
                "l_max" => parse(Int, splitlines[9][1]),
                "mesh_size" => parse(Int, splitlines[10][1]),
                "number_of_wfc" => parse(Int, splitlines[11][1]),
                "number_of_proj" => parse(Int, splitlines[11][2]))
end

function sanitize_beta(dct)
    lines = pop!(dct, "lines")
    index, l, beta, _ = split(lines[1])
    @assert lowercase(beta) == "beta"
    columns = length(split(lines[2]))
    content = [v for line in lines[(begin + 1):end] for v in parse.(Float64, split(line))]
    return Dict{String,Any}("type" => "real",
                "size" => length(content),
                "columns" => columns,
                "index" => parse(Int, index),
                "angular_momentum" => parse(Int, l),
                "content" => content)
end

function sanitize_dij(dct)
    lines = pop!(dct, "lines")
    # Determine size of matrix as the maximum 'j' index found in the second column
    len = maximum([parse(Int, split(line)[2]) for line in lines[(begin + 1):end]])
    content = zeros(Float64, len, len)
    for row in lines[(begin + 1):end]
        i, j, val = split(row)
        content[parse(Int, i), parse(Int, j)] = parse(Float64, val)
    end
    return Dict{String,Any}("columns" => len,
                "content" => content,
                "size" => len^2,
                "type" => "real")
end

function sanitize_pswfc(dct)
    lines = pop!(dct, "lines")
    dct["chi"] = []
    array = []
    for line in lines
        if occursin("Wavefunction", line)
            !isempty(array) && (dct["chi"][end]["content"] = parse.(Float64, array))
            array = []
            label, l, occupation, _ = split(line)
            push!(dct["chi"],
                  Dict{String,Any}("label" => label,
                       "n" => parse(Int, label[begin]),
                       "l" => parse(Int, l),
                       "occupation" => parse(Float64, occupation)))
        else
            append!(array, split(line))
        end
    end
    dct["chi"][end]["content"] = parse.(Float64, array)
    return dct
end

function sanitize_qij(dct)
    lines = pop!(dct, "lines")
    dct["qij"] = []
    dct["nqf"] = parse(Int, first(split(first(lines))))
    array = []
    for line in lines
        if occursin("(l(j))", line)  # Start a new block
            i, j, l_j = parse.(Int, split(line)[begin:(begin + 2)])
            push!(dct["qij"], Dict{String,Any}("i" => i,
                                               "j" => j,
                                               "l(j)" => l_j))
            array = []
        elseif occursin("Q_int", line)
            dct["qij"][end]["q_int"] = parse(Float64, first(split(line)))
        elseif occursin("<PP_QFCOEF>", line) # End of a block
            dct["qij"][end]["content"] = parse.(Float64, array)
        else
            append!(array, split(line))
        end
    end
    return dct
end

function sanitize_rinner(dct)
    lines = pop!(dct, "lines")
    content = parse.(Float64, [split(line)[2] for line in lines])
    return Dict{String,Any}("content" => content)
end

function sanitize_vector(dct)
    lines = pop!(dct, "lines")
    columns = length(split(lines[1]))
    array = [value for line in lines for value in parse.(Float64, split(line))]
    return Dict{String,Any}("columns" => columns,
                "content" => array,
                "size" => length(array),
                "type" => "real")
end

function sanitize(tag, dct)
    if tag == "header"
        dct = sanitize_header(dct)
    elseif tag == "beta"
        dct = sanitize_beta(dct)
    elseif tag == "dij"
        dct = sanitize_dij(dct)
    elseif tag == "pswfc"
        dct = sanitize_pswfc(dct)
    elseif tag == "qij"
        dct = sanitize_qij(dct)
    elseif tag == "rinner"
        dct = sanitize_rinner(dct)
    elseif any(tag .== ["r", "rab", "local", "rhoatom", "nlcc", "qfcoef"])
        dct = sanitize_vector(dct)
    end
    # Remove "lines" key if it exists and is not the only key
    Set(keys(dct)) != Set(("lines",)) && pop!(dct, "lines", nothing)
    return dct
end

function upf1_block_to_dict(lines::Vector{String})
    dct = Dict{String,Any}()
    content_lines = filter(!isempty, strip.(lines))
    !isempty(content_lines) && (dct["lines"] = content_lines)
    # Loop over all (nested) children
    while true
        # Find the next opening tag
        inext = findfirst(startswith.(strip.(lines), "<"))
        # If there are no more tags, we are finished
        isnothing(inext) && break
        # Extract the tag name
        match = lines[inext]
        @assert !isempty(match) "Malformed UPF1 file"
        tag = split(strip(match)[(begin + 1):(end - 1)])[1]
        # Extract and parse that block
        block = extract_block(tag, lines)
        # Remove everything up to the end of the new block
        deleteat!(lines, 1:(inext + length(block) - 1))
        # Parse the contents of the block with a generic parser
        subdct = upf1_block_to_dict(block[(begin + 1):(end - 1)])
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

function upfv1contents_to_dict(str::String)
    dct = upf1_block_to_dict(string.(split(str, "\n")))
    pop!(dct, "lines", nothing)
    return dct
end
