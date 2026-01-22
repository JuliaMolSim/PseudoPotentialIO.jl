function sanitize(value::String)
    if (lowercase(strip(value)) in ("true", ".true.", "t"))
        return true
    elseif (lowercase(strip(value)) in ("false", ".false.", "f"))
        return false
    elseif tryparse(Int, value) !== nothing
        return parse(Int, value)
    elseif tryparse(Float64, value) !== nothing
        return parse(Float64, value)
    else
        return value
    end
end

function upf2_block_to_dict(element::EzXML.Node)
    result = Dict{String, Any}(
        lowercase(attr.name) => sanitize(nodecontent(attr))
        for attr in filter(x -> !(x in ("type", "columns", "size")), attributes(element))
    )
    # Manually adding principal quantum number if it is missing (e.g. result["label"]="3S" => result["n"]=3)
    if startswith(element.name, "PP_CHI") && !haskey(result, "n") && !haskey(result, "label")
        result["n"] = parse(Int, result["label"][begin])
    end
    # If the element has no children, extract the contents of the element
    if countelements(element) == 0
        text = strip(nodecontent(element))
        # Return early if there is no text
        isempty(text) && return result
        # Try to parse the content as a floating-point array
        # Otherwise, keep the content as a string
        value = tryparse.(Float64, split(text))
        if any(isnothing, value)
            value = text
        else
            value = Dict{String,Any}(
                "content" => value,
                "type" => "real",
                "size" => length(value),
                "columns" => length(split(first(split(text, "\n")))),
            )
        end
        isempty(result) && return value
        result["content"] = value
        return result
    end
    # If the element has children, extract the contents of each child
    for child in eachelement(element)
        child_result = upf2_block_to_dict(child)
        # Remove ".{index}" suffix and "PP_" prefix; lowercase the tag name
        tag = child.name |> x -> split(x, ".") |> first |> x -> replace(x, "PP_" => "") |> lowercase
        if haskey(result, tag)
            # Add an outer array if this is the second occurrence of this tag (e.g. BETA, CHI)
            isa(result[tag], Vector{typeof(child_result)}) || (result[tag] = [result[tag]])
            push!(result[tag], child_result)
        else
            result[tag] = child_result
        end
    end
    return result
end

function upfv2contents_to_dict(filecontents::AbstractString)
    root = EzXML.root(EzXML.parsexml(filecontents))
    dct = upf2_block_to_dict(root)
    pop!(dct, "version", nothing)
    return dct
end
