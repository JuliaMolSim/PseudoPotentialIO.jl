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

    # If the element has no children, extract the contents of nodecontent(element)
    if countelements(element) == 0  # TODO: countelements or countnodes?
        text = strip(nodecontent(element))
        # Return early if there is no text
        isempty(text) && return result
        value = try
            parse.(Float64, split(text))
        catch
            text
        end
        if isempty(result)
            result = value
        else
            result["content"] = value
        end
        return result
    end

    # If the element has children, extract the contents of each child
    for child in eachelement(element)
        child_result = upf2_block_to_dict(child)
        tag = child.name |> x -> split(x, ".") |> first |> x -> replace(x, "PP_" => "") |> lowercase
        if haskey(result, tag)
            isa(result[tag], Vector) || (result[tag] = [result[tag]])
            push!(result[tag], child_result)
        else
            result[tag] = child_result
        end
    end

    return result
end

function upfv2contents_to_dict(str::String)
    doc = EzXML.parsexml(str)
    root = EzXML.root(doc)
    dct = upf2_block_to_dict(root)
    pop!(dct, "upf_version", nothing)  # Remove UPF version key
    return dct
end
