"""
Compute the value of a linear mesh `a * i + b` at index `i`.
"""
linear_mesh(i::Int, a, b) = a * i + b

"""
Compute the value of a logarithmic mesh `b * exp(a * (i - 1))` at index `i`.
NB: for UPF pseudopotentials `a = dx`, `b = e^{xmin} / zmesh`.
"""
logarithmic_mesh1(i::Int, a, b) = b * exp(a * (i - 1))
function logarithmic_mesh1(i::Int, xmin, dx, z)
    return exp(xmin) * exp((i - 1) * dx) / z
end

"""
Compute the value of a logarithmic mesh `b * (exp(a * (i - 1)) - 1)` at index `i`.
NB: for UPF pseudopotentials `a = dx`, `b = e^{xmin} / zmesh`.
"""
logarithmic_mesh2(i::Int, a, b) = b * (exp((i - 1) * a) - 1)
function logarithmic_mesh2(i::Int, xmin, dx, z)
    return exp(xmin) * (exp((i - 1) * dx) - 1) / z
end

"""
Guess whether a numerical mesh is linear or one of two types of logarithmic mesh used in 
UPF pseudopotentials.
"""
function guess_mesh_type(r::AbstractVector, rab::AbstractVector)
    nr = length(r)
    # Try linear
    a = r[2] - r[1]
    b = r[1] - a
    rguess = linear_mesh.(1:nr, a, b)
    if all(rguess .≈ r) & all(round.(rab, digits=4) .≈ a)
        # Improve the values
        a = mean(diff(r))
        b = r[1] - a
        return ("linear", a, b)
    end
    # Try log1
    a = log(r[2] / r[1])  # dx
    b = r[2] / exp(a)  # exp(xmin) / zmesh
    rguess = logarithmic_mesh1.(1:nr, a, b)
    if all(rguess .≈ r) && all(rab .≈ a .* r)
        return ("log_1", a, b)
    end
    # Try log2
    b = (r[2]^2 - r[3] * r[1]) / (r[1] + r[3] - 2 * r[2])
    a = log((r[2] + b) / (r[1] + b))
    rguess = logarithmic_mesh2.(1:nr, a, b)
    if all(rguess .≈ r) && all(rab .≈ a .* r .+ a * b)
        return ("log_2", a, b)
    end
    return ("unknown", NaN, NaN)
end
