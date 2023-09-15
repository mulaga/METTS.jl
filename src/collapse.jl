function collapse!(psi::MPS, basis::AbstractString=nothing)
    if basis === nothing || basis == "Z"
        return sample!(psi)
    elseif basis == "X"
        s = siteinds(psi)
        N = length(psi)
        Ry_gates = ops([("Ry", n, (θ=π / 2,)) for n in 1:N], s)
        psi = apply(Ry_gates, psi)
        return sample!(psi)
    else
        error("Invalid basis specified in METTS collapse.")
    end
end


function collapse_with_qn!(psi::MPS, basis::AbstractString=nothing)
    if basis === nothing || basis == "Z"
        return sample!(psi)
    elseif basis == "X"
        psi_noqn = removeqns(psi)
        s = siteinds(psi_noqn)
        N = length(psi_noqn)
        Ry_gates = ops([("Ry", n, (θ=π / 2,)) for n in 1:N], s)
        psi_noqn = apply(Ry_gates, psi_noqn)
        return sample!(psi_noqn)
    else
        error("Invalid basis specified in METTS collapse.")
    end
end
