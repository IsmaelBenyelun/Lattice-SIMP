module Objectives

    using Sensitivity, Statistics, SuiteSparse
    export compute_J_FOSM, compute_σJ_FOSM, compute_J_MC

    ### FOSM
    function compute_J_FOSM(J, ∂J_∂r, Cr, α_S::Float64=0.1, α_E::Float64=1., return_μ_and_σ::Bool=true)
        # Get the mean
        μJ = compute_μJ_FOSM(J)
        # Compute the standard deviation
        σJ = compute_σJ_FOSM(∂J_∂r, Cr)

        # print("Mean and std: $(round(μJ, digits=4)) ($(round(σJ, digits=6))). ")
        # Return
        if return_μ_and_σ
            return α_E * μJ + α_S * σJ, μJ, σJ
        else
            return α_E * μJ + α_S * σJ
        end
    end

    function compute_σJ_FOSM(∂J_∂r, Cr)
        # Compute standard deviation of the objective
        b = Cr * ∂J_∂r
        return sqrt(transpose(∂J_∂r) * b)
    end

    function compute_σJ_FOSM(∂J_∂r, Qchol::SuiteSparse.CHOLMOD.Factor{Float64})
        # Compute standard deviation of the objective
        b = Qchol \ ∂J_∂r
        return sqrt(transpose(∂J_∂r) * b)
    end

    function compute_μJ_FOSM(J)
        # Compute μ_J for the FOSM method
        return J
    end

    ### MC
    function compute_J_MC(J_array, α_S::Float64=.1, α_E::Float64=1., return_μ_and_σ::Bool=true)
        # Get the mean
        μJ = compute_μJ_MC(J_array)
        # Compute the standard deviation
        σJ = compute_σJ_MC(J_array, μJ)

        # print("Mean and std: $(round(μJ, digits=4)) ($(round(σJ, digits=6))). ")
        # Return
        if return_μ_and_σ
            return α_E * μJ + α_S * σJ, μJ, σJ
        else
            return α_E * μJ + α_S * σJ
        end
    end

    function compute_σJ_MC(J_array, μJ)
        # Compute standard deviation of the objective (MC)
        return sqrt(mean((J_array .- μJ) .^ 2))
    end

    function compute_μJ_MC(J_array)
        # Compute μ_J for the MC method
        return mean(J_array)
    end

end #module
