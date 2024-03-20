module Geometry

    using LinearAlgebra, Distances, Interpolations, SparseArrays
    export compute_elements_length, compute_volume, centroids, weights, uniform_cross_sectional_areas_from_Vmax, compute_penalised_property, compute_angles

    function compute_elements_length(coord, connect; start=0)
        """Return a Vector with element lenghts
        """
        m = size(connect)[1]
        elements_length = zeros(m)

        for e in 1:m
            i = (1 - start) + connect[e, 1]
            j = (1 - start) + connect[e, 2]
            xi, yi, zi = coord[i, 1:3]
            xj, yj, zj = coord[j, 1:3]
            elements_length[e] = sqrt((xi - xj)^2 + (yi - yj)^2 + (zi - zj)^2)
        end

        return elements_length
    end

    function compute_angles(coord::Matrix{Float64}, connect::Matrix{Int64}, lengths::Vector{Float64})
        return (coord[connect[:, 2], :] .- coord[connect[:, 1], :]) ./ lengths
    end

    function compute_volume(A::Vector{Float64}, L::Vector{Float64})
        return A ⋅ L
    end

    function centroids(coord, connect, filter::Bool; start=0) #TODO: remove filter in this function and everywhere is called
        """Return a Matrix with centroid coordinates
        """
        m = size(connect)[1]
        centroids = zeros(m, 3)

        for e in 1:m
            i = (1 - start) + connect[e, 1]
            j = (1 - start) + connect[e, 2]

            centroids[e, :] = .5 * (coord[i, :] + coord[j, :])
        end

        return centroids
    end

    function weights(centroids::Matrix{Float64}, R::Float64, filter::Bool)
        if !filter
            return nothing
        end

        distance_matrix = pairwise(Euclidean(), centroids, dims=1)
        w = max.(0, R .- distance_matrix) # broadcast max (max.)
        return sparse(w)
    end

    function uniform_cross_sectional_areas_from_Vmax(vol_max, lengths)
        A0 = vol_max / sum(lengths)
        return A0
    end

    function compute_penalised_property(ρ::Vector{Float64}, p::Float64=1., prop_min::Float64=1e-10, prop_max::Float64=1.)
        return prop_min .+ ρ .^ p .* (prop_max - prop_min)
    end

    heaviside(t) = 0.5 .* (sign.(t) .+ 1.)

    function compute_penalised_property(ρ::Vector{Float64}, p::String, prop_min::Float64=1e-10, prop_max::Float64=1.)
        fun = heaviside(0.5 .- ρ) .* (-16ρ .^ 4 + 12ρ .^ 3) .+ heaviside(ρ .- 0.5) .* ρ
        return prop_min .+ fun .* (prop_max - prop_min)
    end

    function compute_penalised_property(ρ, p::Tuple, Amin::Float64, Amax::Float64)
        return Amin .+ compute_penalised_property.(ρ, [p]) .* (Amax - Amin)
    end

    function compute_penalised_property(ρ, interpolated_funs::Tuple{Interpolations.GriddedInterpolation, Interpolations.GriddedInterpolation})
        if ρ < 0.5
            f_interpolated = interpolated_funs[1]
            return f_interpolated(ρ)
        else
            return ρ
        end
    end

end # module
