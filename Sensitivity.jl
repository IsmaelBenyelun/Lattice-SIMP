module Sensitivity

    using FEM, Interpolations, LinearAlgebra, SparseArrays, Statistics, Base.Threads, ProgressBars, SuiteSparse
    export compliance_sensitivity_FOSM, compliance_sensitivity_MC
    export compliance_sensitivity, compliance_sensitivity_2d, compliance_sensitivity_3d
    export get_sensitivities, get_sensitivities_2d, get_sensitivities_3d, compute_∂J_∂x_2d, compute_∂J_∂x_3d
    export volume_sensitivity

    ### DSO functions
    function compliance_sensitivity(m, coord, connect, L, A, properties, u, ∂filter)
        return compute_∂J_∂x(m, coord, connect, L, A, properties, u, ∂filter)
    end

    function compliance_sensitivity_2d(m, coord, connect, L, A, properties, u, weights)
        return compute_∂J_∂x_2d(m, coord, connect, L, A, properties, u, weights, "area")
    end

    function compliance_sensitivity_2d(m, coord, connect, L, A::Float64, u, weights, ρ, p, Amin, Amax)
        return compute_∂J_∂x_2d(m, coord, connect, L, A, u, weights, "area", ρ, p, Amin, Amax)
    end

    function compliance_sensitivity_2d(m, coord, connect, L, A, properties, u, weights, ρ, p, Emin, Emax)
        return compute_∂J_∂x_2d(m, coord, connect, L, A, properties, u, weights, "area", ρ, p, Emin, Emax)
    end

    function compliance_sensitivity_3d(m, connect, L, angles, properties, u, weights, ρ, p, Amin, Amax)
        return compute_∂J_∂x_3d(m, connect, L, ones(m), angles, properties, u, weights, "area", ρ, p, Amin, Amax)
    end

    function volume_sensitivity(lengths::Array, w, l)
        """∂filter is the identity matrix if no filter is applied
        """
        ∂g_∂shat = lengths
        # ∂g_∂s = transpose(∂filter) * ∂g_∂shat
        ∂g_∂s = w * (∂g_∂shat ./ l) ./ (w * (1 ./ l))
        return ∂g_∂s[:]
    end

    function volume_sensitivity(ρ::Vector{Float64}, p, lengths::Vector{Float64}, weights, Amin::Float64, Amax::Float64; filter_type::String="lengths")
        ∂A_∂ρ = ∂A_∂s(ρ, p, Amin, Amax)
        ∂g_∂shat = lengths .* ∂A_∂ρ
        # ∂g_∂s = transpose(∂filter) * ∂g_∂shat
        if filter_type == "standard"
            ∂g_∂s = weights * (∂g_∂shat ./ sum(weights, dims=2))
        elseif filter_type == "lengths"
            ∂g_∂s = (weights * (∂g_∂shat ./ lengths)) ./ (weights * (1 ./ lengths))
        end

        return ∂g_∂s[:]
    end

    ### RSO -- FOSM function
    function compliance_sensitivity_FOSM(σJ, Cr, ∂J_∂x, ∂J_∂r, ∂2J_∂r∂x, α_S::Float64=.1, α_E::Float64=1.)
        # Mean derivative
        ∂μJ_∂x = ∂J_∂x
        # Variance derivative
        b = 2 * Cr * ∂J_∂r
        ∂σJ2_∂x = transpose(∂2J_∂r∂x) * b # We assume that ∂Cr/∂x = 0

        # println("Sensitivities: $(round(norm(∂μJ_∂x), digits=4)) ($(norm((1 / (2σJ) .* ∂σJ2_∂x))))")
        # Return
        return α_E .* ∂μJ_∂x + (α_S / (2σJ)) .* ∂σJ2_∂x#, ∂σJ2_∂x ./ (2σJ)
    end

    function compliance_sensitivity_FOSM(σJ, Cr, ∂J_∂x, ∂J_∂r, ∂2J_∂r∂x, weights, α_S::Float64=.1, α_E::Float64=1.)
        # Mean derivative
        ∂μJ_∂x = ∂J_∂x
        # Variance derivative
        b = 2 * Cr * ∂J_∂r
        ∂σJ2_∂xhat = transpose(∂2J_∂r∂x) * b # We assume that ∂Cr/∂x = 0

        # Filter variance derivative
        ∂σJ2_∂x = weights * (∂σJ2_∂xhat ./ sum(weights, dims=2))

        # println("Sensitivities: $(round(norm(∂μJ_∂x), digits=4)) ($(norm((1 / (2σJ) .* ∂σJ2_∂x))))")
        # Return
        return α_E * ∂μJ_∂x + (α_S / (2σJ)) * ∂σJ2_∂x#, ∂σJ2_∂x ./ (2σJ)
    end

    function compliance_sensitivity_FOSM(σJ, Qchol::SuiteSparse.CHOLMOD.Factor{Float64}, ∂J_∂x, ∂J_∂r, ∂2J_∂r∂x, weights, α_S::Float64=.1, α_E::Float64=1.)
        # Mean derivative
        ∂μJ_∂x = ∂J_∂x
        # Variance derivative
        b = 2 * (Qchol \ ∂J_∂r)
        ∂σJ2_∂xhat = transpose(∂2J_∂r∂x) * b # We assume that ∂Cr/∂x = 0

        # Filter variance derivative
        ∂σJ2_∂x = weights * (∂σJ2_∂xhat ./ sum(weights, dims=2))

        # println("Sensitivities: $(round(norm(∂μJ_∂x), digits=4)) ($(norm((1 / (2σJ) .* ∂σJ2_∂x))))")
        # Return
        return α_E * ∂μJ_∂x + (α_S / (2σJ)) * ∂σJ2_∂x#, ∂σJ2_∂x ./ (2σJ)
    end

    ### RSO -- MC function
    function compliance_sensitivity_MC(J_array, ∂J_∂x_array, μJ, σJ, α_S::Float64=.1, α_E::Float64=1.)
        n_MC = size(J_array, 1)
        # Mean derivative
        ∂μJ_∂x = mean(∂J_∂x_array, dims=2)
        # Std. derivative
        ∂σJ_∂x = ∂J_∂x_array * (J_array .- μJ) ./ (n_MC * σJ)

        # println("Sensitivities: $(round(norm(∂μJ_∂x), digits=4)) ($(norm(∂σJ_∂x)))")
        # Return
        return α_E * ∂μJ_∂x + α_S * ∂σJ_∂x#, ∂σJ_∂x
    end

    ### Compute all sensitivities
    function get_sensitivities(m, coord, connect, L, A, properties, u, ∂filter, Kglob, dof_f)
        ∂J_∂x = compute_∂J_∂x(m, coord, connect, L, A, properties, u, ∂filter)
        ∂J_∂r = ∂J_∂x
        ∂2J_∂r∂x = compute_∂2J_∂r∂x(m, coord, connect, L, A, properties, u, ∂filter, Kglob, dof_f)
        return ∂J_∂x, ∂J_∂r, ∂2J_∂r∂x
    end

    function get_sensitivities_2d(m, coord, connect, L, A, properties, u, weigths, Kfactorization, dof_f, uncertainty="area")
        ∂J_∂x = compute_∂J_∂x_2d(m, coord, connect, L, A, properties, u, weigths, "area")
        ∂J_∂r = uncertainty == "area" ? ∂J_∂x : compute_∂J_∂x_2d(m, coord, connect, L, A, properties, u, weigths, uncertainty)
        ∂2J_∂r∂x = compute_∂2J_∂r∂x_2d(m, coord, connect, L, A, properties, u, Kfactorization, dof_f, uncertainty)
        return ∂J_∂x, ∂J_∂r, ∂2J_∂r∂x
    end

    function get_sensitivities_2d(m, coord, connect, L, A, Ahat, properties, u, uhat, weights, Kfactorization, dof_f, uncertainty, ρ, ρhat, p, Amin, Amax)
        # ∂J_∂s = compute_∂J_∂x_2d(m, coord, connect, L, A, properties, u, weights, "area", ρ, p, Amin, Amax)
        ∂J_∂s = compute_∂J_∂x_2d(m, coord, connect, L, Ahat, properties, uhat, weights, "area", ρhat, p, Amin, Amax)
        # ∂J_∂r = uncertainty == "area" ? ∂J_∂s : compute_∂J_∂x_2d(m, coord, connect, L, A, properties, u, weights, uncertainty)
        ∂J_∂r = compute_∂J_∂x_2d(m, coord, connect, L, A, properties, u, weights, uncertainty, ρ, p, Amin, Amax)
        ∂2J_∂r∂shat = compute_∂2J_∂r∂x_2d(m, coord, connect, L, Ahat, properties, uhat, Kfactorization, dof_f, uncertainty, ρhat, ρhat, p, Amin, Amax)
        # ∂2J_∂r∂shat = compute_∂2J_∂r∂x_2d(m, coord, connect, L, Ahat, properties, uhat, Kfactorization_hat, dof_f, uncertainty, ρ, ρ, p, Amin, Amax)
        return ∂J_∂s, ∂J_∂r, ∂2J_∂r∂shat
    end

    function get_sensitivities_3d(m, coord, connect, L, A, Ahat, angles, properties, u, uhat, weights, Kfactorization, dof_f, uncertainty, ρ, ρhat, p, Amin, Amax)
        ∂J_∂s = compute_∂J_∂x_3d(m, connect, L, A, angles, properties, u, weights, "area", ρ, p, Amin, Amax)
        ∂J_∂r = compute_∂J_∂x_3d(m, connect, L, A, angles, properties, u, weights, uncertainty, ρ, p, Amin, Amax)
        print("Computing d2J_drds"); @time ∂2J_∂r∂shat = compute_∂2J_∂r∂x_3d(m, coord, connect, L, Ahat, angles, properties, uhat, Kfactorization, dof_f, uncertainty, ρhat, p, Amin, Amax)
        return ∂J_∂s, ∂J_∂r, ∂2J_∂r∂shat
    end

    ########################################################################################
    ### DERIVATIVES ∂J/∂(•)
    ########################################################################################

    ### First order derivatives
    function compute_∂J_∂x(m, coord, connect, L, A, properties, u, ∂filter)
        """∂filter is the identity matrix if no filter is applied
        """
        E = properties[1]
        ν = properties[2]

        ### Calculate W_elem
        sensitivities = zeros(m, 1)
        @threads for i in 1:m
            a = 1 + connect[i, 1]    # NODO INICIAL
            b = 1 + connect[i, 2]    # NODO FINAL
            c = 6*a - 5             # LOS VALORES DE LOS GDL VAN DESDE c HASTA d PARA EL PRIMER NODO
            d = 6*a
            e = 6*b - 5             # DE e HASTA f PARA EL SEGUNDO NODO
            f = 6*b

            dKloc_dA = ∂K_∂A(coord[a, 1:3], coord[b, 1:3], E, ν, A[i], L[i])
            ul = vcat(u[c:d], u[e:f])
            sensitivities[i] =  -ul ⋅ (dKloc_dA * ul)
        end

        out = transpose(vec(sensitivities)) * ∂filter
        return out[:]
    end

    function compute_∂J_∂x_2d(m, coord, connect, L, A, properties, u, weights, uncertainty)
        """Standard sensitivities w.r.t to either area or Young's modulus. Without penalisation
        """
        E = properties[1]

        ### Calculate W_elem
        sensitivities = zeros(m, 1)

        @threads for i in 1:m
            # Get unit vector along truss direction
            node1, node2 = connect[i, :]
            t = coord[node2, :] - coord[node1, :]
            e1a = t / norm(t)
            # Compute ∂K_∂A
            x = uncertainty == "area" ? E : A[i]
            dKloc_dA = ∂K_∂x_2d(e1a, x, L[i])
            # Get displacements affecting local nodes
            dof_local = [2*node1 - 1, 2*node1, 2*node2 - 1, 2*node2]
            ul = u[dof_local]
            sensitivities[i] =  -ul ⋅ (dKloc_dA * ul)
        end
        if uncertainty == "area"
            # out = transpose(vec(sensitivities)) * ∂filter
            out = (weights * (sensitivities ./ L)) ./ (weights * (1 ./ L))
        else
            out = sensitivities
        end
        return out[:]
    end

    function compute_∂J_∂x_2d(m, coord, connect, L, A, properties, u, weights, derivative, ρ, p, Amin, Amax)
        """Penalised areas
        """
        E = properties[1]

        ### Calculate W_elem
        sensitivities = zeros(m, 1)

        @threads for i in 1:m
            # Get unit vector along truss direction
            node1, node2 = connect[i, :]
            t = coord[node2, :] - coord[node1, :]
            e1a = t / norm(t)
            # Compute ∂K_∂A
            x = derivative == "area" ? E : A[i]
            dKloc_dA = ∂K_∂x_2d(e1a, x, L[i])
            dAi_dsi = derivative == "area" ? ∂A_∂s(ρ[i], p, Amin, Amax) : 1.
            dKloc_ds = dAi_dsi .* dKloc_dA
            # Get displacements affecting local nodes
            dof_local = [2*node1 - 1, 2*node1, 2*node2 - 1, 2*node2]
            ul = u[dof_local]
            sensitivities[i] =  -ul ⋅ (dKloc_ds * ul)
        end
        if derivative == "area"
            # out = transpose(vec(sensitivities)) * ∂filter
            out = (weights * (sensitivities ./ L)) ./ (weights * (1 ./ L))
        else
            out = sensitivities
        end
        return out[:]
    end

    function compute_∂J_∂x_2d(m, coord, connect, L, A::Float64, u, weights, derivative, ρ, p, Emin, Emax)
        """SIMP approach (A ≡ constant)
        """

        ### Calculate W_elem
        sensitivities = zeros(m, 1)

        @threads for i in 1:m
            # Get unit vector along truss direction
            node1, node2 = connect[i, :]
            t = coord[node2, :] - coord[node1, :]
            e1a = t / norm(t)
            # Compute ∂K_∂A
            dKloc_dE = ∂K_∂x_2d(e1a, A, L[i])
            dEi_dsi = ∂A_∂s(ρ[i], p, Emin, Emax)
            dKloc_ds = dEi_dsi .* dKloc_dE
            # Get displacements affecting local nodes
            dof_local = [2*node1 - 1, 2*node1, 2*node2 - 1, 2*node2]
            ul = u[dof_local]
            sensitivities[i] =  -ul ⋅ (dKloc_ds * ul)
        end
        if derivative == "area"
            # out = transpose(vec(sensitivities)) * ∂filter
            out = (weights * (sensitivities ./ L)) ./ (weights * (1 ./ L))
        else
            out = sensitivities
        end
        return out[:]
    end

    function compute_∂J_∂x_3d(m, connect, lengths, areas, angles, properties, u, weights, derivative, ρ, p, Amin, Amax)
        """Penalised areas
        """

        E = properties[1]
        sensitivities = zeros(m, 1)

        @threads for i in 1:m
            # Get unit vector along truss direction
            node1, node2 = connect[i, :]
            # Compute ∂K_∂A
            x = derivative == "area" ? E : areas[i]
            dKloc_dA = truss_3d(x, lengths[i], angles[i, 1], angles[i, 2], angles[i, 3])
            dAi_dsi = derivative == "area" ? ∂A_∂s(ρ[i], p, Amin, Amax) : 1.
            dKloc_ds = dAi_dsi .* dKloc_dA
            # Get displacements affecting local nodes
            dof_local = [3*node1 - 2, 3*node1 - 1, 3*node1, 3*node2 - 2, 3*node2 - 1, 3*node2]
            ul = u[dof_local]
            sensitivities[i] =  -ul ⋅ (dKloc_ds * ul)
        end
        if derivative == "area"
            out = weights * (sensitivities ./ sum(weights, dims=2))
        else
            out = sensitivities
        end
        return out[:]
    end

    # Second order derivatives
    function compute_∂2J_∂r∂x(m, coord, connect, L, A, properties, u, ∂filter, Kglob, dof_f)
        # Unpack
        E = properties[1]
        ν = properties[2]

        n = size(coord, 1)
        dof = 6*n

        # Loop to compute and save ∂K/∂x, ∂u/∂x
        ∂K_∂r_array = [zeros(0, 0) for _ in 1:m]
        ∂u_∂x_array = zeros(dof, m)

        print("Computing the inverse... ")
        # Calculate pseudo-inverse of Kr (only once)
        # invK = pinv(Matrix(Kglob[:, dof_f]))
        invK = inv(Matrix(Kglob[dof_f, dof_f]))

        println("Solving 'K du/dx = -dK/dx u' system...")
        @threads for i = 1:m
            # Get the correspondent element (dofs and unit vector)
            a = 1 + connect[i, 1]    # NODO INICIAL
            b = 1 + connect[i, 2]    # NODO FINAL
            c = 6*a - 5             # LOS VALORES DE LOS GDL VAN DESDE c HASTA d PARA EL PRIMER NODO
            d = 6*a
            e = 6*b - 5             # DE e HASTA f PARA EL SEGUNDO NODO
            f = 6*b

            dof_local = [c:d; e:f]
            # Load ∂K/∂Ai
            ∂K_∂Ai = ∂K_∂A(coord[a, 1:3], coord[b, 1:3], E, ν, A[i], L[i])
            # Solve system K (∂u/∂x)_i = -(∂K/∂x)_i u =: b
            b = zeros(dof)
            b[dof_local] = - ∂K_∂Ai * u[dof_local]
            ∂u_∂Ai = invK * b[dof_f]
            # Store ∂K/∂Ai and ∂u/∂Ai
            ∂K_∂r_array[i] = ∂K_∂Ai
            ∂u_∂x_array[dof_f, i] = ∂u_∂Ai
        end

        # Generate the (Hessian) matrix (threaded loop more efficient than vectorized e.g. Eigensum.jl)
        ∂2J_∂r∂x = zeros(m, m)
        @threads for i = 1:m
            # Get the correspondent element (dofs and unit vector)
            a = 1 + connect[i, 1]    # NODO INICIAL
            b = 1 + connect[i, 2]    # NODO FINAL
            c = 6*a - 5             # LOS VALORES DE LOS GDL VAN DESDE c HASTA d PARA EL PRIMER NODO
            d = 6*a
            e = 6*b - 5             # DE e HASTA f PARA EL SEGUNDO NODO
            f = 6*b

            dof_local = [c:d; e:f]
            # Compute ∂2K/∂x^2 (from element i)
            ∂2K_∂x2_i = ∂2K_∂A2(coord[a, 1:3], coord[b, 1:3], E, ν, L[i])

            ∂2J_∂r∂x[i, :] = -2 .* transpose(u[dof_local]) * ∂K_∂r_array[i] * ∂u_∂x_array[dof_local, :] # Vector
            ∂2J_∂r∂x[i, i] += -transpose(u[dof_local]) * ∂2K_∂x2_i * u[dof_local] # Scalar (∂2K/∂xi∂xj = 0 if i != j)
        end

        # Apply filter ∂2J/∂r∂s = (∂x/∂r)^T ∂2J/∂x2 (∂x/∂s)
        out = transpose(∂filter) * ∂2J_∂r∂x * ∂filter
        return out
    end

    function compute_∂2J_∂r∂x_2d(m, coord, connect, L, A, properties, u, Kfactorization, dof_f, derivative)
        # Unpack
        E = properties[1]
        n = size(coord, 1)
        n_dof = 2*n

        # Loop to compute and save ∂K/∂x, ∂u/∂x
        ∂K_∂r_array = [zeros(0, 0) for _ in 1:m]
        ∂u_∂x_array = zeros(n_dof, m)

        # Calculate inverse of Krr (only once)

        @threads for i = 1:m
            # Get the correspondent element (dofs and unit vector)
            node1, node2 = connect[i, :]
            t = coord[node2, :] - coord[node1, :]
            e1a = t / norm(t)
            dof_local = [2*node1 - 1, 2*node1, 2*node2 - 1, 2*node2]
            # Load ∂K_∂Ai
            ∂K_∂Ai = ∂K_∂x_2d(e1a, E, L[i])
            # Solve system K (∂u/∂x)_i = -(∂K/∂x)_i u =: b
            b = zeros(n_dof)
            b[dof_local] = - ∂K_∂Ai * u[dof_local]
            ∂u_∂Ai = Kfactorization \ b[dof_f]
            # ∂u_∂Ai = sparse(Kglob[dof_f, dof_f]) \ b[dof_f]
            ∂K_∂ri = derivative == "area" ? ∂K_∂Ai : ∂K_∂x_2d(e1a, A[i], L[i])
            # Store ∂K_∂ri and ∂u_∂Ai
            ∂K_∂r_array[i] = ∂K_∂ri
            ∂u_∂x_array[dof_f, i] = ∂u_∂Ai
        end

        # Generate the (Hessian) matrix
        ∂2J_∂r∂x = zeros(m, m)

        for i = 1:m
            node1, node2 = connect[i, :]
            t = coord[node2, :] - coord[node1, :]
            e1a = t / norm(t)
            dof_local = [2*node1 - 1, 2*node1, 2*node2 - 1, 2*node2]
            ∂2J_∂ri∂x = -2 .* transpose(u[dof_local]) * ∂K_∂r_array[i] * ∂u_∂x_array[dof_local, :]
            # Compute ∂2K_∂r∂x (diagonal)
            ∂2K_∂ri∂xi = derivative == "area" ? zeros(4, 4) : ∂2K_∂E∂A_2d(e1a, L[i])
            # Assign (row-wise) and diagonal term
            ∂2J_∂r∂x[i, :] = ∂2J_∂ri∂x
            ∂2J_∂r∂x[i, i] += -transpose(u[dof_local]) * ∂2K_∂ri∂xi * u[dof_local] # Scalar (∂2K/∂xi∂xj = 0 if i != j)
        end

        # out = ∂2J_∂r∂x * ∂filter
        out = ∂2J_∂r∂x
        return out
    end

    function compute_∂2J_∂r∂x_2d(m, coord, connect, L, A, properties, u, Kfactorization, dof_f, derivative, ρ, ρhat, p, Amin, Amax)
        # Unpack
        E = properties[1]
        n = size(coord, 1)
        n_dof = 2*n

        # Loop to compute and save ∂K/∂x, ∂u/∂x
        ∂K_∂r_array = [zeros(0, 0) for _ in 1:m]
        ∂u_∂x_array = zeros(n_dof, m)

        # Calculate inverse of Krr (only once)
        # invK = inv(Matrix(Kglob[dof_f, dof_f]))

        @threads for i = 1:m
            # Get the correspondent element (dofs and unit vector)
            node1, node2 = connect[i, :]
            t = coord[node2, :] - coord[node1, :]
            e1a = t / norm(t)
            dof_local = [2*node1 - 1, 2*node1, 2*node2 - 1, 2*node2]
            # Load ∂K_∂Ai
            ∂K_∂Ai = ∂K_∂x_2d(e1a, E, L[i])
            ∂Ai_∂si = ∂A_∂s(ρ[i], p, Amin, Amax)
            # ∂Ai_∂si = ∂A_∂s(ρhat[i], p, Amin, Amax)
            ∂K_∂si = ∂Ai_∂si .* ∂K_∂Ai
            # Solve system K (∂u/∂x)_i = -(∂K/∂x)_i u =: b
            b = zeros(n_dof)
            b[dof_local] = - ∂K_∂si * u[dof_local]
            # b[dof_local] = - ∂K_∂si * uhat[dof_local]
            # ∂u_∂si = invK * b[dof_f]
            ∂u_∂si = Kfactorization \ b[dof_f]
            # ∂u_∂si = sparse(Kglob[dof_f, dof_f]) \ b[dof_f]
            ∂K_∂ri = derivative == "area" ? ∂K_∂si : ∂K_∂x_2d(e1a, A[i], L[i])
            # ∂K_∂ri = derivative == "area" ? ∂K_∂si : ∂K_∂x_2d(e1a, Ahat[i], L[i])
            # Store ∂K_∂ri and ∂u_∂si
            ∂K_∂r_array[i] = ∂K_∂ri
            ∂u_∂x_array[dof_f, i] = ∂u_∂si
        end

        # Generate the (Hessian) matrix
        ∂2J_∂r∂x = zeros(m, m)

        @threads for i = 1:m
            node1, node2 = connect[i, :]
            t = coord[node2, :] - coord[node1, :]
            e1a = t / norm(t)
            dof_local = [2*node1 - 1, 2*node1, 2*node2 - 1, 2*node2]
            ∂2J_∂ri∂x = -2 .* transpose(u[dof_local]) * ∂K_∂r_array[i] * ∂u_∂x_array[dof_local, :]
            # ∂2J_∂ri∂x = -2 .* transpose(uhat[dof_local]) * ∂K_∂r_array[i] * ∂u_∂x_array[dof_local, :]

            # Compute ∂2K_∂r∂x (diagonal)
            ∂2K_∂ri∂Ai = derivative == "area" ? zeros(4, 4) : ∂2K_∂E∂A_2d(e1a, L[i])
            # ∂Ai_∂si = ∂A_∂s(ρ[i], p, Amin, Amax)
            ∂Ai_∂si = ∂A_∂s(ρhat[i], p, Amin, Amax)
            ∂2K_∂ri∂si = ∂Ai_∂si .* ∂2K_∂ri∂Ai
            # Assign (row-wise) and diagonal term
            ∂2J_∂r∂x[i, :] = ∂2J_∂ri∂x
            ∂2J_∂r∂x[i, i] += -transpose(u[dof_local]) * ∂2K_∂ri∂si * u[dof_local] # Scalar (∂2K/∂xi∂xj = 0 if i != j)
            # ∂2J_∂r∂x[i, i] += -transpose(uhat[dof_local]) * ∂2K_∂ri∂si * uhat[dof_local] # Scalar (∂2K/∂xi∂xj = 0 if i != j)
        end

        # Apply filter ∂2J/∂r∂s = ∂2J/∂r∂shat (∂shat/∂s)
        # out = ∂2J_∂r∂x * ∂filter
        out = ∂2J_∂r∂x
        return out
    end

    function compute_∂2J_∂r∂x_3d(m, coord, connect, L, A, angles, properties, u, Kfactorization, dof_f, derivative, ρ, p, Amin, Amax)
        # Unpack
        E = properties[1]
        n = size(coord, 1)
        n_dof = 3*n

        # Loop to compute and save ∂K/∂x, ∂u/∂x
        ∂K_∂r_array = [zeros(0, 0) for _ in 1:m]
        ∂u_∂x_array = zeros(n_dof, m)

        @threads for i = ProgressBar(1:m)
            # Get the correspondent element (dofs and unit vector)
            node1, node2 = connect[i, :]
            dof_local = [3*node1 - 2, 3*node1 - 1, 3*node1, 3*node2 - 2, 3*node2 - 1, 3*node2]
            # Load ∂K_∂Ai
            ∂K_∂Ai = truss_3d(E, L[i], angles[i, 1], angles[i, 2], angles[i, 3])
            ∂Ai_∂si = ∂A_∂s(ρ[i], p, Amin, Amax)
            ∂K_∂si = ∂Ai_∂si .* ∂K_∂Ai
            # Solve system K (∂u/∂x)_i = -(∂K/∂x)_i u =: b
            b = zeros(n_dof)
            b[dof_local] = - ∂K_∂si * u[dof_local]
            ∂u_∂si = Kfactorization \ b[dof_f]
            # ∂u_∂si = sparse(Kglob[dof_f, dof_f]) \ b[dof_f]
            ∂K_∂ri = derivative == "area" ? ∂K_∂si : truss_3d(A[i], L[i], angles[i, 1], angles[i, 2], angles[i, 3])
            # Store ∂K_∂ri and ∂u_∂si
            ∂K_∂r_array[i] = ∂K_∂ri
            ∂u_∂x_array[dof_f, i] = ∂u_∂si
        end

        # Generate the (Hessian) matrix
        ∂2J_∂r∂x = zeros(m, m)

        @threads for i = ProgressBar(1:m)
            node1, node2 = connect[i, :]
            dof_local = [3*node1 - 2, 3*node1 - 1, 3*node1, 3*node2 - 2, 3*node2 - 1, 3*node2]
            ∂2J_∂ri∂x = -2 .* transpose(u[dof_local]) * ∂K_∂r_array[i] * ∂u_∂x_array[dof_local, :]

            # Compute ∂2K_∂r∂x (diagonal)
            ∂2K_∂ri∂Ai = derivative == "area" ? zeros(6, 6) : truss_3d(1., L[i], angles[i, 1], angles[i, 2], angles[i, 3])
            ∂Ai_∂si = ∂A_∂s(ρ[i], p, Amin, Amax)
            ∂2K_∂ri∂si = ∂Ai_∂si .* ∂2K_∂ri∂Ai
            # Assign (row-wise) and diagonal term
            ∂2J_∂r∂x[i, :] = ∂2J_∂ri∂x
            ∂2J_∂r∂x[i, i] += -transpose(u[dof_local]) * ∂2K_∂ri∂si * u[dof_local] # Scalar (∂2K/∂xi∂xj = 0 if i != j)
        end

        out = ∂2J_∂r∂x
        return out
    end

    ########################################################################################
    ### Stiffness Matrix derivatives
    ########################################################################################
    function ∂K_∂A(coordi::Vector{Float64}, coordj::Vector{Float64}, E::Float64, ν::Float64, A::Float64, L::Float64)

        G = E / (2*(1 + ν))
        EA = E * A

        dKlocal_dA = [
            [E/L 0 0 0 0 0 -E/L 0 0 0 0 0];
            [0 6EA/(π*L^3) 0 0 0 3EA / (π * L^2) 0 -6EA/(π*L^3) 0 0 0 3EA / (π * L^2)];
            [0 0 6EA/(π*L^3) 0 -3EA / (π * L^2) 0 0 0 -6EA/(π*L^3) 0 -3EA / (π * L^2) 0];
            [0 0 0 G*A/(π * L) 0 0 0 0 0 -G*A/(π * L) 0 0];
            [0 0 -3EA / (π * L^2) 0 2EA / (π * L) 0 0 0 3EA / (π * L^2) 0 EA / (π * L) 0];
            [0 3EA / (π * L^2) 0 0 0 2EA / (π * L) 0 -3EA / (π * L^2) 0 0 0 EA / (π * L)];
            [-E/L 0 0 0 0 0 E/L 0 0 0 0 0];
            [0 -6EA/(π*L^3) 0 0 0 -3EA / (π * L^2) 0 6EA/(π*L^3) 0 0 0 -3EA / (π * L^2)];
            [0 0 -6EA/(π*L^3) 0 3EA / (π * L^2) 0 0 0 6EA/(π*L^3) 0 3EA / (π * L^2) 0];
            [0 0 0 -G*A/(π * L) 0 0 0 0 0 G*A/(π * L) 0 0];
            [0 0 -3EA / (π * L^2) 0 EA / (π * L) 0 0 0 3EA / (π * L^2) 0 2EA / (π * L) 0];
            [0 3EA / (π * L^2) 0 0 0 EA / (π * L) 0 -3EA / (π * L^2) 0 0 0 2EA / (π * L)]
        ]

        u = [(coordj[1]-coordi[1]), (coordj[2]-coordi[2]), (coordj[3]-coordi[3])]/L    #ES MI EJE X DE LAS COORDENADAS LOCALES

        ang_z = asin(u[2])
        ang_y = atan(-u[3], u[1])
        ang_x = 0

        cos_x = cos(ang_x)
        cos_y = cos(ang_y)
        cos_z = cos(ang_z)

        sin_x = sin(ang_x)
        sin_y = sin(ang_y)
        sin_z = sin(ang_z)

        T = [
            [(cos_y*cos_z) (sin_y*sin_x - cos_y*sin_z*cos_x) (cos_y*sin_z*sin_x + sin_y*cos_x)];
            [sin_z (cos_z*cos_x) (-cos_z*sin_x)];
            [(-sin_y*cos_z) (sin_y*sin_z*cos_x + cos_y*sin_x) (cos_y*cos_x - sin_y*sin_z*sin_x)]
        ]


        R = zeros(12,12)
        R[1:3,1:3] = copy(T)
        R[4:6,4:6] = copy(T)
        R[7:9,7:9] = copy(T)
        R[10:12,10:12] = copy(T)

        dKloc_dA = R * dKlocal_dA * transpose(R)

        return dKloc_dA
    end

    function ∂2K_∂A2(coordi::Vector{Float64}, coordj::Vector{Float64}, E::Float64, ν::Float64, L::Float64)

        G = E / (2*(1 + ν))

        d2Klocal_dA2 = [
            [0 0 0 0 0 0 0 0 0 0 0 0];
            [0 6E/(π*L^3) 0 0 0 3E / (π * L^2) 0 -6E/(π*L^3) 0 0 0 3E / (π * L^2)];
            [0 0 6E/(π*L^3) 0 -3E / (π * L^2) 0 0 0 -6E/(π*L^3) 0 -3E / (π * L^2) 0];
            [0 0 0 G/(π * L) 0 0 0 0 0 -G/(π * L) 0 0];
            [0 0 -3E / (π * L^2) 0 2E / (π * L) 0 0 0 3E / (π * L^2) 0 E / (π * L) 0];
            [0 3E / (π * L^2) 0 0 0 2E / (π * L) 0 -3E / (π * L^2) 0 0 0 E / (π * L)];
            [0 0 0 0 0 0 0 0 0 0 0 0];
            [0 -6E/(π*L^3) 0 0 0 -3E / (π * L^2) 0 6E/(π*L^3) 0 0 0 -3E / (π * L^2)];
            [0 0 -6E/(π*L^3) 0 3E / (π * L^2) 0 0 0 6E/(π*L^3) 0 3E / (π * L^2) 0];
            [0 0 0 -G/(π * L) 0 0 0 0 0 G/(π * L) 0 0];
            [0 0 -3E / (π * L^2) 0 E / (π * L) 0 0 0 3E / (π * L^2) 0 2E / (π * L) 0];
            [0 3E / (π * L^2) 0 0 0 E / (π * L) 0 -3E / (π * L^2) 0 0 0 2E / (π * L)]
        ]

        u = [(coordj[1]-coordi[1]), (coordj[2]-coordi[2]), (coordj[3]-coordi[3])]/L    #ES MI EJE X DE LAS COORDENADAS LOCALES

        ang_z = asin(u[2])
        ang_y = atan(-u[3], u[1])
        ang_x = 0

        cos_x = cos(ang_x)
        cos_y = cos(ang_y)
        cos_z = cos(ang_z)

        sin_x = sin(ang_x)
        sin_y = sin(ang_y)
        sin_z = sin(ang_z)

        T = [
            [(cos_y*cos_z) (sin_y*sin_x - cos_y*sin_z*cos_x) (cos_y*sin_z*sin_x + sin_y*cos_x)];
            [sin_z (cos_z*cos_x) (-cos_z*sin_x)];
            [(-sin_y*cos_z) (sin_y*sin_z*cos_x + cos_y*sin_x) (cos_y*cos_x - sin_y*sin_z*sin_x)]
        ]


        R = zeros(12,12)
        R[1:3,1:3] = copy(T)
        R[4:6,4:6] = copy(T)
        R[7:9,7:9] = copy(T)
        R[10:12,10:12] = copy(T)

        d2Kloc_dA2 = R * d2Klocal_dA2 * transpose(R)

        return d2Kloc_dA2
    end

    # function ∂K_∂A_2d(e1a::Vector{Float64}, E::Float64, L::Float64)
    #     # Declare unit vectors
    #     e1 = [1, 0, 0]
    #     e2 = [0, 1, 0]
    #     e3 = [0, 0, 1]
    #     # Unit vector local y' axis
    #     e2a = e3 × e1a

    #     dKlocal_dA = E / L * [
    #          1 0 -1 0;
    #          0 0  0 0;
    #         -1 0  1 0;
    #          0 0  0 0
    #     ]
    #     # Rotation matrix
    #     T = [e1⋅e1a e2⋅e1a; e1⋅e2a e2⋅e2a]
    #     # Q matrix
    #     Q = [T zeros(2, 2); zeros(2, 2) T]
    #     # Rotate local stiffness matrix
    #     dKloc_dA = transpose(Q) * dKlocal_dA * Q
    #     return dKloc_dA
    # end

    function ∂K_∂x_2d(e1a::Vector{Float64}, x::Float64, L::Float64)
        # Declare unit vectors
        e1 = [1, 0, 0]
        e2 = [0, 1, 0]
        e3 = [0, 0, 1]
        # Unit vector local y' axis
        e2a = e3 × e1a

        dKlocal_dA = x / L * [
             1 0 -1 0;
             0 0  0 0;
            -1 0  1 0;
             0 0  0 0
        ]
        # Rotation matrix
        T = [e1⋅e1a e2⋅e1a; e1⋅e2a e2⋅e2a]
        # Q matrix
        Q = [T zeros(2, 2); zeros(2, 2) T]
        # Rotate local stiffness matrix
        dKloc_dA = transpose(Q) * dKlocal_dA * Q
        return dKloc_dA
    end

    function ∂2K_∂E∂A_2d(e1a::Vector{Float64}, L::Float64)
        return ∂K_∂x_2d(e1a, 1., L)
    end

    ########################################################################################
    ### Area derivative (penalised densities method)
    ########################################################################################
    function ∂A_∂s(ρ, p::Float64, Amin, Amax)
        return p .* (ρ .^ (p - 1)) .* (Amax - Amin)
    end

    heaviside(t) = 0.5 .* (sign.(t) .+ 1.)

    function ∂A_∂s(ρ, p::String, Amin, Amax)
        dfun = heaviside(0.5 .- ρ) .* (-64ρ .^ 3 + 36ρ .^ 2) .+ heaviside(ρ .- 0.5)
        return dfun .* (Amax - Amin)
    end

    function ∂A_∂s(ρ, p::Tuple, Amin::Float64, Amax::Float64)
        return ∂A_∂s.(ρ, [p]) .* (Amax - Amin)
    end

    function ∂A_∂s(ρ, interpolated_funs::Tuple{Interpolations.GriddedInterpolation, Interpolations.GriddedInterpolation})
        if ρ < 0.5
            df_interpolated = interpolated_funs[2]
            return df_interpolated(ρ)
        else
            return 1.
        end
    end

end # module
