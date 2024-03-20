module FEM

    using LinearAlgebra, Statistics, SparseArrays, LinearSolve, Geometry, Base.Threads

    export _solver, solver, solver_2d, build_Kglobal, solver_3d, truss_3d, get_member_forces_2d, get_member_forces_3d

    function matriz_rigidez(coordi::Vector{Float64}, coordj::Vector{Float64}, E::Float64, ν::Float64, A::Float64, L::Float64)
        I = [A^2 /(4π) A^2 /(4π)]
        J = A^2 / (2π)

        G = E / (2*(1 + ν))

        Klocal = [
            [E*A/L 0 0 0 0 0 -E*A/L 0 0 0 0 0];
            [0 12*E*I[2]/L^3 0 0 0 6*E*I[2]/L^2 0 -12*E*I[2]/L^3 0 0 0 6*E*I[2]/L^2];
            [0 0 12*E*I[1]/L^3 0 -6*E*I[1]/L^2 0 0 0 -12*E*I[1]/L^3 0 -6*E*I[1]/L^2 0];
            [0 0 0 G*J/L 0 0 0 0 0 -G*J/L 0 0];
            [0 0 -6*E*I[1]/L^2 0 4*E*I[1]/L 0 0 0 6*E*I[1]/L^2 0 2*E*I[1]/L 0];
            [0 6*E*I[2]/L^2 0 0 0 4*E*I[2]/L 0 -6*E*I[2]/L^2 0 0 0 2*E*I[2]/L];
            [-E*A/L 0 0 0 0 0 E*A/L 0 0 0 0 0];
            [0 -12*E*I[2]/L^3 0 0 0 -6*E*I[2]/L^2 0 12*E*I[2]/L^3 0 0 0 -6*E*I[2]/L^2];
            [0 0 -12*E*I[1]/L^3 0 6*E*I[1]/L^2 0 0 0 12*E*I[1]/L^3 0 6*E*I[1]/L^2 0];
            [0 0 0 -G*J/L 0 0 0 0 0 G*J/L 0 0];
            [0 0 -6*E*I[1]/L^2 0 2*E*I[1]/L 0 0 0 6*E*I[1]/L^2 0 4*E*I[1]/L 0];
            [0 6*E*I[2]/L^2 0 0 0 2*E*I[2]/L 0 -6*E*I[2]/L^2 0 0 0 4*E*I[2]/L]
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

        Kloc = R*Klocal*transpose(R)

        return Kloc
    end

    function giro(coordi::Vector{Float64}, coordj::Vector{Float64}, desplazamientos::Array{Float64}, fuerzas::Array{Float64})

        L = sqrt((coordj[1]-coordi[1])^2 + (coordj[2]-coordi[2])^2 + (coordj[3]-coordi[3])^2)

        u = [(coordj[1]-coordi[1]), (coordj[2]-coordi[2]), (coordj[3]-coordi[3])]/L

        ang_z = asin(u[2])
        ang_y = atan(-u[3],u[1])
        ang_x = 0

        cos_x = cos(ang_x)
        cos_y = cos(ang_y)
        cos_z = cos(ang_z)

        sin_x = sin(ang_x)
        sin_y = sin(ang_y)
        sin_z = sin(ang_z)

        T = [[(cos_y*cos_z) (sin_y*sin_x - cos_y*sin_z*cos_x) (cos_y*sin_z*sin_x + sin_y*cos_x)];
            [sin_z (cos_z*cos_x) (-cos_z*sin_x)];
            [(-sin_y*cos_z) (sin_y*sin_z*cos_x + cos_y*sin_x) (cos_y*cos_x - sin_y*sin_z*sin_x)]]

        R = zeros(12,12)
        R[1:3,1:3] = copy(T)
        R[4:6,4:6] = copy(T)
        R[7:9,7:9] = copy(T)
        R[10:12,10:12] = copy(T)

        desplazamientos_girados = transpose(R)*transpose(desplazamientos)
        fuerzas_giradas = transpose(R)*transpose(fuerzas)

        return fuerzas_giradas, desplazamientos_girados # DEVUELVE LA MATRIZ DE RIGIDEZ Y EL VECTOR DESPLAZAMIENTOS EN COORDENADAS LOCALES DEL ELEMENTO
    end

    function build_Kglobal(gdl::Int64, m::Int64, coord::Matrix{Float64}, conect::Matrix{Int64}, propiedades::Array{Float64}, L::Vector{Float64}, A::Vector{Float64})

        Kglobal = zeros(gdl, gdl)

        E = propiedades[1]
        ν = propiedades[2]

        for i in 1:m
            a = 1 + conect[i,1]     #NODO INICIAL
            b = 1 + conect[i,2]     #NODO FINAL
            c = 6*a - 5             #LOS VALORES DE LOS GDL VAN DESDE c HASTA d PARA EL PRIMER NODO
            d = 6*a
            e = 6*b - 5             #DE e HASTA f PARA EL SEGUNDO NODO
            f = 6*b

            Kl = matriz_rigidez(coord[a,1:3], coord[b,1:3], E, ν, A[i], L[i])
            #SE INTRODUCE CADA TERMINO DE LA MATRIZ DEL ELEMENTO EN LA MATRIZ DE LA ESTRUCTURA EN LA POSICION CORRESPONDIENTE
            Kglobal[c:d,c:d] += copy(Kl[1:6,1:6])
            Kglobal[c:d,e:f] += copy(Kl[1:6,7:12])
            Kglobal[e:f,c:d] += copy(Kl[7:12,1:6])
            Kglobal[e:f,e:f] += copy(Kl[7:12,7:12])
        end

        return Kglobal
    end

    function _solver(gdl::Int64, n::Int64, m::Int64, conect::Matrix{Int64}, coord::Matrix{Float64}, propiedades::Array{Float64}, index::Array{Any}, desp::Array{Float64}, force::Array{Float64}, Kglobal::Array{Float64}, A::Vector{Float64}, longitud_elementos::Array{Float64})
        ##########################################
        ################# SOLVER #################
        ##########################################
        displacement = copy(desp)

        dof_r = index
        dof_f = setdiff(1:gdl,dof_r)  # SON LOS GRADOS DE LIBERTAD LIBRES

        #REDUCCION DEL SISTEMA Y RESOLUCION
        f_f = force[dof_f]
        K_fr = Kglobal[dof_f, dof_r]  #TODAS LAS FILAS Y SOLO LAS COLUMNAS DE LAS CONDICIONES DE CONTORNO
        u_r = displacement[dof_r]
        K_ff = Kglobal[dof_f,dof_f]

        # Condensación estática
        u_f = sparse(K_ff) \ (f_f - K_fr*u_r)

        displacement[dof_f] = u_f

        energy = 0.5 * force ⋅ displacement

        #SACAMOS LOS DESPLAZAMIENTOS LINEALES Y ANGULARES POR SEPARADO PARA LEERLO MAS FACIL
        desp_lin = zeros(n, 3)
        desp_ang = zeros(n, 3)
        for i in 1:n
            desp_lin[i, :] = displacement[(6*i - 5):(6*i - 3)]
            desp_ang[i, :] = displacement[(6*i - 2):6*i]
        end

        #HACEMOS LOS MISMO CON FUERZAS Y MOMENTOS
        load = zeros(n, 3)
        momento = zeros(n, 3)

        for i in 1:n
            inicial_f = 6*i - 5
            final_f = 6*i -3
            load[i, :] = force[inicial_f:final_f]

            inicial_m = 6*i - 2
            final_m = 6*i
            momento[i, :] = force[inicial_m:final_m]
        end

        #########################################
        ######  DEFORMACIONES Y TENSIONES  ######
        #########################################

        coordn1 = zeros(1,3)
        coordn2 = zeros(1,3)
        desp1 = zeros(1,3)
        desp2 = zeros(1,3)
        coordn1prima = zeros(1,3)
        coordn2prima = zeros(1,3)
        indices= zeros(12)
        desp_locales = zeros(12)
        fuerza_locales = zeros(12)
        eps_elem = zeros(m,3)
        sigma_elem = zeros(m,3)

        for i in 1:m
            n1 = 1 + conect[i,1]
            n2 = 1 + conect[i,2]
            # n3 = 1 + conect[i,3]

            E = propiedades[1]
            ν = propiedades[2]

            L = longitud_elementos[i]
            R = √(A[i] / π)

            J = A[i]^2 / (2π)

            #PARTE AXIL
            coordn1 = coord[n1,1:3]
            coordn2 = coord[n2,1:3]
            L = sqrt((coordn2[1]-coordn1[1])^2 + (coordn2[2]-coordn1[2])^2 + (coordn2[3]-coordn1[3])^2)

            desp1 = desp_lin[n1,1:3]
            desp2 = desp_lin[n2,1:3]
            coordn1prima = coordn1 + desp1
            coordn2prima = coordn2 + desp2
            Lprima = sqrt((coordn2prima[1]-coordn1prima[1])^2 + (coordn2prima[2]-coordn1prima[2])^2 + (coordn2prima[3]-coordn1prima[3])^2)

            eps_axil = (Lprima-L)/L
            sigma_axil = eps_axil*E

            # GIRO DESPLAZAMIENTOS Y FUERZAS
            gdl1 = (6*n1)-5
            gdl2 = 6*n1
            gdl3 = (6*n2)-5
            gdl4 = 6*n2
            G = E/(2*(1 + ν))
            indices = [gdl2-5 gdl2-4 gdl2-3 gdl2-2 gdl2-1 gdl2 gdl4-5 gdl4-4 gdl4-3 gdl4-2 gdl4-1 gdl4]

            desp_aux = displacement[indices]
            fuerza_aux = force[indices]

            fuerza_locales, desp_locales = giro(coord[n1,:], coord[n2,:], desp_aux, fuerza_aux)

            # ELECCION DE GIRO EN Y
            if abs(desp_locales[5])>abs(desp_locales[11])
                θ_y = abs(desp_locales[11])
            else
                θ_y = abs(desp_locales[5])
            end

            #ELECCION DE GIRO EN Z
            if abs(desp_locales[6])>abs(desp_locales[12])
                θ_z = abs(desp_locales[12])
            else
                θ_z = abs(desp_locales[6])
            end

            #CALCULO DE TENSIONES EN EL ELEMENTO
            sigma_elem[i,1] = sigma_axil
            sigma_elem[i,2] = (sigma_elem[i,1]/abs(sigma_elem[i, 1])) * abs((θ_y*R*4*E)/longitud_elementos[i] + (θ_z*R*4*E)/longitud_elementos[i])
            sigma_elem[i,3] = fuerza_locales[10]*R / J

            #CALCULO DE DEFORMACIONES EN EL ELEMENTO
            eps_elem[i,1] = eps_axil
            eps_elem[i,2] = sigma_elem[i, 2] / E
            eps_elem[i,3] = sigma_elem[i, 3] / G

            # Compute strain and total energy of element
            # σ_e = [sigma_elem[i, 1] + sigma_elem[i, 2] 0 sigma_elem[i, 3]]
            # ε_e = [eps_elem[i, 1] + eps_elem[i, 2] 0 eps_elem[i, 3]]

        end

        #CREAR CONECTIn3IDAD DE NODOS
        # conect_nodos = Array{Int64}[]
        # #LO LLENAMOS DE TANTOS VECTORES COMO NODOS TIENE LA ESTRUCTURA
        # for _ in 1:n
        #     h = []
        #     push!(conect_nodos,h)
        # end
        conect_nodos = [[] for _ in 1:n]

        #METEMOS, EN EL VECTOR DE CADA NODO, TODOS LOS ELEMENTOS QUE LO TOCAN
        for i in 1:m
            n1 = 1 + conect[i,1]
            n2 = 1 + conect[i,2]
            push!(conect_nodos[n1],i)
            push!(conect_nodos[n2],i)
        end

        #ASIGNAMOS LAS DEFORMACIONES Y TENSIONES A LOS NODOS QUE PERTENEZCAN A CADA ELEMENTO
        #LOS DATOS DE TENSIONES Y DEFORMACIONES SE ACUMULAN EN DOS MATRICES
        #CADA VECTOR TIENE 9 COMPONENTES: 3 DE AXIAL, 3 DE FLECTOR, 3 DE CORTANTE  (MEDIA, MAX, MIN)

        eps_nod = zeros(n,9)
        sigma_nod = zeros(n,9)
        elementos = []
        for i in 1:n
            elementos = conect_nodos[i]

            #PARTE DE DEFORMACIONES
            eps_axial = eps_elem[elementos,1]
            media1 = mean(eps_axial)
            max1 = maximum(eps_axial)
            min1 = minimum(eps_axial)
            eps_nod[i,1:3] = [media1 max1 min1]

            eps_flec = eps_elem[elementos,2]
            media2 = mean(eps_flec)
            max2 = maximum(eps_flec)
            min2 = minimum(eps_flec)
            eps_nod[i,4:6] = [media2 max2 min2]

            eps_cort = eps_elem[elementos,3]
            media3 = mean(eps_cort)
            max3 = maximum(eps_cort)
            min3 = minimum(eps_cort)
            eps_nod[i,7:9] = [media3 max3 min3]

            #PARTE DE TENSIONES
            sigma_axial = sigma_elem[elementos,1]
            media4 = mean(sigma_axial)
            max4 = maximum(sigma_axial)
            min4 = minimum(sigma_axial)
            sigma_nod[i,1:3] = [media4 max4 min4]

            sigma_flec = sigma_elem[elementos,2]
            media5 = mean(sigma_flec)
            max5 = maximum(sigma_flec)
            min5 = minimum(sigma_flec)
            sigma_nod[i,4:6] = [media5 max5 min5]

            sigma_cort = sigma_elem[elementos,3]
            media6 = mean(sigma_cort)
            max6 = maximum(sigma_cort)
            min6 = minimum(sigma_cort)
            sigma_nod[i,7:9] = [media6 max6 min6]
        end

        return displacement, desp_lin, desp_ang, force, load, momento, sigma_elem, eps_elem, sigma_nod, eps_nod, energy
    end

    function solver(
        gdl::Int64, n::Int64, m::Int64, conect::Matrix{Int64}, coord::Matrix{Float64},
        propiedades::Array{Float64}, A::Vector{Float64}, index::Array{Any},
        desp::Array{Float64}, force::Array{Float64}, verbose::Bool=false,
        calculate_poisson::Bool=false, u::Float64=0.
    )

        elements_length = compute_elements_length(coord, conect)
        volume = compute_volume(A, elements_length)
        Kglobal = build_Kglobal(gdl, m, coord, conect, propiedades, elements_length, A)

        displacement_vector, desp_lin, desp_ang, fuerza, force, momento, σ_elem, ε_elem, _, _, energy = _solver(
            gdl, n, m, conect, coord, propiedades, index, desp, force, Kglobal, A, elements_length
        )
        if verbose
            println("Successful run. J = $(round(2 * energy, digits=4)) (Ψ = $(round(energy / volume, digits=2)) N/m^2)")
        end

        if calculate_poisson
            ν_xz, ν_yz = poisson_approximation(n, coord, desp_lin, u)
        else
            ν_xz = 0; ν_yz = 0;
        end

        return Dict(
            "displacement_vector" => displacement_vector,
            "force_vector" => fuerza,
            "disp_lin_elem" => desp_lin,
            "disp_ang_elem" => desp_ang,
            "force_elem" => force,
            "moment_elem" => momento,
            "σ_elem" => σ_elem,
            "ε_elem" => ε_elem,
            "energy" => energy,
            "compliance" => 2 * energy,
            "volumen" => volume,
            "lengths" => elements_length,
            "Ψ" => energy / volume,
            "Kglobal" => sparse(Kglobal),
            "E_z" => 0,
            "ν_xz" => ν_xz,
            "ν_yz" => ν_yz
        )
    end

    function solver_2d(
        dof::Int64, n::Int64, m::Int64, connect::Matrix{Int64}, coord::Matrix{Float64},
        properties::Array{Float64}, A::Vector{Float64}, dof_r::Vector{Int64},
        disp::Array{Float64}, force::Array{Float64}, filter::Bool=true, verbose::Bool=false
    )
    """Solve 2d lattice (truss structure).
    """
        elements_length = compute_elements_length(coord, connect, start=1)

        E = properties[1]
        # Declare unit vectors
        e1 = [1, 0, 0]
        e2 = [0, 1, 0]
        e3 = [0, 0, 1]
        # Allocate global sitffness matrix
        Kglobal = zeros(dof, dof)
        # Loop for stiffness matrix assembly
        for i = 1:m
            # Get the nodes for element i
            node1, node2 = connect[i, :]
            # Vector local x' axis
            t = coord[node2, :] - coord[node1, :]
            # Normalize it
            e1a = t / norm(t)
            # Unit vector local y' axis
            e2a = e3 × e1a
            Ke_local = E * A[i] / elements_length[i] * [
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
            Ke = transpose(Q) * Ke_local * Q
            dof_local = [2*node1 - 1, 2*node1, 2*node2 - 1, 2*node2]
            # Assembly
            Kglobal[dof_local, dof_local] += Ke
        end

        ##########################
        ### SOLVE (once assembled)
        ##########################
        displacement = copy(disp)

        dof_f = setdiff(1:dof, dof_r)

        # System reduction
        f_f = force[dof_f]
        # K_fr = Kglobal[dof_f, dof_r]
        # K_rf = transpose(K_fr)

        K_ff = Kglobal[dof_f, dof_f]
        # K_rr = Kglobal[dof_r, dof_r]
        # u_r = displacement[dof_r]

        # Condensación estática
        C = cholesky(sparse(K_ff))
        b = f_f# - K_fr*u_r
        u_f = C \ b
        # f_r = K_rf*u_f + K_rr * u_r

        displacement[dof_f] = u_f
        # force[dof_r] = f_r

        # Energy calculation
        compliance = force ⋅ displacement

        ### POST-PROCESS
        nodal_displacements = hcat(displacement[1:2:end], displacement[2:2:end], zeros(n))
        nodal_forces = hcat(force[1:2:end], force[2:2:end])

        element_forces = get_member_forces_2d(
            m, connect, coord, E, A, elements_length, displacement
        )

        if verbose
            println("Successful run. J = $(round(compliance, digits=4))")
        end

        str = filter ? "_filtered" : ""

        return Dict(
            "displacement_vector$str" => displacement,
            "nodal_displacements$str" => nodal_displacements,
            "nodal_forces$str" => nodal_forces,
            "element_forces$str" => element_forces,
            "compliance$str" => compliance,
            "Kfactorization$str" => C
        )
    end

    function solver_2d(
        dof::Int64, n::Int64, m::Int64, connect::Matrix{Int64}, coord::Matrix{Float64},
        E::Vector{Float64}, A::Float64, dof_r::Vector{Int64},
        disp::Array{Float64}, force::Array{Float64}, filter::Bool=true, verbose::Bool=false
    )
    """Solve 2d lattice (truss structure).
    """
        elements_length = compute_elements_length(coord, connect, start=1)

        # Declare unit vectors
        e1 = [1, 0, 0]
        e2 = [0, 1, 0]
        e3 = [0, 0, 1]
        # Allocate global sitffness matrix
        Kglobal = zeros(dof, dof)
        # Loop for stiffness matrix assembly
        for i = 1:m
            # Get the nodes for element i
            node1, node2 = connect[i, :]
            # Vector local x' axis
            t = coord[node2, :] - coord[node1, :]
            # Normalize it
            e1a = t / norm(t)
            # Unit vector local y' axis
            e2a = e3 × e1a
            Ke_local = E[i] * A / elements_length[i] * [
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
            Ke = transpose(Q) * Ke_local * Q
            dof_local = [2*node1 - 1, 2*node1, 2*node2 - 1, 2*node2]
            # Assembly
            Kglobal[dof_local, dof_local] += Ke
        end

        ##########################
        ### SOLVE (once assembled)
        ##########################
        displacement = copy(disp)

        dof_f = setdiff(1:dof, dof_r)

        # System reduction
        f_f = force[dof_f]
        # K_fr = Kglobal[dof_f, dof_r]
        # K_rf = transpose(K_fr)

        K_ff = Kglobal[dof_f, dof_f]
        # K_rr = Kglobal[dof_r, dof_r]
        # u_r = displacement[dof_r]

        # Condensación estática
        C = cholesky(sparse(K_ff))
        b = f_f# - K_fr*u_r
        u_f = C \ b
        # f_r = K_rf*u_f + K_rr * u_r

        displacement[dof_f] = u_f
        # force[dof_r] = f_r

        # Energy calculation
        compliance = force ⋅ displacement

        ### POST-PROCESS
        nodal_displacements = hcat(displacement[1:2:end], displacement[2:2:end], zeros(n))
        nodal_forces = hcat(force[1:2:end], force[2:2:end])

        if verbose
            println("Successful run. J = $(round(compliance, digits=4))")
        end

        str = filter ? "_filtered" : ""

        return Dict(
            "displacement_vector$str" => displacement,
            "nodal_displacements$str" => nodal_displacements,
            "nodal_forces$str" => nodal_forces,
            "element_forces$str" => zeros(m, 3),
            "compliance$str" => compliance,
            "Kfactorization$str" => C
        )
    end

    function get_member_forces_2d(
        m::Int64, connect::Matrix{Int64}, coord::Matrix{Float64},
        E::Float64, A::Vector{Float64}, lengths::Vector{Float64},
        u::Array{Float64}
    )
    # Allocate forces and stresses vectors
        forces = zeros(m)
        # Declare unit vectors
        e1 = [1, 0, 0]
        e2 = [0, 1, 0]
        e3 = [0, 0, 1]
        # Loop for stiffness matrix assembly
        for e = 1:m
            # Get the nodes for element e
            node1, node2 = connect[e, :]
            # Vector local x' axis
            t = coord[node2, :] - coord[node1, :]
            # Normalize it
            e1a = t / norm(t)
            # Unit vector local y' axis
            e2a = e3 × e1a
            # Rotation matrix
            T = [e1⋅e1a e2⋅e1a; e1⋅e2a e2⋅e2a]
            # Q matrix
            Q = [T zeros(2, 2); zeros(2, 2) T]
            # Get displacements affecting local nodes
            dof_local = [2*node1 - 1, 2*node1, 2*node2 - 1, 2*node2]
            ue = u[dof_local]
            # Rotate member displacements to local reference frame
            ue_local = Q * ue
            forces[e] = E * A[e] * (ue_local[3] - ue_local[1]) / lengths[e]
        end
        return forces
    end

    function solver_3d(
        dof::Int64, m::Int64, connect::Matrix{Int64}, properties::Vector{Float64}, A::Vector{Float64},
        lengths::Vector{Float64}, angles::Matrix{Float64}, dof_r::Vector{Int64}, disp::Array{Float64}, force::Array{Float64},
        filter::Bool=true, verbose::Bool=false
    )
    """Solve 3d lattice (truss structure).
    """
        E = properties[1]
        # Allocate global sitffness matrix
        Kglobal = spzeros(dof, dof)

        # Loop for stiffness matrix assembly
        print("Assembly..."); @time for i = 1:m
            # Get the nodes for element i
            node1, node2 = connect[i, :]
            # Compute Kloc in global axis
            Ke = truss_3d(E * A[i], lengths[i], angles[i, 1], angles[i, 2], angles[i, 3])
            dof_local = [3*node1 - 2, 3*node1 - 1, 3*node1, 3*node2 - 2, 3*node2 - 1, 3*node2]
            # Assembly
            Kglobal[dof_local, dof_local] += Ke
        end

        ##########################
        ### SOLVE (once assembled)
        ##########################
        displacement = copy(disp)

        dof_f = setdiff(1:dof, dof_r)

        # System reduction
        f_f = force[dof_f]
        K_fr = Kglobal[dof_f, dof_r]
        K_rf = transpose(K_fr)

        K_ff = Kglobal[dof_f,dof_f]
        K_rr = Kglobal[dof_r, dof_r]
        u_r = displacement[dof_r]

        # Condensación estática
        b = f_f #- K_fr * u_r
        C = cholesky(K_ff)
        u_f = C \ b

        f_r = K_rf*u_f + K_rr * u_r

        displacement[dof_f] = u_f
        force[dof_r] = f_r

        # Energy calculation
        compliance = force ⋅ displacement

        ### POST-PROCESS
        nodal_displacements = hcat(displacement[1:3:end], displacement[2:3:end], displacement[3:3:end])
        nodal_forces = hcat(force[1:3:end], force[2:3:end], force[3:3:end])
        # elem_forces = zeros(m, 3)
        elem_forces = get_member_forces_3d(m, connect, E, A, lengths, angles, displacement)

        if verbose
            println("Successful run. J = $(round(compliance, digits=4))")
        end

        str = filter ? "_filtered" : ""

        return Dict(
            "displacement_vector$str" => displacement,
            "nodal_displacements$str" => nodal_displacements,
            "nodal_forces$str" => nodal_forces,
            "element_forces$str" => elem_forces,
            "compliance$str" => compliance,
            "Kfactorization$str" => C
        )
    end

    function get_member_forces_3d(
        m::Int64, connect::Matrix{Int64}, E::Float64, A::Vector{Float64},
        lengths::Vector{Float64}, angles::Matrix{Float64}, u::Array{Float64}
    )
        # Allocate forces and stresses vectors
        forces = zeros(m)
        # Loop for stiffness matrix assembly
        for e = 1:m
            # Get the nodes for element e
            node1, node2 = connect[e, :]
            # Get displacements affecting local nodes
            dof_local = [3*node1 - 2, 3*node1 - 1, 3*node1, 3*node2 - 2, 3*node2 - 1, 3*node2]
            ue = u[dof_local]
            # Rotate difference of displacements between node2 and node1
            Δue_local = (ue[4:6] - ue[1:3]) ⋅ angles[m, :]
            #
            forces[e] = E * A[e] * Δue_local / lengths[e]
        end
        return forces
    end

    function get_member_stresses(
        connect::Matrix{Int64}, coord::Matrix{Float64},
        E::Float64, lengths::Vector{Float64}, nodal_displacements::Matrix{Float64}
    )
        new_coord = coord + nodal_displacements
        new_lengths = compute_elements_length(new_coord, connect, start=1)

        strains = new_lengths ./ lengths .- 1.
        stresses = E .* strains

        return stresses
    end


    function elongation_calculation(n, coord, deformed, limits::Vector; axis::Int=1)
        plano_inf_Z = []
        plano_sup_Z = []

        for i in 1:n
            if coord[i, axis] ≤ limits[1]
                push!(plano_inf_Z, i)
            elseif coord[i, axis] ≥ limits[2]
                push!(plano_sup_Z, i)
            end
        end

        coord_plano_inf_Z = sort(deformed[plano_inf_Z, axis])
        coord_plano_sup_Z = sort(deformed[plano_sup_Z, axis])

        limite1_inf_Z = trunc(Int, length(coord_plano_inf_Z) * 0.) + 1
        limite2_inf_Z = trunc(Int, length(coord_plano_inf_Z) * 1.)

        limite1_sup_Z = trunc(Int, length(coord_plano_sup_Z) * 0.) + 1
        limite2_sup_Z = trunc(Int, length(coord_plano_sup_Z) * 1.)

        mean_displacement_inf = mean(coord_plano_inf_Z[limite1_inf_Z:limite2_inf_Z])
        mean_displacement_sup = mean(coord_plano_sup_Z[limite1_sup_Z:limite2_sup_Z])

        return mean_displacement_sup - mean_displacement_inf
    end

    function poisson_approximation(n, coord, disp_lin, u)
        """Poisson ratios estimation
        """
        # Deformed configuration
        deformed = coord + disp_lin

        # Lengths
        xmin, ymin, zmin = minimum(coord, dims=1)
        xmax, ymax, zmax = maximum(coord, dims=1)

        dx = xmax - xmin
        dy = ymax - ymin
        dz = zmax - zmin

        # Elongations
        elongation_x = elongation_calculation(n, coord, deformed, [xmin, xmax], axis=1)
        elongation_y = elongation_calculation(n, coord, deformed, [ymin, ymax], axis=2)

        # Strains
        ε_xx_Z = (elongation_x - dx) / dx
        ε_yy_Z = (elongation_y - dy) / dy
        ε_zz_Z = u / dz

        # Poisson ratios
        ν_xz = -ε_xx_Z / ε_zz_Z
        ν_yz = -ε_yy_Z / ε_zz_Z

        return ν_xz, ν_yz
    end

    function truss_3d(EA::Float64, L::Float64, cx::Float64, cy::Float64, cz::Float64)
        subK = (EA/L) .* [
            cx^2 cx*cy cx*cz;
            cx*cy cy^2 cy*cz;
            cx*cz cy*cz cz^2
        ]
        return [subK -subK; -subK subK]
    end

end
