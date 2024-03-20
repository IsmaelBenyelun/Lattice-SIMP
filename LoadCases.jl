module LoadCases

    using Statistics
    push!(LOAD_PATH, pwd())
    using FEM

    export load_case_X, load_case_Y, load_case_Z, load_case_Z_param, uniaxial_Z, uniaxial_Z_local, uniaxial_Z_chair, load_case_uniaxial_2d, load_case_cantilever_2d, load_case_test_3d, load_case_GE_bracket

    function restricted_gdls(n, coord)
        """(Common) Restricted nodes applicable to 3 uniaxial load cases
        Symmetry assumption
        """
        ##########################################################
        ################  CONDICIONES DE CONTORNO  ###############
        ##########################################################

        index_aux = []
        # index_reacciones_Z = []        #ESTE VECTOR CONTIENE LOS INDICES DE LOS NODOS DEL PLANO INFERIOR Z (POR SI ES NECESARIO EN ALGUN MOMENTO)
        # index_reacciones_X = []        #ESTE VECTOR CONTIENE LOS INDICES DE LOS NODOS DEL PLANO INFERIOR X
        # index_reacciones_Y = []        #ESTE VECTOR CONTIENE LOS INDICES DE LOS NODOS DEL PLANO INFERIOR Y

        # condicion_planoz = 0.3
        condicion_planoz = 0.1
        index_aux_Z = []    #ESTE VECTOR CONTIENE LOS GDL QUE ESTAN RESTRINGIDOS DEL PLANO Z INFERIOR

        # condicion_planoy = 0.4
        condicion_planoy = 0.1
        index_aux_Y = []    #ESTE VECTOR CONTIENE LOS GDL QUE ESTAN RESTRINGIDOS DEL PLANO Y INFERIOR

        # condicion_planox = 0.5
        condicion_planox = 0.1
        index_aux_X = []    #ESTE VECTOR CONTIENE LOS GDL QUE ESTAN RESTRINGIDOS DEL PLANO X INFERIOR


        for i in 1:n
            gdlx = 6*i - 5
            gdly = 6*i - 4
            gdlz = 6*i - 3

            if coord[i, 3] ≤ condicion_planoz
                push!(index_aux_Z, gdlz)
                # push!(index_reacciones_Z, gdlz)
            end

            if coord[i, 1] ≤ condicion_planox
                push!(index_aux_X, gdlx)
                # push!(index_reacciones_X, gdlx)
            end

            if coord[i, 2] ≤ condicion_planoy
                push!(index_aux_Y, gdly)
                # push!(index_reacciones_Y, gdly)
            end
        end

        append!(index_aux, index_aux_Z)
        append!(index_aux, index_aux_X)
        append!(index_aux, index_aux_Y)

        return index_aux
    end

    function restricted_gdls_disk(n, coord)
        """(Common) Restricted nodes applicable to 3 uniaxial load cases
        Symmetry assumption
        """
        ##########################################################
        ################  CONDICIONES DE CONTORNO  ###############
        ##########################################################

        disk_radius = 2
        x0 = 6.5
        y0 = 6.5

        index_aux = []
        # index_reacciones_Z = []        #ESTE VECTOR CONTIENE LOS INDICES DE LOS NODOS DEL PLANO INFERIOR Z (POR SI ES NECESARIO EN ALGUN MOMENTO)
        # index_reacciones_X = []        #ESTE VECTOR CONTIENE LOS INDICES DE LOS NODOS DEL PLANO INFERIOR X
        # index_reacciones_Y = []        #ESTE VECTOR CONTIENE LOS INDICES DE LOS NODOS DEL PLANO INFERIOR Y

        # condicion_planoz = 0.3
        condicion_planoz = 0.1
        index_aux_Z = []    #ESTE VECTOR CONTIENE LOS GDL QUE ESTAN RESTRINGIDOS DEL PLANO Z INFERIOR

        # condicion_planoy = 0.4
        condicion_planoy = 0.1
        index_aux_Y = []    #ESTE VECTOR CONTIENE LOS GDL QUE ESTAN RESTRINGIDOS DEL PLANO Y INFERIOR

        # condicion_planox = 0.5
        condicion_planox = 0.1
        index_aux_X = []    #ESTE VECTOR CONTIENE LOS GDL QUE ESTAN RESTRINGIDOS DEL PLANO X INFERIOR


        for i in 1:n
            gdlx = 6*i - 5
            gdly = 6*i - 4
            gdlz = 6*i - 3

            x, y, z = coord[i, :]

            if x ≤ condicion_planox
                push!(index_aux_X, gdlx)
            end

            if y ≤ condicion_planoy
                push!(index_aux_Y, gdly)
            end

            if (z ≤ condicion_planoz) && ((x - x0)^2 + (y - y0)^2 ≤ disk_radius^2)
                push!(index_aux_Z, gdlx)
                push!(index_aux_Z, gdly)
                push!(index_aux_Z, gdlz)
            end
        end

        append!(index_aux, index_aux_Z)
        append!(index_aux, index_aux_X)
        append!(index_aux, index_aux_Y)

        return index_aux
    end

    function load_case_X(n, m, gdl, coord, conect, propiedades, Kglobal, longitud_elementos, modulo_torsion)
        """Load case in X
        """
        ##########################################################
        ###########  CARGAS Y CONDICIONES DE CONTORNO  ###########
        ##########################################################
        index_aux = restricted_gdls(n, coord)
        #INICIALIZAMOS VECTORES DE FUERZAS Y DESPLAZAMIENTOS
        desp_X = zeros(gdl,1)
        fuerza_X = zeros(gdl,1)

        #CASO DE CARGA
        condicion_desp_x = maximum(coord[:,1])
        condicion_desp_z1_X = -1
        condicion_desp_z2_X = 100
        condicion_desp_y1_X = -1
        condicion_desp_y2_X = 100
        desplazamiento_x = -0.05
        impuesto_X = []            #SIRVE PARA VER EN DONDE SE HAN COLOCADO LOS DESPLAZAMIENTOS EN CASO DE QUERER COMPROBARLO
        for i in 1:n
            gdlx = (6*i)-5
            if coord[i,1]>=condicion_desp_x && condicion_desp_z1_X<=coord[i,3]<=condicion_desp_z2_X && condicion_desp_y1_X<=coord[i,2]<=condicion_desp_y2_X
                desp_X[gdlx] = desplazamiento_x
                push!(impuesto_X, gdlx)
            end
        end

        append!(impuesto_X, index_aux)

        index_X = sort(impuesto_X)     #ORDENA EL VECTOR DE INDICES DE MENOR A MAYOR

        ##########################################################################
        ##########  RESOLUCION Y CALCULO DEFORMACIONES Y TENSIONES  ##############
        ##########################################################################

        desp_lin_X, desp_ang_X, fuerza_X, force_X, momento_X, sigma_elem_X, eps_elem_X, sigma_nod_X, eps_nod_X, energy = solver(gdl, n, m, conect, coord, propiedades, index_X, desp_X, Kglobal, longitud_elementos, modulo_torsion)

        ##########################################################
        ############  CALCULO YOUNG Y POISSON TOTALES  ###########
        ##########################################################

        dx = maximum(coord[:,1]) - minimum(coord[:,1])
        dy = maximum(coord[:,2]) - minimum(coord[:,2])
        dz = maximum(coord[:,3]) - minimum(coord[:,3])

        coord_prima_X = coord + desp_lin_X

        #DEFORMACION EN Z
        planoz_inferior_X = []
        planoz_superior_X = []
        for i in 1:n
            if coord[i,3]<=minimum(coord[:,3])
                push!(planoz_inferior_X, i)
            elseif coord[i,3]>=maximum(coord[:,3])
                push!(planoz_superior_X, i)
            end
        end

        coord_planoz_inferior_X = sort(coord_prima_X[planoz_inferior_X,3])
        coord_planoz_superior_X = sort(coord_prima_X[planoz_superior_X,3])

        limite1zinferior_X = trunc(Int, length(coord_planoz_inferior_X) * 0.25)
        limite2zinferior_X = trunc(Int, length(coord_planoz_inferior_X) * 0.75)
        limite1zsuperior_X = trunc(Int, length(coord_planoz_superior_X) * 0.25)
        limite2zsuperior_X = trunc(Int, length(coord_planoz_superior_X) * 0.75)

        dz_prima_X = mean(coord_planoz_superior_X[limite1zsuperior_X:limite2zsuperior_X]) - mean(coord_planoz_inferior_X[limite1zinferior_X:limite2zinferior_X])

        eps_z_X = (dz_prima_X - dz)/dz

        #DEFORMACION EN Y
        planoy_inferior_X = []
        planoy_superior_X = []
        for i in 1:n
            if coord[i,2]<=minimum(coord[:,2])
                push!(planoy_inferior_X, i)
            elseif coord[i,2]>=maximum(coord[:,2])
                push!(planoy_superior_X, i)
            end
        end

        coord_planoy_inferior_X = sort(coord_prima_X[planoy_inferior_X,2])
        coord_planoy_superior_X = sort(coord_prima_X[planoy_superior_X,2])

        limite1yinferior_X = trunc(Int, length(coord_planoy_inferior_X) * 0.25)
        limite2yinferior_X = trunc(Int, length(coord_planoy_inferior_X) * 0.75)
        limite1ysuperior_X = trunc(Int, length(coord_planoy_superior_X) * 0.25)
        limite2ysuperior_X = trunc(Int, length(coord_planoy_superior_X) * 0.75)

        dy_prima_X = mean(coord_planoy_superior_X[limite1ysuperior_X:limite2ysuperior_X]) - mean(coord_planoy_inferior_X[limite1yinferior_X:limite2yinferior_X])

        eps_y_X = (dy_prima_X - dy)/dy

        #DEFORMACION EN X
        eps_x_X = desplazamiento_x/dx

        poisson_zx = -eps_z_X/eps_x_X
        poisson_yx = -eps_y_X/eps_x_X
        # poisson_total_X = (poisson_zx+poisson_yx)/2

        fuerzas_planox_inferior = []
        for i in 1:n
            if coord[i,1]<=minimum(coord[:,1])
                gdlx = (6*i)-5
                push!(fuerzas_planox_inferior, gdlx)
            end
        end

        reacciones_X = fuerza_X[fuerzas_planox_inferior]
        fuerza_total_X = sum(reacciones_X)

        tension_media_X = fuerza_total_X/(dz*dy)
        E_x = tension_media_X/abs(eps_x_X)

        println(eps_x_X)
        println(reacciones_X)
        println(fuerza_total_X)

        return E_x, poisson_zx, poisson_yx, tension_media_X, eps_x_X, energy, desp_lin_X, desp_ang_X, force_X, sigma_elem_X, eps_elem_X
    end

    function load_case_Y(n, m, gdl, coord, conect, propiedades, Kglobal, longitud_elementos, modulo_torsion)
        """Load case in Y
        """
        ##########################################################
        ###########  CARGAS Y CONDICIONES DE CONTORNO  ###########
        ##########################################################
        index_aux = restricted_gdls(n, coord)
        #INICIALIZAMOS VECTORES DE FUERZAS Y DESPLAZAMIENTOS
        desp_Y = zeros(gdl,1)
        fuerza_Y = zeros(gdl,1)

        #CASO DE CARGA
        condicion_desp_y = maximum(coord[:,2])
        condicion_desp_x1_Y = -1
        condicion_desp_x2_Y = 100
        condicion_desp_z1_Y = -1
        condicion_desp_z2_Y = 100
        desplazamiento_y = -0.05
        impuesto_Y = []            #SIRVE PARA VER EN DONDE SE HAN COLOCADO LOS DESPLAZAMIENTOS EN CASO DE QUERER COMPROBARLO
        for i in 1:n
            gdly = (6*i)-4
            if coord[i,2]>=condicion_desp_y && condicion_desp_x1_Y<=coord[i,1]<=condicion_desp_x2_Y && condicion_desp_z1_Y<=coord[i,3]<=condicion_desp_z2_Y
                desp_Y[gdly] = desplazamiento_y
                push!(impuesto_Y, gdly)
            end
        end

        append!(impuesto_Y, index_aux)

        index_Y = sort(impuesto_Y)      #ORDENA EL VECTOR DE INDICES DE MENOR A MAYOR

        ##########################################################################
        ##########  RESOLUCION Y CALCULO DEFORMACIONES Y TENSIONES  ##############
        ##########################################################################

        desp_lin_Y, desp_ang_Y, fuerza_Y, force_Y, momento_Y, sigma_elem_Y, eps_elem_Y, sigma_nod_Y, eps_nod_Y, energy = solver(gdl, n, m, conect, coord, propiedades, index_Y, desp_Y, Kglobal, longitud_elementos, modulo_torsion)

        ##########################################################
        ############  CALCULO YOUNG Y POISSON TOTALES  ###########
        ##########################################################

        dx = maximum(coord[:,1]) - minimum(coord[:,1])
        dy = maximum(coord[:,2]) - minimum(coord[:,2])
        dz = maximum(coord[:,3]) - minimum(coord[:,3])

        coord_prima_Y = coord + desp_lin_Y

        #DEFORMACION EN X
        planox_inferior_Y = []
        planox_superior_Y = []
        for i in 1:n
            if coord[i,1]<=minimum(coord[:,1])
                push!(planox_inferior_Y, i)
            elseif coord[i,1]>=maximum(coord[:,1])
                push!(planox_superior_Y, i)
            end
        end

        coord_planox_inferior_Y = sort(coord_prima_Y[planox_inferior_Y,1])
        coord_planox_superior_Y = sort(coord_prima_Y[planox_superior_Y,1])

        limite1xinferior_Y = trunc(Int, length(coord_planox_inferior_Y) * 0.25)
        limite2xinferior_Y = trunc(Int, length(coord_planox_inferior_Y) * 0.75)
        limite1xsuperior_Y = trunc(Int, length(coord_planox_superior_Y) * 0.25)
        limite2xsuperior_Y = trunc(Int, length(coord_planox_superior_Y) * 0.75)

        dx_prima_Y = mean(coord_planox_superior_Y[limite1xsuperior_Y:limite2xsuperior_Y]) - mean(coord_planox_inferior_Y[limite1xinferior_Y:limite2xinferior_Y])

        eps_x_Y = (dx_prima_Y - dx)/dx

        #DEFORMACION EN Z
        planoz_inferior_Y = []
        planoz_superior_Y = []
        for i in 1:n
            if coord[i,3]<=minimum(coord[:,3])
                push!(planoz_inferior_Y, i)
            elseif coord[i,3]>=maximum(coord[:,3])
                push!(planoz_superior_Y, i)
            end
        end

        coord_planoz_inferior_Y = sort(coord_prima_Y[planoz_inferior_Y,3])
        coord_planoz_superior_Y = sort(coord_prima_Y[planoz_superior_Y,3])

        limite1zinferior_Y = trunc(Int, length(coord_planoz_inferior_Y) * 0.25)
        limite2zinferior_Y = trunc(Int, length(coord_planoz_inferior_Y) * 0.75)
        limite1zsuperior_Y = trunc(Int, length(coord_planoz_superior_Y) * 0.25)
        limite2zsuperior_Y = trunc(Int, length(coord_planoz_superior_Y) * 0.75)

        dz_prima_Y = mean(coord_planoz_superior_Y[limite1zsuperior_Y:limite2zsuperior_Y]) - mean(coord_planoz_inferior_Y[limite1zinferior_Y:limite2zinferior_Y])

        eps_z_Y = (dz_prima_Y - dz)/dz

        #DEFORMACION EN Y
        eps_y_Y = desplazamiento_y/dy

        poisson_xy = -eps_x_Y/eps_y_Y
        poisson_zy = -eps_z_Y/eps_y_Y
        # poisson_total_Y = (poisson_xy+poisson_zy)/2

        fuerzas_planoy_inferior = []
        for i in 1:n
            if coord[i,2]<=minimum(coord[:,2])
                gdly = (6*i)-4
                push!(fuerzas_planoy_inferior, gdly)
            end
        end

        reacciones_Y = fuerza_Y[fuerzas_planoy_inferior]
        fuerza_total_Y = sum(reacciones_Y)

        tension_media_Y = fuerza_total_Y/(dx*dz)
        E_y = tension_media_Y/abs(eps_y_Y)

        return E_y, poisson_xy, poisson_zy, tension_media_Y, eps_y_Y, energy, desp_lin_Y, desp_ang_Y, force_Y, sigma_elem_Y, eps_elem_Y
    end

    function load_case_Z(n, m, gdl, coord, conect, propiedades, Kglobal, longitud_elementos, modulo_torsion)
        """Load case in Z
        """
        return load_case_Z_param(-0.05, n, m, gdl, coord, conect, propiedades, Kglobal, longitud_elementos, modulo_torsion)

    end

    function uniaxial_Z(load, n, gdl, coord, load_nature)
        # index_aux = restricted_gdls(n, coord)
        index_aux = restricted_gdls_disk(n, coord)
        # INICIALIZAMOS VECTOR DE DESPLAZAMIENTOS
        value_u = zeros(gdl)
        value_f = zeros(gdl)

        # CASO DE CARGA
        condicion_desp_z = maximum(coord[:,3])
        condicion_desp_x1_Z = -1
        condicion_desp_x2_Z = 100
        condicion_desp_y1_Z = -1
        condicion_desp_y2_Z = 100

        impuesto = [] # SIRVE PARA VER EN DONDE SE HAN COLOCADO LOS DESPLAZAMIENTOS EN CASO DE QUERER COMPROBARLO

        for i in 1:n
            gdlz = 6*i - 3
            # For z ≡ top face of the cube, and -1 ≤ x ≤ 100, and -1 ≤ y ≤ 100
            if (coord[i, 3] ≥ condicion_desp_z) && (condicion_desp_x1_Z ≤ coord[i, 1] ≤ condicion_desp_x2_Z) && (condicion_desp_y1_Z ≤ coord[i, 2] ≤ condicion_desp_y2_Z)
                if load_nature == "force"
                    value_f[gdlz] = load
                elseif load_nature == "displacement"
                    value_u[gdlz] = load
                    push!(impuesto, gdlz)
                end
            end
        end

        append!(impuesto, index_aux)

        restricted_dof = sort(impuesto)      #ORDENA EL VECTOR DE INDICES DE MENOR A MAYOR
        free_dof = setdiff(1:gdl, restricted_dof)
        return restricted_dof, value_u, value_f, load, free_dof
    end

    function uniaxial_Z_local(load, n, gdl, coord, load_nature)
        index_aux = restricted_gdls(n, coord)
        # INICIALIZAMOS VECTOR DE DESPLAZAMIENTOS
        value_u = zeros(gdl)
        value_f = zeros(gdl)

        # CASO DE CARGA
        condicion_desp_z = maximum(coord[:, 3])
        radius = 2

        impuesto = [] # SIRVE PARA VER EN DONDE SE HAN COLOCADO LOS DESPLAZAMIENTOS EN CASO DE QUERER COMPROBARLO

        for i in 1:n
            gdlz = 6*i - 3
            # For z ≡ top face of the cube, and x^2 + y^2 ≤ radius^2 (create a circle)
            if (coord[i, 3] ≥ condicion_desp_z) && (coord[i, 1]^2 + coord[i, 2]^2 ≤ radius^2)
                if load_nature == "force"
                    value_f[gdlz] = load
                elseif load_nature == "displacement"
                    value_u[gdlz] = load
                    push!(impuesto, gdlz)
                end
            end
        end

        append!(impuesto, index_aux)

        restricted_dof = sort(impuesto)      #ORDENA EL VECTOR DE INDICES DE MENOR A MAYOR
        free_dof = setdiff(1:gdl, restricted_dof)
        return restricted_dof, value_u, value_f, load, free_dof
    end

    function uniaxial_Z_chair(load, n, gdl, coord, load_nature)
        # Initialize displacement and force vectors
        value_u = zeros(gdl)
        value_f = zeros(gdl)
        impuesto = []

        # Get restricted displacement dofs
        index_aux = restricted_gdls_disk(n, coord)

        # Load geometry description
        condicion_desp_z = maximum(coord[:, 3])
        radius = 2

        for i in 1:n
            gdlz = 6*i - 3
            x, y, z = coord[i, :]
            # For z ≡ top face of the cube, and x^2 + y^2 ≤ radius^2 (create a circle)
            if (z ≥ condicion_desp_z) && ((x / radius)^2 + (y / radius)^2 ≤ 1)
                if load_nature == "force"
                    value_f[gdlz] = load
                elseif load_nature == "displacement"
                    value_u[gdlz] = load
                    push!(impuesto, gdlz)
                end

            end
        end

        append!(impuesto, index_aux)

        # Sort restricted index
        restricted_dof = sort(impuesto)
        free_dof = setdiff(1:gdl, restricted_dof)
        return restricted_dof, value_u, value_f, load, free_dof
    end

    function load_case_Z_param(u, n, m, gdl, coord, conect, propiedades, Kglobal, longitud_elementos, modulo_torsion)
        """Load case in Z
        """
        ####ª#####################################################
        ###########  CARGAS Y CONDICIONES DE CONTORNO  ###########
        ##########################################################
        index_Z, desp_Z = uniaxial_Z(u, 0, n, gdl, coord)

        ##########################################################################
        ############  RESOLUCION Y CALCULO DEFORMACIONES Y TENSIONES  ############
        ##########################################################################

        desp, desp_lin_Z, desp_ang_Z, fuerza_Z, force_Z, momento_Z, sigma_elem_Z, eps_elem_Z, sigma_nod_Z, eps_nod_Z, energy = _solver(gdl, n, m, conect, coord, propiedades, index_Z, desp_Z, Kglobal, longitud_elementos, modulo_torsion)

        ##########################################################
        ############  CALCULO YOUNG Y POISSON TOTALES  ###########
        ##########################################################

        ### Poisson ratios estimation
        xmin, ymin, zmin = minimum(coord, dims=1)
        xmax, ymax, zmax = maximum(coord, dims=1)

        dx = xmax - xmin
        dy = ymax - ymin
        dz = zmax - zmin

        coord_prima_Z = coord + desp_lin_Z

        # Displacement/strain in x
        planox_inferior_Z = []
        planox_superior_Z = []

        for i in 1:n
            if coord[i,1] ≤ xmin
                push!(planox_inferior_Z, i)
            elseif coord[i,1] ≥ xmax
                push!(planox_superior_Z, i)
            end
        end

        coord_planox_inferior_Z = sort(coord_prima_Z[planox_inferior_Z, 1])
        coord_planox_superior_Z = sort(coord_prima_Z[planox_superior_Z, 1])

        limite1xinferior_Z = trunc(Int, length(coord_planox_inferior_Z) * 0.)
        limite2xinferior_Z = trunc(Int, length(coord_planox_inferior_Z) * 0.3)

        limite1xsuperior_Z = trunc(Int, length(coord_planox_superior_Z) * 0.)
        limite2xsuperior_Z = trunc(Int, length(coord_planox_superior_Z) * 0.3)

        dx_prima_Z = mean(coord_planox_superior_Z[limite1xsuperior_Z:limite2xsuperior_Z]) - mean(coord_planox_inferior_Z[limite1xinferior_Z:limite2xinferior_Z])

        eps_x_Z = (dx_prima_Z - dx)/dx

        # Displacement/strain in y
        planoy_inferior_Z = []
        planoy_superior_Z = []

        for i in 1:n
            if coord[i,2] ≤ ymin
                push!(planoy_inferior_Z, i)
            elseif coord[i,2] ≥ ymax
                push!(planoy_superior_Z, i)
            end
        end

        coord_planoy_inferior_Z = sort(coord_prima_Z[planoy_inferior_Z, 2])
        coord_planoy_superior_Z = sort(coord_prima_Z[planoy_superior_Z, 2])

        limite1yinferior_Z = trunc(Int, length(coord_planoy_inferior_Z) * 0.25)
        limite2yinferior_Z = trunc(Int, length(coord_planoy_inferior_Z) * 0.75)

        limite1ysuperior_Z = trunc(Int, length(coord_planoy_superior_Z) * 0.25)
        limite2ysuperior_Z = trunc(Int, length(coord_planoy_superior_Z) * 0.75)

        dy_prima_Z = mean(coord_planoy_superior_Z[limite1ysuperior_Z:limite2ysuperior_Z]) - mean(coord_planoy_inferior_Z[limite1yinferior_Z:limite2yinferior_Z])

        eps_y_Z = (dy_prima_Z - dy)/dy

        # Displacement/strain in z
        eps_z_Z = u/dz

        poisson_xz = -eps_x_Z/eps_z_Z
        poisson_yz = -eps_y_Z/eps_z_Z
        # poisson_total_Z = (poisson_xz+poisson_yz)/2

        ### Young's modulus estimation
        # Reaction forces calculation
        fuerzas_planoz_inferior = []
        for i in 1:n
            if coord[i,3] ≤ zmin
                gdlz = (6*i)-3
                push!(fuerzas_planoz_inferior, gdlz)
            end
        end

        reacciones_Z = fuerza_Z[fuerzas_planoz_inferior]
        fuerza_total_Z = sum(reacciones_Z)

        tension_media_Z = fuerza_total_Z / (dx*dy)
        E_z = tension_media_Z / abs(eps_z_Z)

        return E_z, poisson_xz, poisson_yz, tension_media_Z, eps_z_Z, energy, desp, desp_lin_Z, desp_ang_Z, force_Z, sigma_elem_Z, eps_elem_Z
    end

    function restricted_dofs_2d(coord)
        restricted_nodes = findall(==(0), coord[:, 1])

        h_dofs = 2 * restricted_nodes .- 1
        v_dofs = 2 * restricted_nodes

        restricted_dofs = vcat(h_dofs, v_dofs)
        return sort(restricted_dofs)
    end

    function get_imposed_load_dofs_uniaxial_2d(coord, nodes_h, nodes_w)
        # ind1 = findall(==(nodes_w - 1), coord[:, 1])
        x_max = maximum(coord[:, 1])
        ind1 = findall(==(x_max), coord[:, 1])
        half_h = nodes_h ÷ 2
        # ind2 = [half_h] if nodes_h % 2 == 1 else [half_h - 1, half_h]
        ind2 = (nodes_h % 2 == 1) ? [half_h + 1] : [half_h, half_h + 1]

        selected_nodes = ind1[ind2]

        return 2 * selected_nodes .- 1
    end

    function get_imposed_load_dofs_cantilever_2d(coord, nodes_h, nodes_w)
        # ind1 = findall(==(nodes_w - 1), coord[:, 1])
        x_max = maximum(coord[:, 1])
        ind1 = findall(==(x_max), coord[:, 1])
        half_h = nodes_h ÷ 2
        # ind2 = [half_h] if nodes_h % 2 == 1 else [half_h - 1, half_h]
        ind2 = (nodes_h % 2 == 1) ? [half_h + 1] : [half_h, half_h + 1]

        selected_nodes = ind1[ind2]

        return 2 .* selected_nodes
    end

    function load_case_uniaxial_2d(load, dof, dimensions, coord, load_nature)
        # Unpack
        nodes_h, nodes_w = dimensions
        # Allocate u, f vectors
        value_u = zeros(dof)
        value_f = zeros(dof)

        # Restricted dofs (displacement = 0). The remaining are free
        restricted_dofs = restricted_dofs_2d(coord)
        free_dofs = setdiff(1:dof, restricted_dofs)
        # Imposed dof (either force or displacement, depending on load_nature)
        imposed_dofs = get_imposed_load_dofs_uniaxial_2d(coord, nodes_h, nodes_w)

        if load_nature == "force"
            value_f[imposed_dofs] .= load
        elseif load_nature == "displacement"
            value_u[imposed_dofs] .= load
            append!(restricted_dofs, imposed_dofs)
        end

        return sort(restricted_dofs), value_u, value_f, load, free_dofs
    end

    function load_case_cantilever_2d(load, dof, dimensions, coord, load_nature)
        # Unpack
        nodes_h, nodes_w = dimensions
        # Allocate u, f vectors
        value_u = zeros(dof)
        value_f = zeros(dof)

        # Restricted dofs (displacement = 0). The remaining are free
        restricted_dofs = restricted_dofs_2d(coord)
        free_dofs = setdiff(1:dof, restricted_dofs)
        # Imposed dof (either force or displacement, depending on load_nature)
        imposed_dofs = get_imposed_load_dofs_cantilever_2d(coord, nodes_h, nodes_w)

        if load_nature == "force"
            value_f[imposed_dofs] .= load
        elseif load_nature == "displacement"
            value_u[imposed_dofs] .= load
            append!(restricted_dofs, imposed_dofs)
        end

        return sort(restricted_dofs), value_u, value_f, load, free_dofs
    end

    function load_case_test_3d(load, dof, coord, load_nature)
        value_u = zeros(dof)
        value_f = zeros(dof)

        restricted_dofs = 1:12
        free_dofs = setdiff(1:dof, restricted_dofs)
        imposed_dofs = [15]

        if load_nature == "force"
            value_f[imposed_dofs] .= load
        elseif load_nature == "displacement"
            value_u[imposed_dofs] .= load
            append!(restricted_dofs, imposed_dofs)
        end

        return Vector(restricted_dofs), value_u, value_f, load, free_dofs
    end

    function restricted_dofs_GE_bracket(coord)
        x_leq_20 = findall(<=(20.), coord[:, 1])
        x_geq_160 = findall(>=(160.), coord[:, 1])
        y_leq_20 = findall(<=(20.), coord[:, 2])
        y_geq_100 = findall(>=(100.), coord[:, 2])

        corner1 = intersect(x_geq_160, y_leq_20)
        corner2 = intersect(x_geq_160, y_geq_100)
        corner3 = intersect(x_leq_20, y_leq_20)
        corner4 = intersect(x_leq_20, y_geq_100)

        restricted_nodes = vcat(corner1, corner2, corner3, corner4)
        x_dofs = 3 .* restricted_nodes .- 2
        y_dofs = 3 .* restricted_nodes .- 1
        z_dofs = 3 .* restricted_nodes

        return sort(vcat(x_dofs, y_dofs, z_dofs))
    end

    function get_imposed_load_dofs_GE_bracket(coord, discretisation)
        z_eq = 70.
        z_geq = 55.
        z_leq = 65.

        if (discretisation == "coarse") | (discretisation == "coarsest")
            y_eq = 15.
            y_geq = 20.
            y_leq = 30.

            corner_nodes = []
        elseif discretisation == "finer"
            y_eq = 10.
            y_geq = 15.
            y_leq = 25.

            y_c = 12.5
            z_c = 67.5

            y_eq_2 = findall(==(y_c), coord[:, 2])
            z_eq_2 = findall(==(z_c), coord[:, 3])
            corner_nodes = intersect(y_eq_2, z_eq_2)
        end


        y_eq_ = findall(==(y_eq), coord[:, 2])
        z_geq_ = findall(>=(z_geq), coord[:, 3])
        z_leq_ = findall(<=(z_leq), coord[:, 3])

        z_eq_ = findall(==(z_eq), coord[:, 3])
        y_geq_ = findall(>=(y_geq), coord[:, 2])
        y_leq_ = findall(<=(y_leq), coord[:, 2])

        vertical_nodes = intersect(y_eq_, z_geq_, z_leq_)
        horizontal_nodes = intersect(z_eq_, y_geq_, y_leq_)

        imposed_nodes = sort(vcat(vertical_nodes, horizontal_nodes, corner_nodes))

        y_dofs = 3 .* imposed_nodes .- 1
        z_dofs = 3 .* imposed_nodes
        return y_dofs, z_dofs
    end

    function load_case_GE_bracket(load, dof, coord, discretisation, load_nature)
        # Allocate u, f vectors
        value_u = zeros(dof)
        value_f = zeros(dof)

        # Restricted dofs (displacement = 0). The remaining are free
        restricted_dofs = restricted_dofs_GE_bracket(coord)
        free_dofs = setdiff(1:dof, restricted_dofs)
        # Imposed dof (either force or displacement, depending on load_nature)
        imposed_dofs_y, imposed_dofs_z = get_imposed_load_dofs_GE_bracket(coord, discretisation)

        if load_nature == "force"
            value_f[imposed_dofs_y] .= -load / sqrt(2)
            value_f[imposed_dofs_z] .= load / sqrt(2)
        elseif load_nature == "displacement"
            value_u[imposed_dofs] .= load / sqrt(2)
            append!(restricted_dofs, imposed_dofs)
        end

        return sort(restricted_dofs), value_u, value_f, load, free_dofs
    end

end
