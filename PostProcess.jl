module PostProcess

    using WriteVTK
    using DelimitedFiles

    export save_results, save_results_2d, create_vtk, create_vtk_2d

    ### Paraview
    function write_paraview(filename, points, cells, label, desp_lin, desp_ang, force, sigma_elem, eps_elem, area)#, W_elem)
        """Function to save in ParaView the results
        """
        # Create vtk file
        vtkfile = vtk_grid(filename, points, cells)
        # Append all the field information
        vtkfile["displacements", VTKPointData()] = transpose(desp_lin)
        vtkfile["giros", VTKPointData()] = transpose(desp_ang)
        vtkfile["label", VTKPointData()] = transpose(label)
        vtkfile["forces", VTKPointData()] = transpose(force)
        vtkfile["stress (axial)", VTKCellData()] = sigma_elem[:,1]
        vtkfile["stress (flector)", VTKCellData()] = sigma_elem[:,2]
        vtkfile["stress (torsion)", VTKCellData()] = sigma_elem[:,3]
        vtkfile["strain (axial)", VTKCellData()] = eps_elem[:,1]
        vtkfile["strain (flector)", VTKCellData()] = eps_elem[:,2]
        vtkfile["strain (torsion)", VTKCellData()] = eps_elem[:,3]
        vtkfile["cross-sectional area", VTKCellData()] = area
        # Save it
        return vtk_save(vtkfile)
    end

    function create_vtk(coord, connect, disp_lin, disp_ang, force, σ_elem, ε_elem, area, counter, data_folder)
        ####################################################
        #################  PRINT PARAVIEW  #################
        ####################################################
        celltype = VTKCellTypes.VTK_LINE
        cells = MeshCell[]

        n = size(coord)[1]
        m = size(connect)[1]

        label = 1:n

        for i in 1:m
            c = MeshCell(celltype, [1 + connect[i, 1], 1 + connect[i, 2]])
            push!(cells, c)
        end

        filename = data_folder * "\\LC-it-$counter"

        out = write_paraview(filename, transpose(coord), cells, label, disp_lin, disp_ang, force, σ_elem, ε_elem, area)
    end

    function write_paraview_2d(filename, points, cells, label, displacement, force, force_elem, area, densities=nothing)
        """Function to save in ParaView the results
        """
        # Create vtk file
        vtkfile = vtk_grid(filename, points, cells)
        # Append all the field information
        vtkfile["label", VTKPointData()] = transpose(label)
        vtkfile["nodal displacements", VTKPointData()] = transpose(displacement)
        vtkfile["nodal forces", VTKPointData()] = transpose(force)
        vtkfile["element forces", VTKCellData()] = force_elem
        vtkfile["cross-sectional area", VTKCellData()] = area
        vtkfile["densities", VTKCellData()] = typeof(densities) == Nothing ? zero(area) : densities
        # Save it
        return vtk_save(vtkfile)
    end

    function create_vtk_2d(coord, connect, displacement, force, force_elem, area, counter, data_folder)
        ####################################################
        #################  PRINT PARAVIEW  #################
        ####################################################
        celltype = VTKCellTypes.VTK_LINE
        cells = MeshCell[]

        n = size(coord)[1]
        m = size(connect)[1]

        label = 1:n

        for i in 1:m
            c = MeshCell(celltype, [connect[i, 1], connect[i, 2]])
            push!(cells, c)
        end

        filename = data_folder * "\\LC-it-$counter"

        # out = write_paraview_2d(filename, transpose(coord), cells, label, displacement, force, force_elem, area, densities)
        points = transpose(coord)

        # Create vtk file
        vtkfile = vtk_grid(filename, points, cells)
        # Append all the field information
        vtkfile["label", VTKPointData()] = transpose(label)
        vtkfile["nodal displacements", VTKPointData()] = transpose(displacement)
        vtkfile["nodal forces", VTKPointData()] = transpose(force)
        vtkfile["member forces", VTKCellData()] = force_elem
        vtkfile["member stresses", VTKCellData()] = force_elem ./ area
        vtkfile["cross-sectional area", VTKCellData()] = area
        # Save it
        return vtk_save(vtkfile)
    end

    function create_vtk_2d(coord, connect, displacement, force, force_elem, area, filtered_area, densities, filtered_densities, counter, data_folder)
        ####################################################
        #################  PRINT PARAVIEW  #################
        ####################################################
        celltype = VTKCellTypes.VTK_LINE
        cells = MeshCell[]

        n = size(coord)[1]
        m = size(connect)[1]

        label = 1:n

        for i in 1:m
            c = MeshCell(celltype, [connect[i, 1], connect[i, 2]])
            push!(cells, c)
        end

        filename = data_folder * "\\LC-it-$counter"

        # out = write_paraview_2d(filename, transpose(coord), cells, label, displacement, force, force_elem, area, densities)
        points = transpose(coord)

        # Create vtk file
        vtkfile = vtk_grid(filename, points, cells)
        # Append all the field information
        vtkfile["label", VTKPointData()] = transpose(label)
        vtkfile["nodal displacements", VTKPointData()] = transpose(displacement)
        vtkfile["nodal forces", VTKPointData()] = transpose(force)
        vtkfile["member forces", VTKCellData()] = force_elem
        vtkfile["member_stresses", VTKCellData()] = abs.(force_elem ./ (area))

        vtkfile["cross-sectional area", VTKCellData()] = area
        vtkfile["cross-sectional area (filtered)", VTKCellData()] = filtered_area
        vtkfile["densities", VTKCellData()] = densities
        vtkfile["densities (filtered)", VTKCellData()] = filtered_densities

        # Save it
        return vtk_save(vtkfile)
    end

    function create_vtk_2d_v8(coord, connect, displacement, force, force_elem, area, filtered_area, densities, filtered_densities, counter, data_folder)
        ####################################################
        #################  PRINT PARAVIEW  #################
        ####################################################
        celltype = VTKCellTypes.VTK_LINE
        cells = MeshCell[]

        n = size(coord)[1]
        m = size(connect)[1]

        label = 1:n

        for i in 1:m
            c = MeshCell(celltype, [connect[i, 1], connect[i, 2]])
            push!(cells, c)
        end

        filename = data_folder * "\\LC-it-$counter"

        # out = write_paraview_2d(filename, transpose(coord), cells, label, displacement, force, force_elem, area, densities)
        points = transpose(coord)

        # Create vtk file
        vtkfile = vtk_grid(filename, points, cells)
        # Append all the field information
        vtkfile["label", VTKPointData()] = transpose(label)
        vtkfile["nodal displacements", VTKPointData()] = transpose(displacement)
        vtkfile["nodal forces", VTKPointData()] = transpose(force)
        # vtkfile["member forces", VTKCellData()] = force_elem
        # vtkfile["member_stresses", VTKCellData()] = abs.(force_elem ./ (area))
        # vtkfile["member stresses", VTKCellData()] = abs.(force_elem ./ (filtered_area))

        vtkfile["member_forces", VTKCellData()] = force_elem .* area
        vtkfile["member forces", VTKCellData()] = force_elem .* filtered_area
        vtkfile["member stresses", VTKCellData()] = force_elem
        vtkfile["member stresses (abs)", VTKCellData()] = abs.(force_elem)

        vtkfile["cross-sectional area", VTKCellData()] = area
        vtkfile["cross-sectional area (filtered)", VTKCellData()] = filtered_area
        vtkfile["densities", VTKCellData()] = densities
        vtkfile["densities (filtered)", VTKCellData()] = filtered_densities

        # Save it
        return vtk_save(vtkfile)
    end

    ### Write CSV
    function save_results(
        geom_param, E_z, ν_xz, ν_yz, energy, Ψ, μJ, σJ, totalJ, κ, σf, ℓ, ν, n_MC, data_folder
    )
        ##########################################################
        #################  SALIDA DE RESULTADOS  #################
        ##########################################################

        results_filename = data_folder * "\\results.csv"
        header = "H,DJOINT,HEIGHT_STAR,DSTAR,meanJ,stdJ,totalJ,alpha,sigma_f,length_scale,v,n_MC,E_z,nu_xz,nu_yz,W,Psi\n"

        if !isfile(results_filename)
            open(results_filename, "w") do io
                write(io, header)
            end
        end

        result = transpose([
            geom_param; μJ; σJ; totalJ; κ; σf; ℓ; ν; n_MC; E_z; ν_xz; ν_yz; energy; Ψ
        ])
        # result = [transpose(geom_param) E_z ν_xz ν_yz energy Ψ]

        open(results_filename, "a") do io
            writedlm(io, result, ",")
        end

        return nothing
    end

    function save_results_2d(μJ, σJ, totalJ, κ, σf, ℓ, ν, n_MC, data_folder)
        results_filename = data_folder * "\\results.csv"
        header = "meanJ,stdJ,totalJ,K,sigma_f,length_scale,v,n_MC\n"

        if !isfile(results_filename)
            open(results_filename, "w") do io
                write(io, header)
            end
        end

        result = transpose([μJ, σJ, totalJ, κ, σf, ℓ, ν, n_MC])

        open(results_filename, "a") do io
            writedlm(io, result, ",")
        end

        return nothing
    end

end
