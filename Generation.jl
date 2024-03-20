module Generation

    using DelimitedFiles
    export generate_lattice

    #################################################
    ### LATTICE FUNCTIONS
    #################################################
    function get_nodes_and_elements(nodes_h, nodes_w)
        # Number of nodes (n) and elements (m)
        elems_w = nodes_w - 1
        elems_h = nodes_h - 1
        n = nodes_w * nodes_h
        m = elems_w * elems_h * 4 + elems_h + elems_w
        return n, m
    end

    function create_coord(nodes_h, nodes_w, l=1.; scale_w=1., scale_h=1.)
        # Create coord
        coord = zeros(nodes_h * nodes_w, 3)
        k = 1
        for i in 1:nodes_h
            for j in 1:nodes_w
                # push!(coord, [k, (j - 1)*l, (i - 1)*l, 0])
                coord[k, :] = [(j - 1)* scale_w * l, (i - 1)* scale_h * l, 0.]
                k += 1
            end
        end
        return coord
    end

    function create_connect(nodes_h, nodes_w)
        # Create connect
        elems_h = nodes_h - 1
        elems_w = nodes_w - 1

        _, m = get_nodes_and_elements(nodes_h, nodes_w)

        connect = zeros(Int64, m, 2)
        e = 1
        for i in 1:elems_h
            for j in 1:elems_w
                a = (i - 1) * nodes_w + j
                b = (i - 1) * nodes_w + j + 1
                c = i * nodes_w + j
                d = i * nodes_w + j + 1
                connect[e, :] = [a, b]
                connect[e + 1, :] = [a, c]
                connect[e + 2, :] = [b, c]
                connect[e + 3, :] = [a, d]

                e += 4
                if j == elems_w
                    connect[e, :] = [b, d]
                    e += 1
                end
                if i == elems_h
                    connect[e, :] = [c, d]
                    e += 1
                end
            end
        end

        return connect
    end

    function generate_lattice(dimensions, l, data_path=".\\test"; scale_w=1., scale_h=1.)
        # Unpack
        nodes_h, nodes_w = dimensions
        # Calculate
        n, m = get_nodes_and_elements(nodes_h, nodes_w)
        coord = create_coord(nodes_h, nodes_w, l, scale_w=scale_w, scale_h=scale_h)
        connect = create_connect(nodes_h, nodes_w)
        # Return
        return Dict(
            "coord" => coord, #convert(Matrix{Float64}, coord),
            "connect" => connect, #convert(Matrix{Int64}, connect),
            "n" => n,
            "m" => m,
            "dof" => 2 * n,
            "cube_vol" => prod(maximum(coord[:, 1:2], dims=1) - minimum(coord[:, 1:2], dims=1))
        )
    end

    function generate_lattice(coord_file::String, connect_file::String)
        coord = readdlm(coord_file, ',', Float64)
        connect = readdlm(connect_file, ',', Int64)
        n = size(coord, 1); m = size(connect, 1); dof = 3*n
        return Dict(
            "coord" => coord,
            "connect" => connect,
            "n" => n,
            "m" => m,
            "dof" => dof
        )
    end

end # module
