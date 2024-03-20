module BsplinesUtils

    using NURBS, StaticArrays, Interpolations
    export get_interpolations

    function bspline_curve(degree, ctrlpoints)
        knots = vcat(zeros(degree + 1), 0.5, ones(degree + 1))
        return BsplineCurve(Bspline(degree, knots), ctrlpoints)
    end

    function bspline_curve1()
        ctrlpts = [
            SA[0.0, 0.0, 0.0],
            SA[0.15, 0.005, 0.0],
            SA[0.55, 0.01, 0.0],
            SA[0.45, 0.45, 0.0],
            SA[.5, 0.5, 0.0]
        ]
        return bspline_curve(3, ctrlpts)
    end

    function bspline_curve2()
        ctrlpts = [
            SA[0.0, 0.0, 0.0],
            SA[0.15, 0.0, 0.0],
            SA[0.3, 0.005, 0.0],
            SA[0.525, 0.01, 0.0],
            SA[0.45, 0.45, 0.0],
            SA[.5, 0.5, 0.0]
        ]
        return bspline_curve(4, ctrlpts)
    end

    function bspline_curve3()
        ctrlpts = [
            SA[0.0, 0.0, 0.0],
            SA[0.15, 0.0, 0.0],
            SA[0.3, 0.005, 0.0],
            SA[0.55, 0.01, 0.0],
            SA[0.45, 0.45, 0.0],
            SA[.5, 0.5, 0.0]
        ]
        return bspline_curve(4, ctrlpts)
    end

    function bspline_curve5()
        ctrlpts = [
            SA[0.0, 0.0, 0.0],
            SA[0.35, 0.0, 0.0],
            SA[0.45, 0.45, 0.0],
            SA[.5, 0.5, 0.0]
        ]
        return bspline_curve(2, ctrlpts)
    end

    function _get_interpolations(NBspline::BsplineCurve)
        evalpoints = collect(0:0.0001:1.)
        C = curveDerivativesPoints(NBspline, evalpoints, 1)
        fun = transpose(reinterpret(Float64, C[1]')) # parametrized spline
        dfun = transpose(reinterpret(Float64, C[2]')) # hodograph of spline

        x = fun[:, 1]; y = fun[:, 2]
        dx_du = dfun[:, 1]; dy_du = dfun[:, 2]

        f_interpolated = interpolate((x,), y, Gridded(Linear()))
        df_interpolated = interpolate((x,), dy_du ./ dx_du, Gridded(Linear()))

        return (f_interpolated, df_interpolated)
    end

    function get_interpolations(option::String="normal")
        if option == "normal"
            NBspline = bspline_curve2()
        elseif option == "smoother"
            NBspline = bspline_curve5()
        end

        return _get_interpolations(NBspline)
    end

end #module
