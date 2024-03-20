module Filters

    using LinearAlgebra, SparseArrays
    export weight_filter, filter_derivative

    function weight_filter(x, w, lengths, filter::Bool)
        if filter
            numerator = w * (x ./ lengths)
            denominator = w * (1 ./ lengths)
            return numerator ./ denominator
        else
            return x
        end
    end

    function weight_filter(x, w, filter::Bool)
        if filter
            numerator = w * x
            denominator = vec(sum(w, dims=1))
            return numerator ./ denominator
        else
            return x
        end
    end

    function filter_derivative(w, lengths, filter::Bool)
        """The derivatives are constant (i.e. not dependant on x). Only computed once
        """
        if !filter
            m = length(lengths)
            return I(m)
        end

        numerator = w ./ transpose(lengths)
        denominator = vec(sum(w ./ lengths, dims=1))
        return numerator ./ denominator
    end

    function filter_derivative(w)
        """The derivatives are constant (i.e. not dependant on x). Only computed once
        """
        return Matrix(w ./ sum(w, dims=2))
    end

end #module
