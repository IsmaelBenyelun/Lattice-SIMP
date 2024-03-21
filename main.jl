"""Script for Lattice Framework (suitable for compute, optimization, etc.)

    - Geometrical parameters configuration
    - Load case configuration
    - Materials set
    - Domain (dimensions) set
"""

using NLopt
using Dates
using Random
using TickTock
using Distributions
using LinearAlgebra
using DelimitedFiles
using BenchmarkTools

push!(LOAD_PATH, pwd())
using Generation, LoadCases
using Geometry, Filters, FEM
using Objectives, Sensitivity
using PostProcess
using BsplinesUtils

DATA_FOLDER = ".\\output-files"

if !isdir(DATA_FOLDER)
    mkdir(DATA_FOLDER)
end

f = 1e2 # 1e3
E = 7e7 # 2e9
# f = 1.
# E = 100.
L = .5

volfrac = .4
# Vmax = 4.376487985345041 #6.0
# Vmax = L^2 / 2

# Filtering and penalisation
FILTER = true
R = 1.
p = 3.
# p = "poly"
# p = get_interpolations("smoother")

Amin = 1e-10
# Amax = 1.
Amax = 5.1e-3

# Initial point and tolerance
s0 = volfrac
# s0 = 1e-3
TOL = 1e-4

# Mean - variance hyperparameter
# μ_star = 0.028; σ_star = 0.00038
μ_star = 1.; σ_star = 1.

# μ_star = 0.5337583952657886; σ_star = 0.002186224069578384 # p = 1
# μ_star = 0.6208861130873831; σ_star = 1 # p = 2
# μ_star = 0.6380135644974949; σ_star = 1 # p = 3
# μ_star = 0.6326254525851839; σ_star = 1 # p = 4
# μ_star = .6086383797555627; σ_star = 1. # spline (smoother)
# μ_star = 0.6589174145979452; σ_star = 1. # spline

α = 1.
α_E = α / μ_star
α_S = (1 - α) / σ_star

# Matérn GP hyperparameters
σ = 0.1E
ℓ = .5
ν = 1.5

mutable struct Lattice
    data_folder::String
    dimensions::Vector{Any}
    mesh_parameters::Dict
    geom_parameters::Dict
    material_properties::Vector{Float64}
    s::Vector{Float64}
    shat::Vector{Float64}
    load_case::Tuple # restricted_dof, displacement_value, force_value, free_dof
    static_solution::Dict
    counter::Int64
    covariance_matrix::Matrix{Float64}
end

function build_lattice(lattice::Lattice, load_case_name::String, filter::Bool=FILTER, load::Float64=f, load_nature::String="force")
    """Runs metamaterial being fixed all the inputs
    """
    # Generate metamaterial
    lattice.mesh_parameters = generate_lattice(
        lattice.dimensions,
        L,
        lattice.data_folder,
        scale_h=1. #.75
    )
    # Define load case
    if load_case_name == "uniaxial"
        load_case_function = load_case_uniaxial_2d
    elseif load_case_name == "cantilever"
        load_case_function = load_case_cantilever_2d
    end

    lattice.load_case = load_case_function(
        load,
        lattice.mesh_parameters["dof"],
        lattice.dimensions,
        lattice.mesh_parameters["coord"],
        load_nature
    )

    # Compute geometry parameters (fixed)
    lattice.geom_parameters["lengths"] = compute_elements_length(
        lattice.mesh_parameters["coord"], lattice.mesh_parameters["connect"], start=1
    )

    # Apply filter (if selected)
    lattice.geom_parameters["centroids"] = centroids(
        lattice.mesh_parameters["coord"], lattice.mesh_parameters["connect"], filter, start=1
    )
    lattice.geom_parameters["weights"] = weights(
        lattice.geom_parameters["centroids"], R, filter
    )
    # lattice.geom_parameters["∂filter"] = filter_derivative(
    #     lattice.geom_parameters["weights"], lattice.geom_parameters["lengths"], filter
    # )
    return nothing
end

function compute_obj(lattice::Lattice, verbose::Bool=false, filter::Bool=FILTER)
    # Run/compute obj
    print("Iteration $(lattice.counter). ")
    # Filter variables
    lattice.shat = weight_filter(lattice.s, lattice.geom_parameters["weights"], lattice.geom_parameters["lengths"], filter)
    # Areas
    areas = compute_penalised_property(lattice.s, p, Amin, Amax)
    areas_hat = compute_penalised_property(lattice.shat, p, Amin, Amax)
    # Solve model
    dict_normal = solver_2d(
        lattice.mesh_parameters["dof"],
        lattice.mesh_parameters["n"],
        lattice.mesh_parameters["m"],
        lattice.mesh_parameters["connect"],
        lattice.mesh_parameters["coord"],
        lattice.material_properties,
        areas,
        lattice.load_case[1], # dof_r
        lattice.load_case[2], # displacement array
        lattice.load_case[3], # force array
        false,
        verbose
    )
    dict_filtered = solver_2d(
        lattice.mesh_parameters["dof"],
        lattice.mesh_parameters["n"],
        lattice.mesh_parameters["m"],
        lattice.mesh_parameters["connect"],
        lattice.mesh_parameters["coord"],
        lattice.material_properties,
        areas_hat,
        lattice.load_case[1], # dof_r
        lattice.load_case[2], # displacement array
        lattice.load_case[3], # force array
        true,
        verbose
    )
    merge!(lattice.static_solution, dict_normal)
    merge!(lattice.static_solution, dict_filtered)

    J = lattice.static_solution["compliance"]
    # Compute sensitivities (required for computing J and dJ_dx)
    dJ_dx = compliance_sensitivity_2d(
        lattice.mesh_parameters["m"],
        lattice.mesh_parameters["coord"],
        lattice.mesh_parameters["connect"],
        lattice.geom_parameters["lengths"],
        areas_hat,
        lattice.material_properties,
        lattice.static_solution["displacement_vector_filtered"],
        lattice.geom_parameters["weights"],
        lattice.shat,
        p,
        Amin,
        Amax
    )

    # println(norm(dJ_dx))
    # Post-process/save csv
    save_results_2d(
        J, 0, J, 1., 0., 0., 0., 0, lattice.data_folder
    )
    # Post-process/save vtk
    create_vtk_2d(
        lattice.mesh_parameters["coord"],
        lattice.mesh_parameters["connect"],
        lattice.static_solution["nodal_displacements"],
        lattice.static_solution["nodal_forces"],
        lattice.static_solution["element_forces"],
        areas,
        areas_hat,
        lattice.s,
        lattice.shat,
        lattice.counter,
        lattice.data_folder
    )
    return α_E * J, α_E * dJ_dx, lattice.static_solution["compliance"]
end

function objective(x::Vector, grad::Vector, lattice::Lattice)
    """Function to be wrapped as objective
    """
    # Increment counter
    lattice.counter += 1
    # Assign design space
    lattice.s = x
    J, dJ_dx, compliance = compute_obj(lattice)
    # Return suitable values for NLopt algorithm
    grad[:] = dJ_dx
    # Print value
    print("Obj:  $(round(J, digits=4)). (compliance = $(round(compliance, digits=4))). ")
    return J
end

function constraint(x::Vector, grad::Vector, lattice::Lattice, filter::Bool=FILTER)
    """Function to be wrapped as constraint
    """
    areas = compute_penalised_property(x, 1., Amin, Amax) # p = 1 trick
    # areas = compute_penalised_property(x, p, Amin, Amax) # p = 1 trick
    volume = compute_volume(areas, lattice.geom_parameters["lengths"])
    g = volume - sum(Amax .* lattice.geom_parameters["lengths"]) * volfrac
    # g = volume - Vmax

    dg_dx = volume_sensitivity(lattice.shat, 1., lattice.geom_parameters["lengths"], lattice.geom_parameters["weights"], Amin, Amax) # p = 1 trick
    # dg_dx = volume_sensitivity(lattice.shat, p, lattice.geom_parameters["lengths"], lattice.geom_parameters["weights"], Amin, Amax) # p = 1 trick
    # Return suitable values for NLopt algorithm
    grad[:] = dg_dx
    println("Constr: $(round(g, digits=4))")
    return g
end

function optimize_subproblem(x0::Vector, lattice::Lattice, tol::Float64=TOL)
    # Dimensions
    m = lattice.mesh_parameters["m"]

    # Create opt instance
    opt = Opt(:LD_MMA, m)
    opt.lower_bounds = zeros(m)
    opt.upper_bounds = ones(m)
    opt.xtol_rel = tol

    # Assign objective and constraints
    opt.min_objective = (x, g) -> objective(x, g, lattice)
    inequality_constraint!(opt, (x, g) -> constraint(x, g, lattice), 1e-8)

    # Run optimization
    minf, minx, _ = optimize(opt, x0)

    ### Save info about last iteration (x, ∂J/∂x)
    lattice.counter += 1000
    J = lattice.static_solution["compliance"]
    save_results_2d(
        J, 0, J, 1., .1E, 0., 0., 0, lattice.data_folder
    )
    return minf, minx
end

# Create data_folder
# load_case_name = "uniaxial"
load_case_name = "cantilever"

timestamp = Dates.format(now(UTC), "YYYY-mm-dd_HH_MM_SS")
data_folder = DATA_FOLDER * "\\$timestamp-LC-$load_case_name"
mkdir(data_folder)

# Define parameters
# nodes_h = 31
# nodes_w = 21

nodes_h = 21 #41
nodes_w = 41 #81

# nodes_h = 2
# nodes_w = 2

# nodes_h = 3
# nodes_w = 5

dimensions = [nodes_h, nodes_w]
mesh_params = Dict([])
geom_params = Dict([])
properties = [E]

k = 0

# Instantiate Lattice
lattice = Lattice(
    data_folder,
    dimensions,
    mesh_params,
    geom_params,
    properties,
    [],
    [],
    (),
    Dict([]),
    k,
    zeros(0, 0)
)
build_lattice(lattice, load_case_name)

# Create input set and compute
x0 = s0 * ones(lattice.mesh_parameters["m"])
# x0 = rand(Normal(s0, 1e-1), lattice.mesh_parameters["m"])

### Run optimization
# lattice.s = x0; out = compute_obj(lattice);
out = optimize_subproblem(x0, lattice);
println("Optimization finished.")
