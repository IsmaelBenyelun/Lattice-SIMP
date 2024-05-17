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
using CPUTime

push!(LOAD_PATH, pwd())
using Generation
using Geometry, Filters
using LoadCases, FEM
using BsplinesUtils
using Objectives, Sensitivity
using PostProcess

DATA_FOLDER = ".\\output-files"

f = 1. # 1e3
E = 100. # 2e9

volfrac = .2

# Filtering and penalisation
FILTER = true
R = 20.
# p = 3.
p = get_interpolations()

Amin = 1e-10
Amax = 1.

# Initial point and tolerance
s0 = volfrac
TOL = 1e-2

mutable struct Lattice3d_
    data_folder::String
    mesh_parameters::Dict
    geom_parameters::Dict
    material_properties::Vector{Float64}
    s::Vector{Float64}
    shat::Vector{Float64}
    load_case::Tuple # restricted_dof, displacement_value, force_value, free_dof
    static_solution::Dict
    counter::Int64
end

function build_lattice(lattice::Lattice3d_, coord_file::String, connect_file::String, load_case_name::String, discretisation::String, filter::Bool=FILTER, load::Float64=f, load_nature::String="force")
    """Runs metamaterial being fixed all the inputs
    """
    # Generate metamaterial
    lattice.mesh_parameters = generate_lattice(
        coord_file, connect_file
    )
    # Define load case
    if load_case_name == "uniaxial"
        load_case_function = load_case_uniaxial_2d
    elseif load_case_name == "cantilever"
        load_case_function = load_case_cantilever_2d
    elseif load_case_name == "test"
        load_case_function = load_case_test_3d
    elseif load_case_name == "GE_bracket"
        load_case_function = load_case_GE_bracket
    end

    # print("Loading load case...")
    lattice.load_case = load_case_function(
        load,
        lattice.mesh_parameters["dof"],
        lattice.mesh_parameters["coord"],
        discretisation,
        load_nature
    )

    # Compute geometry parameters (fixed)
    # print("Computing elements length...")
    lattice.geom_parameters["lengths"] = compute_elements_length(
        lattice.mesh_parameters["coord"], lattice.mesh_parameters["connect"], start=1
    )

    # print("Computing angles...")
    lattice.geom_parameters["angles"] = compute_angles(
        lattice.mesh_parameters["coord"], lattice.mesh_parameters["connect"], lattice.geom_parameters["lengths"]
    )

    # Apply filter (if selected)
    # print("Computing centroids...")
    lattice.geom_parameters["centroids"] = centroids(
        lattice.mesh_parameters["coord"], lattice.mesh_parameters["connect"], filter, start=1
    )

    print("Computing weights..."); @time lattice.geom_parameters["weights"] = weights(
        lattice.geom_parameters["centroids"], R, filter
    )

    return nothing
end

function compute_obj(lattice::Lattice3d_, verbose::Bool=false, filter::Bool=FILTER)
    # Run/compute obj
    print("Iteration $(lattice.counter). ")
    # Filter variables
    lattice.shat = weight_filter(lattice.s, lattice.geom_parameters["weights"], filter)
    # Areas
    areas = compute_penalised_property(lattice.s, p, Amin, Amax)
    areas_hat = compute_penalised_property(lattice.shat, p, Amin, Amax)
    # Solve model
    # print("Solving system...")
    lattice.static_solution = solver_3d(
        lattice.mesh_parameters["dof"],
        lattice.mesh_parameters["m"],
        lattice.mesh_parameters["connect"],
        lattice.material_properties,
        areas_hat,
        lattice.geom_parameters["lengths"],
        lattice.geom_parameters["angles"],
        lattice.load_case[1], # dof_r
        lattice.load_case[2], # displacement array
        lattice.load_case[3], # force array
        false,
        verbose
    )

    # Retrive compliance
    J = lattice.static_solution["compliance"]

    # Calculate gradients
    dJ_dx = compliance_sensitivity_3d(
        lattice.mesh_parameters["m"],
        lattice.mesh_parameters["connect"],
        lattice.geom_parameters["lengths"],
        lattice.geom_parameters["angles"],
        lattice.material_properties,
        lattice.static_solution["displacement_vector"],
        lattice.geom_parameters["weights"],
        lattice.shat,
        p,
        Amin,
        Amax
    )

    # Post-process/save csv
    save_results_2d(
        J, 0, J, 1., 0, 0, 0, 0, lattice.data_folder
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
    return J, dJ_dx, lattice.static_solution["compliance"]
end


function objective(x::Vector, grad::Vector, lattice::Lattice3d_)
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

function constraint(x::Vector, grad::Vector, lattice::Lattice3d_, filter::Bool=FILTER)
    """Function to be wrapped as constraint
    """
    areas = compute_penalised_property(x, 1., Amin, Amax) # p = 1 trick
    volume = compute_volume(areas, lattice.geom_parameters["lengths"])

    g = volume - sum(Amax .* lattice.geom_parameters["lengths"]) * volfrac
    dg_dx = volume_sensitivity(lattice.shat, 1., lattice.geom_parameters["lengths"], lattice.geom_parameters["weights"], Amin, Amax, filter_type="standard") # p = 1 trick

    # Return suitable values for NLopt algorithm
    grad[:] = dg_dx
    println("Constr: $(round(g, digits=4))")
    return g
end

function optimize_subproblem(x0::Vector, lattice::Lattice3d_, tol::Float64=TOL)
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
    return minf, minx
end

# Create data_folder
# load_case_name = "uniaxial"
# load_case_name = "cantilever"
load_case_name = "GE_bracket"

timestamp = Dates.format(now(UTC), "YYYY-mm-dd_HH_MM_SS")
data_folder = DATA_FOLDER * "\\$timestamp-LC-$load_case_name"
mkdir(data_folder)

discretisation = "coarsest"
# discretisation = "coarse"
# discretisation = "finer"

# Define files
# coord_file = "data/nodes_test.csv"; connect_file = "data/connect_test.csv"
coord_file = "data/nodes_$discretisation.csv"; connect_file = "data/connect_$discretisation.csv"

k = 0

# Instantiate Lattice
lattice = Lattice3d_(
    data_folder,
    Dict([]),
    Dict([]),
    [E],
    [],
    [],
    (),
    Dict([]),
    k
)
build_lattice(lattice, coord_file, connect_file, load_case_name, discretisation)

# Create input set and compute
x0 = s0 * ones(lattice.mesh_parameters["m"]);
# x0 = rand(Normal(s0, 1e-1), lattice.mesh_parameters["m"])

### Run optimization
# lattice.s = x0; @time out = compute_obj(lattice);
out = optimize_subproblem(x0, lattice);
println("Optimization finished.")
