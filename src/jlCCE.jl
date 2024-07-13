module jlCCE

# physical constants
const hbar = 1.054571817e-34 # J s -- CODATA 2022 (rounded value)

# struct to hold parameters for the run
mutable struct SpinSystem
    # file with coordinates (anything that the file read can handle)
    coord_file::String
    # name of the atom defining the spin center ("Cu", "V", ...)
    spin_center::String
    # index within unit cell to select spin center, if name is not unique
    spin_center_index::Int
    # g factor at the center (x,y,z)
    g_factor::Vector{Float64}(missing,3)
    # magnetic axes of spin center (column vectors for x, y, z direction)
    magnetic_axes::Matrix{Float64}(missing,3,3)
    # nuclei defining spin bath ("H", etc.)
    nuc_spin_bath::String
    # corresponding nuclear g factor
    gn_spin_bath::Float64
    # magnetic field
    B0::Vector{Float64}(missing,3)
    # minimum interaction radius (usually 0.)
    r_min::Float64
    # maximum interaction radius
    r_max::Float64
    # min and maxiumum of time interval
    t_min::Float64
    t_max::Float64
    # number of time steps in interval
    n_time_step::Int
end

"""
    cce(system::SpinSystem)

    explanatory text goes here
"""
function cce(system::SpinSystem)
    # load coordinates
    # identify spin center
    # get list of spin bath nuclei
    
    # precompute A values for electron nucleus pairs (function call)

    # initialize Intensity to 1. for all times
    # double loop m>n over nuclear pairs
        # call a function to:
        # compute b, c, omega for each pair, no need to store on array
        # compute v for all t values of simulation (for this m,n)
        # probabaly best: compute Intensity = Intensity .* exp(v)

    # return intensity and time
end

end