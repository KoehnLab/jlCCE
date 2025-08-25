module ParameterSpaceScreening

using LinearAlgebra

export get_distance_vectors,pss

# physical constants
# reduced Planck constant 
const hbar = 1.054571817e-27 # erg s -- J s -- CODATA 2022 (rounded value)
# Bohr magneton 
const mu_b = 9.274010066e-21 # erg G^-1 J T^-1 -- CODATA 2022 (rounded value)
# nuclear magneton
const mu_n = 5.050783739e-24 # erg G^-1  J T^-1 -- CODATA 2022 (rounded value)
# conversion AA to cm (cgs unit system)
const aacm = 1e-8

"""
    get_distance_vectors(R,R_12,alpha)

Determine the distance vectors of one electron spin and two nuclear spins using 
    Jacobi coordinates

input: R - distance between electron spin and center of mass of nuclear spins
        R_12 - distance between the nuclear spins
        alpha - angle between the two distance vectors

the function returns distance vectors between the electron spin and the two nuclear spins 
    as well as a distance vector between the nuclear spins 
"""
function get_distance_vectors(R::Float64,R_12::Float64,alpha::Float64)
    #println("")
    #println("Determine distance vectors using Jacobi coordinates")
    #println("---------------------------------------------------\n")
    vec_r = [R, 0., 0.]
    vec_r12 = [R_12 * cos(alpha), R_12 * sin(alpha), 0.]
    vec_r1 = vec_r + 1/2 * vec_r12
    vec_r2 = vec_r - 1/2 * vec_r12

    #println("Input distances (R, R_12): ",R," ",R_12)
    #println("Determined distance vectors: ")
    #println("Distance vector e-n1: ", vec_r1)
    #println("Distance vector e-n2: ", vec_r2)
    #println("Distance vector n1-n2: ", vec_r12)

    return vec_r1,vec_r2,vec_r12
end

""" 
    pss(vec_r1::Vector{Float64},vec_r2,vec_r12,g_e,g_n, time, B0,alpha)

Determine the modulation depth λ, the nuclear zero-quantum frequency ω and the Hahn echo intensity
    using Jacobe coordinates.

input: vec_r1 - distance coordinate of nuc. 1
        vec_r2 - distance coordinate of nuc. 2
        vec_r12 - distance coordinates between nuc. 1,2
        g_e - g factor of the free electron
        g_n - nuclear g factor
        time - time for the simulation of the Hahn echo intensity
        B0 - magnetic field

the function returns the modulation depth, the nuclear zero-quantum frequency and the intensity 
    for a given pair of Jacobi coordinates
"""
function pss(vec_r1::Vector{Float64},vec_r2::Vector{Float64},vec_r12::Vector{Float64},g_e::Float64,g_n::Float64, time::Float64, B0::Vector{Float64})
    #println("")
    #println("Parameter Space Screening")
    #println("-------------------------\n")
    
    #println("Running on ",Threads.nthreads()," threads\n")

    # determine the gryomagnetic ratio
    gamma_electron = (g_e * mu_b)/ hbar
    gamma_n = (g_n * mu_n) / hbar

    # rescale distance coordinates from AA to cm (cgs unit system)
    r1 = vec_r1 * aacm
    r2 = vec_r2 * aacm
    r12 = vec_r12 * aacm

    # precompute values of the hyperfine coupling constant A for electron nucleus pairs
    # nuc. 1
    r1_cross_B0 = cross(r1,B0)
    r1_dot_B0 = dot(r1,B0)
    theta_1 = atan(norm(r1_cross_B0),r1_dot_B0)
    r1_norm = max(norm(r1),1e-8)
    A_1 = -gamma_n * gamma_electron * hbar * (1 - 3 * cos(theta_1)^2) / r1_norm^3
    
    # nuc. 2
    r2_cross_B0 = cross(r2,B0)
    r2_dot_B0 = dot(r2,B0)
    theta_2 = atan(norm(r2_cross_B0),r2_dot_B0)
    r2_norm = max(norm(r2),1e-8)
    A_2 = -gamma_n * gamma_electron * hbar * (1 - 3 * cos(theta_2)^2) / r2_norm^3
    #println("Hyperfine coupling constants A of nuc. 1,2 in Hz: ", A_1," ",A_2)
    
    # precompute value of the coupling constant b 
    r12_cross_B0 = cross(r12,B0);
    r12_dot_B0 = dot(r12,B0);
    theta_12 = atan(norm(r12_cross_B0),r12_dot_B0);
    b_12 = (- 1/4) * gamma_n^2 * hbar * ((1 - 3* cos(theta_12)^2)/(norm(r12))^3); 
    #println("Coupling constant b in Hz: ",b_12)

    # compute modulation depth
    c_12 = (A_1-A_2)/(4*b_12)
    modulation_depth = c_12^2/((1 + c_12^2)^2)
    #println("Modulation depth: ", modulation_depth)

    # compute nuclear zero quantum frequency
    nuclear_zq_frequency = 2 * b_12 * sqrt(1 + c_12^2)
    #println("Nuclear zero-quantum frequency: ", nuclear_zq_frequency)

    # compute pair contribution v_12 to the echo envelope 
    pair_contribution = - modulation_depth * (cos(nuclear_zq_frequency * time) - 1)^2
    #println("Pair contribution: ", pair_contribution)

    # compute echo envelope
    intensity = 1
    intensity = intensity * (1 + pair_contribution)
    #println("Simulated intensity: ", intensity)

    return modulation_depth,nuclear_zq_frequency,pair_contribution,intensity 
end


end
