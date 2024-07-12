begin
using LinearAlgebra
using Printf
using DelimitedFiles
using Serialization
using Plots
using MKL
using BenchmarkTools
using CSV
using DataFrames
end 

begin
function parameters()
    nuc::String = "H"
    B0::Vector{Float64} = [0.0, 0.0, 1.0]
    t::Vector{Float64} = collect(range(0, 1.5*10^-5, length=30))
    radius_min::Float64 = 0.0
    radius_max::Float64 = 5.0
    hbar::Float64 = 1.054571817*10^-27
    mu_n::Float64 = 5.0507866*10^-24
    g_n::Float64 = 5.585694680000000
    gamma_n::Float64 = (g_n * mu_n) / hbar
    mu_b::Float64 = 9.27400968*10^-21
    g_e::Float64 = 2.25
    gamma_electron::Float64 = (g_e * mu_b) / hbar
    return nuc, B0, t, radius_min, radius_max, hbar, gamma_n, gamma_electron
end

function def_molecule()
    molecule = "vmnt3"
    if molecule == "pddbm2"
        # xyzfile = "/home/suchaneck/masterarbeit/molecular_structures/pddbm2.xyz"
        xyzfile = "C:\\Users\\sucha\\MATLAB Drive\\molecular_structures\\pddbm2.xyz"
        coord_electron = [14.8993 13.0050 29.5112]
        metal = "Pd"
    elseif molecule == "vodpm2"
        # xyzfile = "/home/suchaneck/masterarbeit/molecular_structures/vodpm2.xyz"
        xyzfile = "C:\\Users\\sucha\\MATLAB Drive\\molecular_structures\\vodpm2.xyz"
        coord_electron = [16.6316 28.7243 28.0580]
        metal = "V"
    elseif molecule == "vodbm2"
        # xyzfile = "/home/suchaneck/masterarbeit/molecular_structures/vodbm2.xyz"
        xyzfile = "C:\\Users\\sucha\\MATLAB Drive\\molecular_structures\\vodbm2.xyz"
        coord_electron = [30.5801 48.0953 30.6971]
        metal = "V"
    elseif molecule == "voacac2"
        # xyzfile = "/home/suchaneck/masterarbeit/molecular_structures/voacac2.xyz"
        xyzfile = "C:\\Users\\sucha\\MATLAB Drive\\molecular_structures\\voacac2.xyz"
        coord_electron = [47.8550 32.1037 39.4976]
        metal = "V"
    elseif molecule == "vdmit3"
        # xyzfile = "/home/suchaneck/masterarbeit/molecular_structures/vdmit3.xyz"
        xyzfile = "C:\\Users\\sucha\\MATLAB Drive\\molecular_structures\\vdmit3.xyz"
        coord_electron = [42.3870 29.3239 34.4058]
        metal = "V"
    elseif molecule == "vodmit2"
        # xyzfile = "/home/suchaneck/masterarbeit/molecular_structures/vodmit2.xyz"
        xyzfile = "C:\\Users\\sucha\\MATLAB Drive\\molecular_structures\\vodmit2.xyz"
        coord_electron = [37.6208 24.4527 35.9016]
        metal = "V"
    elseif molecule == "vdbddto3"
        # xyzfile = "/home/suchaneck/masterarbeit/molecular_structures/vdbddto3.xyz"
        xyzfile = "C:\\Users\\sucha\\MATLAB Drive\\molecular_structures\\vdbddto3.xyz"
        coord_electron = [15.2913 21.2973 38.0247]
        metal = "V"
    elseif molecule == "vmnt3"
        # xyzfile = "/home/suchaneck/masterarbeit/molecular_structures/vmnt3.xyz"
        xyzfile = "C:\\Users\\sucha\\MATLAB Drive\\molecular_structures\\vmnt3.xyz"
        coord_electron = [39.7660 33.4080 31.8062]
        metal = "V"
    else molecule == "vopc"
        # xyzfile = "/home/suchaneck/masterarbeit/molecular_structures/vopc.xyz"
        xyzfile = "C:\\Users\\sucha\\MATLAB Drive\\molecular_structures\\vopc.xyz"
        coord_electron = [45.1115 32.6557 27.1924]
        metal = "V"
    end
    return xyzfile, coord_electron, metal
end


function read_geom(xyzfile::String,nuc::String)
    coords = readdlm(xyzfile,skipstart=2)[:,2:4]
    coords = convert(Array{Float64,2},coords)
    atoms = readdlm(xyzfile,String,skipstart=2)[:,1]
    atoms_nuc = findall(x -> x == nuc, atoms)
    coords_nuc::Matrix{Float64} = coords[atoms_nuc,:]
    return coords_nuc
end

    
function calculate_r_e(coords_nuc::Array{Float64},coord_electron::Array{Float64})
    # calculation of the distances of the electron and the nuclei
    r_e::Matrix{Float64} = zeros(Float64,size(coords_nuc))
    r_e_norm::Vector{Float64} = zeros(Float64,size(coords_nuc)[1])
    for i in 1:size(coords_nuc)[1]
         r_e[i,:] = coords_nuc[i,:] - coord_electron[:]
         r_e_norm[i] = norm(r_e[i,:]) 
    end
    return r_e_norm 
end

function determine_coords_nuc_restricted(radius_min::Float64,radius_max::Float64,r_e_norm::Vector{Float64},coords_nuc::Matrix{Float64})
    # nuclei within a restricted radius
    coords_nuc_restricted = Matrix{Float64}(undef,0,3)
    for i in 1:size(coords_nuc)[1]
        if r_e_norm[i] < radius_max && r_e_norm[i] > radius_min
             coords_nuc_restricted = vcat(coords_nuc_restricted, transpose(coords_nuc[i,:]))
        end
    end
    return coords_nuc_restricted
end

function calculate_r_i(coords_nuc_restricted::Array{Float64},coord_electron::Array{Float64})
    r_i::Matrix{Float64} = zeros(Float64,size(coords_nuc_restricted))
    for i in 1:size(coords_nuc_restricted)[1]
         r_i[i,:] = (coords_nuc_restricted[i,:] - coord_electron[:]) .* 10^-8
    end
    return r_i
end

function calculate_A_n(coords_nuc_restricted::Array{Float64},r_i::Array{Float64},B0::Array{Float64},gamma_n::Float64, gamma_electron::Float64,hbar::Float64) 
    # calculation of the hyperfine coupling constant A_n
    A_n::Vector{Float64} = zeros(size(coords_nuc_restricted)[1])
    r_i_x_B0::Matrix{Float64} = zeros(size(coords_nuc_restricted))
    r_i_dot_B0::Vector{Float64} = zeros(size(coords_nuc_restricted)[1])
    theta_i::Vector{Float64} = zeros(size(coords_nuc_restricted)[1])
    r_i_norm::Vector{Float64} = zeros(Float64,size(coords_nuc_restricted)[1])
    for i in 1:size(coords_nuc_restricted)[1]
         r_i_x_B0[i, :] = cross(r_i[i, :], B0[:])
         r_i_dot_B0[i] = dot(r_i[i, :], B0[:])
         theta_i[i] = atan(norm(r_i_x_B0[i, :]), r_i_dot_B0[i])
         r_i_norm[i] = norm(r_i[i,:]) 
         A_n[i] = -gamma_n * gamma_electron * hbar * (1 - 3 * cos(theta_i[i])^2) / r_i_norm[i]^3
    end
    return A_n
end


function calculate_b_nm(coords_nuc_restricted,gamma_n,B0,hbar)
    # calculation of b_nm
    b_nm::Matrix{Float64} = zeros(size(coords_nuc_restricted)[1], size(coords_nuc_restricted)[1])
    for n in 1:size(coords_nuc_restricted)[1]-1
        for m in n+1:size(coords_nuc_restricted)[1]
            r_nm = (coords_nuc_restricted[m, :] .- coords_nuc_restricted[n, :]) * 10^-8
            r_nm_x_B0 = cross(r_nm, B0)
            r_nm_dot_B0 = dot(r_nm, B0)
            theta_nm = atan(norm(r_nm_x_B0), r_nm_dot_B0)
            b_nm[n, m] = -0.25 .* gamma_n.^2 .* hbar .* (1 .- 3 .* cos(theta_nm).^2) ./ norm(r_nm).^3
        end
    end
    return b_nm
end


function calculate_v_nm_sim(A_n::Vector{Float64},t::Vector{Float64},b_nm::Matrix{Float64}) # Most expensive part
    v_nm::Matrix{Float64} = zeros(size(A_n)[1], size(A_n)[1])
    sim::Vector{Float64} = zeros(size(t)[1])
    c_nm::Float64 = 0.0
    w_nm::Float64 = 0.0
    for j in 1:size(t)[1]
        for n in 1:size(A_n)[1]-1
            for m in n+1:size(A_n)[1]
                c_nm = (A_n[n] - A_n[m]) / (4 * b_nm[n, m])
                w_nm = 2 * b_nm[n, m] * sqrt(1 + c_nm.^2)
                v_nm[n, m] = -((c_nm^2) / (1 + c_nm^2)^2) * (cos(w_nm * t[j]) - 1)^2
            end
        end
        sim[j] = sum(v_nm) 
    end
    return sim
end

function calculate_intensity(sim)
    # calculation of the intensities
    intensity = exp.(sim)
    return intensity
end

end

function main()
    nuc::String, B0::Vector{Float64}, t::Vector{Float64}, radius_min::Float64, radius_max::Float64, hbar::Float64, gamma_n::Float64, gamma_electron::Float64 = parameters()
    xyzfile::String, coord_electron::Matrix{Float64}, metal::String = def_molecule()
    coords_nuc::Matrix{Float64} = read_geom(xyzfile, nuc)
    r_e_norm::Vector{Float64} = calculate_r_e(coords_nuc,coord_electron)
    coords_nuc_restricted::Matrix{Float64} = determine_coords_nuc_restricted(radius_min,radius_max,r_e_norm,coords_nuc)
    display("Number of nuclei:")
    display(size(coords_nuc_restricted)[1])
    r_i::Matrix{Float64} = calculate_r_i(coords_nuc_restricted,coord_electron)
    A_n::Vector{Float64} = calculate_A_n(coords_nuc_restricted,r_i,B0,gamma_n,gamma_electron,hbar)
    b_nm::Matrix{Float64} = calculate_b_nm(coords_nuc_restricted,gamma_n,B0,hbar)
    sim::Vector{Float64} = calculate_v_nm_sim(A_n,t,b_nm)
    intensity::Vector{Float64} = calculate_intensity(sim)
    return intensity,t
end

@time intensity, t = main();

# safe data as csv file
# writedlm("/home/suchaneck/masterarbeit/decoherence_decay_curves/int_vodpm2_30.csv",intensity)
# writedlm("C:/Users/sucha/MATLAB Drive/decoherence_decay_curves/pddbm2/sim_int_pddbm2_radius_15_julia.csv",intensity)

# plot(t*2*10^6,intensity)
# plot!(xlabel="time 2t", ylabel="intensity")




