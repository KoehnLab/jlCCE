module jlCCE

#push!(LOAD_PATH,".")

export SpinSystem, cce

using AtomsIO
using Unitful
using LinearAlgebra
using Printf

# import the readCIF.jl file to read cif files 
using readCIF
using SpinBase

# physical constants
# reduced Planck constant 
const hbar = 1.054571817e-27 # erg s -- J s -- CODATA 2022 (rounded value)
# Bohr magneton 
const mu_b = 9.274010066e-21 # erg G^-1 J T^-1 -- CODATA 2022 (rounded value)
# nuclear magneton
const mu_n = 5.050783739e-24 # erg G^-1  J T^-1 -- CODATA 2022 (rounded value)
# conversion AA to cm (cgs unit system)
const aacm = 1e-8

# struct to hold parameters for the run
mutable struct SpinSystem
    # file with coordinates (anything that the file read can handle)
    coord_file::String
    # name of the atom defining the spin center ("Cu", "V", ...)
    spin_center::String
    # index within unit cell to select spin center, if name is not unique
    spin_center_index::Int
    # type of simulation
    simulation_type::String
    use_exp::Bool
    do_cce1::Bool
    # spin quantum number
    s_el::Float64
    # g factor at the center (x,y,z)
    g_factor::Vector{Float64}
    # magnetic axes of spin center (column vectors for x, y, z direction)
    magnetic_axes::Matrix{Float64}
    # nuclei defining spin bath ("H", etc.)
    nuc_spin_bath::String
    # spin quantum number of nuclei
    s_nuc::Float64
    # corresponding nuclear g factor
    gn_spin_bath::Float64
    # magnetic field
    B0::Vector{Float64}
    # minimum interaction radius (usually 0.)
    r_min::Float64
    # maximum interaction radius
    r_max::Float64
    # maximum distance of two bath spins
    r_max_bath::Float64
    # minimum and maxiumum of time interval
    t_min::Float64
    t_max::Float64
    # number of time steps in interval
    n_time_step::Int
    # determine magnetic axes --> false for the simulation of the intensity
    det_mag_axes::Bool
end

# convenient constructor with defaults for all but the first 3 parameters
SpinSystem(coord_file,spin_center,spin_center_index) = SpinSystem(
    coord_file,spin_center,spin_center_index,
    "highfield_analytic",false,true,
    0.5,[2.0,2.0,2.0],[1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0],
    "H",0.5,5.58569468,
    [0.,0.,1.],0.0,50.0,100.0,0.0,1e-3,25,false)

"""
    cce(system::SpinSystem)

input: Spinsystem - cif file of the spin system, name of the metal spin center,
        index of the metal atom, g factor of the central spin, magnetic axes of 
        the spin center, nuclei defining spin center, corresponding nuclear g 
        factor, magnetic field, maximun und manimum radius around the spin center

simulate the Hahn echo intensity using the second-order cluster expantion 

Lit.: Witzel, W. M.; Das Sarma, S. Quantum theory for electron spin decoherence in-
duced by nuclear spin dynamics in semiconductor quantum computer architectures:
Spectral diffusion of localized electron spins in the nuclear solid-state environment.
Phys. Rev. B 2006, 74, 035322

"""
function cce(system::SpinSystem)

    println("")
    println("This is jlCCE")
    println("=============\n")

    println("Running on ",Threads.nthreads()," threads\n")
 
    if system.coord_file != "test" && system.coord_file != "test_rotated"
        # check cif file
        print("Coordinates will be read from: ",system.coord_file,"\n")

        # identify spin center
        if system.spin_center == "V"
            atomic_number_metal = 23
        elseif system.spin_center == "Cu"
            atomic_number_metal = 29
        elseif system.spin_center == "Pd"
            atomic_number_metal = 46
        else 
            print("Error currently only V, Cu and Pd \n")
            exit()
        end 
        print("\n")
        print("Current metal of the spin center: ",system.spin_center,"\n")

        # identify nuclear spin bath 
        if system.nuc_spin_bath == "H"
            atomic_number_nuclei = 1
        else 
            print("Error only proton bath \n")
            exit()
        end
        print("Considered nuclear spins of the spin bath: ",system.nuc_spin_bath,"\n")

        # get list of spin bath nuclei using the module readCIF

        # call function get_coordinates: determine lattice of the spin system, the coordinates of the 
            # electron spin center (x,y,z) and coordinates of the nuclear spins of the unit cell     
        lattice,coord_electron_spin,coords_nuclear_spins_unit_cell,coords_oxygen_unit_cell= 
            get_coordinates(system.coord_file,atomic_number_metal,system.spin_center_index,atomic_number_nuclei,false)

        println("\n")
        println("lattice vectors:")
        @printf " a = [%20.6f %20.6f  %20.6f]\n" lattice[1][1] lattice[1][2] lattice[1][3]
        @printf " b = [%20.6f %20.6f  %20.6f]\n" lattice[2][1] lattice[2][2] lattice[2][3]
        @printf " c = [%20.6f %20.6f  %20.6f]\n\n" lattice[3][1] lattice[3][2] lattice[3][3]
 
        println("Selected nucleus ",system.spin_center," no. ",system.spin_center_index," in unit cell")
        println("coordinates:")
        @printf " x  %20.6f Å\n" coord_electron_spin[1]
        @printf " y  %20.6f Å\n" coord_electron_spin[2]
        @printf " z  %20.6f Å\n\n" coord_electron_spin[3]

        # call function set_supercell_list: determine the cell list
        cell_list = set_supercell_list(system.r_max,lattice)

        n_cell = size(cell_list,1)
        println("r_max = ",system.r_max)
        println("Number of replicated unit cells for spin bath: ",n_cell,"\n")
        #print("cell list: \n")
        #print(cell_list,"\n")

        # call function get_bath_list: determine the distance coordinates between the electron spin center
        # and the nuclear spins and the number of considered nuclear spins in the spin bath
        distance_coordinates_el_nucs,n_nuc= 
            get_bath_list(system.r_min,system.r_max,lattice,coords_nuclear_spins_unit_cell,coord_electron_spin,coords_oxygen_unit_cell,system.det_mag_axes)
            #print("oxygen coords: \n")
	    #print(distance_coordinates_el_spin_oxygen)
	    
    else   # system.coords_file = test
    	println("Entered test mode")

    	if system.coord_file == "test"

            distance_coordinates_el_nucs = []

            # test - keep these (used for test suite)
            n_nuc = 2

	    I1 = [-10.,0.,10.]
            I2 = [10.,0.,20.]
	
	    push!(distance_coordinates_el_nucs,I1)
            push!(distance_coordinates_el_nucs,I2)

	elseif system.coord_file == "test_rotated" 

            # rotate the system
            #println("\n Investigation of rotated test system")
	   
	    distance_coordinates_el_nucs = []
            
	    n_nuc = 2

	    theta = 90
	    phi = 0

	    rot_mat_y = [cos(theta) 0 sin(theta); 0 1 0; -sin(theta) 0 cos(theta)]
	    rot_mat_z = [cos(phi) -sin(phi) 0; sin(phi) cos(phi) 0; 0 0 1]
	
	    I1 = [-10.,0.,10.]
	    I2 = [10.,0.,20.]

	    I1_rot = rot_mat_z * (rot_mat_y*I1)
	    I2_rot = rot_mat_z * (rot_mat_y*I2)

            push!(distance_coordinates_el_nucs,I1_rot)
            push!(distance_coordinates_el_nucs,I2_rot)
	end
    end

    println("Number of bath nuclei: ",n_nuc,"\n")   

    # rescale distance coordinates from AA to cm (cgs unit system)
    distance_coordinates_el_nucs = distance_coordinates_el_nucs .* aacm

    #print("Distance coordinates between the electron spin center and the nuclear spins: \n")
    #print(distance_coordinates_el_nucs,"\n") 
    
    # rotate the applied magnetic field
    #rotation_matrix_y,rotation_matrix_z = rotation_matrices(system.theta,system.phi)
    #B0 = system.magnetic_axes * (rotation_matrix_z * (rotation_matrix_y * system.B0)) 
    
    #println("Rotation of the inital magnetic field: ", system.B0)
     
    #println("Rotation around:")
    #@printf "theta (Y) %10.2f ° \n" rad2deg(system.theta)
    #@printf "phi (Z)   %12.2f ° \n" rad2deg(system.phi)
    #println("Rotation angle around Y axis (theta): ", rad2deg(system.theta))
    #println("Rotation angle around Z axis (phi): ", rad2deg(system.phi))
    #print("\n")

    println("Applied magnetic field:")
    @printf " x  %20.6f Gauss\n" system.B0[1]
    @printf " y  %20.6f Gauss\n" system.B0[2]
    @printf " z  %20.6f Gauss\n\n" system.B0[3]

    #print_matrix("Distances: ",distance_coordinates_el_nucs)
    
    # considering the anisotropy of the g factor --> determine an effective g factor
    Bnorm = system.B0/norm(system.B0)
    g_eff = Bnorm' * system.magnetic_axes * diagm(system.g_factor) * system.magnetic_axes' * Bnorm 

    println("Magnetic axes:")
    @printf "x [%10.6f %10.6f %10.6f] Å\n" system.magnetic_axes[1,1] system.magnetic_axes[2,1] system.magnetic_axes[3,1]
    @printf "y [%10.6f %10.6f %10.6f] Å\n" system.magnetic_axes[1,2] system.magnetic_axes[2,2] system.magnetic_axes[3,2]
    @printf "z [%10.6f %10.6f %10.6f] Å\n\n" system.magnetic_axes[1,3] system.magnetic_axes[2,3] system.magnetic_axes[3,3]

    println("g factor of ",system.coord_file,": ")
    @printf " x  %20.6f Å\n" system.g_factor[1]
    @printf " y  %20.6f Å\n" system.g_factor[2]
    @printf " z  %20.6f Å\n\n" system.g_factor[3]

    println("Determined effective g factor: ",g_eff)
 
    # calculation of the gryomagnetic ratios of the central electron spin center and the nucle ar spins of the spin bath
    gamma_electron = g_eff .* (mu_b / hbar)
    gamma_n = (system.gn_spin_bath * mu_n) / hbar

    # set time for the simulation
    time_hahn_echo = collect(range(system.t_min,system.t_max,system.n_time_step))

    if system.simulation_type == "highfield_analytic"

        # test settings
        err = false
        if system.s_el != 0.5 || system.s_nuc != 0.5
            print("Error: Analytic approach is only for spin-1/2 nuclei\n")
            print("   Found: s_el = ",system.s_el,"  s_nuc = ",system.s_nuc,"\n")
            err = true
        end
        if norm(system.B0)<1e-10
            print("Error: Too small field applied (zero field?)\n")
            print("   Found: B0 = ",system.B0,"\n")
            err = true
        end
        #if system.g_factor[1] != system.g_factor[2] || system.g_factor[1] != system.g_factor[3]
        #    print("Error: Analytic approach only for isotropic spin\n")
        #    print("    Found: g_factor = ",system.g_factor,"\n")
        #    err = true
        #end
        if err
            error("Unsuitable settings found (see above)")
        end

        intensity = cce_hf_analytic(distance_coordinates_el_nucs,n_nuc,system.r_max_bath*aacm,gamma_n,gamma_electron[1],system.B0,time_hahn_echo,system.use_exp)
        # dummies for iCCE1 and iCCE2
        iCCE1 = ones(size(time_hahn_echo))
        iCCE2 = intensity
    elseif system.simulation_type == "exact"
        intensity,iCCE1,iCCE2 = cce_exact(distance_coordinates_el_nucs,n_nuc,system.r_max_bath*aacm,
             system.s_nuc,gamma_n,system.s_el,gamma_electron,system.magnetic_axes,system.B0,system.do_cce1,time_hahn_echo)
    else
        print("Unkonwn simulation type: ",system.simulation_type,"\n")
        intensity,iCCE1,iCCE2 = zeros(size(time_hahn_echo))
    end

    # return intensity and time
    return time_hahn_echo,intensity,iCCE1,iCCE2
end

"""
    cce_hf_analytic(distance_coordinates_el_nucs,n_nuc,r_max_bath,gamma_n,gamma_electron,B0,time_hahn_echo)

should be called by cce (rather than directly); distance_coordinates_el_nucs are the distances 
of the bath nuclei from the spin center, gamma_n and gamma_electron the respective gyromagnetic 
rations required for computing the dipolar interaction, B0 is the magnetic field direction,
time_hahn_echo the time set for which the signal is to be computed

the function returns the intensity for the given time set
"""
function cce_hf_analytic(distance_coordinates_el_nucs,n_nuc,r_max_bath,gamma_n,gamma_electron,B0,time_hahn_echo,use_exp)

    n_time_step = size(time_hahn_echo)

    # precompute values of the hyperfine coupling constant A for electron nucleus pairs
    A_n = zeros(n_nuc)
    #print("Distance between el spin and nuc spins: \n") 
    for i in 1:size(distance_coordinates_el_nucs)[1]
        r_i_x_B0 = cross(distance_coordinates_el_nucs[i], B0)
        r_i_dot_B0 = dot(distance_coordinates_el_nucs[i], B0)
        theta_i = atan(norm(r_i_x_B0), r_i_dot_B0)
        r_i_norm = norm(distance_coordinates_el_nucs[i])
        #print(r_i_norm,"\n") 
        #print(gamma_n,gamma_electron,hbar,theta_i,r_i_norm)
        A_n[i] = -gamma_n * gamma_electron * hbar * (1 - 3 * cos(theta_i)^2) / r_i_norm^3
    end
    minA = minimum(abs.(A_n))
    maxA = maximum(abs.(A_n))
    print("\n")
    print("Hyperfine coupling constants (min, max) in Hz: ",minA,"  ",maxA,"\n")
    print("\n")
    
    # not directly used below, but call this to report the screening
    n_pairs, pair_list, n_pair_contr = make_pair_list(distance_coordinates_el_nucs,r_max_bath)

    print("Total number of pairs in bath: ",n_nuc*(n_nuc-1)÷2,"\n")
    print("Screened number of pairs:      ",n_pairs,"\n")
    print("Screening distance was: ",r_max_bath/aacm," Å\n\n")

    # initialize Intensity to 1. for all times
    intensity = ones(n_time_step)
    
    # compute the coupling constant b and the values of c and omega for each pair of bath nuclei 
    # compute the pair intensity for all t values of the simulation for each pair of bath nuclei 
    # compute the intensity for each t value using the pair intensities
    for n in 1:n_nuc-1
        for m in n+1:n_nuc
            r_nm = (distance_coordinates_el_nucs[m] - distance_coordinates_el_nucs[n]) 
            if norm(r_nm) > r_max_bath
                continue
            end
            r_nm_x_B0 = cross(r_nm, B0)
            r_nm_dot_B0 = dot(r_nm, B0)
            theta_nm = atan(norm(r_nm_x_B0), r_nm_dot_B0)
            b_nm = -0.25 * gamma_n^2 * hbar * (1 - 3 * cos(theta_nm)^2) / norm(r_nm)^3
            c_nm = (A_n[n] - A_n[m]) / (4. * b_nm)
            w_nm = 2. * b_nm * sqrt(1 + c_nm^2)
            if use_exp
                for j in 1:size(time_hahn_echo)[1]
                    v_nm = -((c_nm^2) / (1 + c_nm^2)^2) * (cos(w_nm * time_hahn_echo[j]) - 1)^2
                    intensity[j] = intensity[j] * exp(v_nm)
                end
            else
                for j in 1:size(time_hahn_echo)[1]
                    v_nm = -((c_nm^2) / (1 + c_nm^2)^2) * (cos(w_nm * time_hahn_echo[j]) - 1)^2
                    intensity[j] = intensity[j] * (1 + v_nm)
                end
            end 
        end
    end        
     
    return intensity
end

"""
    hyperfine(gamma_ten,gamma_n,r)

compute the hyperfine interaction tensor (dipole approximation) based on the gyromagnetic 
tensor of the electron (gamma_ten = mu_b/hbar * g_tensor) and the (isotropic) nuclear 
magnetic moment gamma_n. r is the distance.
NOTE: non isotropic electron moments are not yet debugged!
"""
function hyperfine(gamma_ten,gamma_n,r)


    dist = norm(r)
    rn = r ./ dist
    fact = gamma_n*hbar / dist^3
    
    A = fact*(I(3) - 3 * rn*transpose(rn))

    At = transpose(gamma_ten) * A 

    return At 
end

"""
    print_matrix(title,mat)

print a title and the matrix (using display function)
"""
function print_matrix(title,mat)

    print(title,"\n")
    display(mat)
    print("\n")

end


"""
    n_e_contribution(various args)

core routine for "exact" propagation of 2-particle system (internal use only)
"""
function n_e_contribution(dim_el,dim_nuc,gamma_n,r12,A_1,Hmag,Mmat,Smat,Imat,rho0,time_he)

    dim = dim_nuc*dim_el
    nt = size(time_he)[1]

    # spin H for 2 spin system
    S1 = unit_mat(dim_el)
    I1 = unit_mat(dim_nuc)

    # interaction contribution
    Hmat = zeros(dim,dim)
    for i in 1:3
        for j in 1:3
            Hmat += kron(Smat[i],Imat[j]).*A_1[i,j]
        end
    end

    #print_matrix("Hmat:",Hmat)

    # add magnetic contribution
    Hmat += Hmag

    #print_matrix("Hmag:",Hmag)

    # Sx operator for spin (emulating ideal π pulse)
    sigmaX = kron(Smat[1],I1)*2

    mx = tr(Mmat[1]*rho0)
    my = tr(Mmat[2]*rho0)

    #print_matrix("rho0",rho0)

    # initial signal
    signal_0 = 2*real(mx+1im*my)

    signal = zeros(nt)
    for idx in 1:nt

        Arg = - (1im * time_he[idx]) * Hmat
        Ut = exp(Arg)

        U = Ut * sigmaX * Ut

        rho_tau = U * rho0 * adjoint(U)

        mx = tr(Mmat[1]*rho_tau)
        my = tr(Mmat[2]*rho_tau)

        signal[idx] = 2*real(mx+1im*my)/signal_0
        
    end

    return signal

end



"""
    n_n_e_contribution(various args)

core routine for "exact" propagation of 3-particle system (internal use only)
"""
function n_n_e_contribution(dim_el,dim_nuc,gamma_n,r12,A_1,A_2,Hmag,Mmat,Smat,Imat,rho0,time_he)

    dim = dim_nuc^2*dim_el
    nt = size(time_he)[1]

    # spin H for 3 spin system
    S1 = unit_mat(dim_el)
    I1 = unit_mat(dim_nuc)

    # spin-spin coupling
    dist = norm(r12)
    rn = r12 ./ dist
    fact = gamma_n^2*hbar/dist^3

    B = fact * (I(3) - 3 * rn*transpose(rn))

    #print_matrix("A1",A_1)
    #print_matrix("A2",A_2)
    #print_matrix("B",B)

    # interaction contribution
    Hmat = zeros(dim,dim)
    for i in 1:3
        for j in 1:3
            Hmat += kron(Smat[i],kron(Imat[j],I1)).*A_1[i,j]
            Hmat += kron(Smat[i],kron(I1,Imat[j])).*A_2[i,j]
            Hmat += kron(S1,kron(Imat[i],Imat[j])).*B[i,j]
        end
    end

    #print_matrix("Hmat:",Hmat)

    # add magnetic contribution
    Hmat += Hmag

    #print_matrix("Hmag:",Hmag)

    # Sx operator for spin (emulating ideal π pulse)
    sigmaX = kron(Smat[1],kron(I1,I1))*2

    mx = tr(Mmat[1]*rho0)
    my = tr(Mmat[2]*rho0)

    #print_matrix("rho0",rho0)

    # initial signal
    signal_0 = 2*real(mx+1im*my)

    signal = zeros(nt)
    for idx in 1:nt

        Arg = - (1im * time_he[idx]) * Hmat
        Ut = exp(Arg)

        U = Ut * sigmaX * Ut

        rho_tau = U * rho0 * adjoint(U)

        mx = tr(Mmat[1]*rho_tau)
        my = tr(Mmat[2]*rho_tau)

        signal[idx] = 2*real(mx+1im*my)/signal_0
        
    end

    return signal

end


function is_unit(mat,thr)

    nd = ndims(mat)
    if nd != 2 
        return false
    end
    nd1 = size(mat,1)
    nd2 = size(mat,2)
    if nd1 != nd2
        return false
    end
    OK = true
    for i in 1:nd1
        for j in 1:nd2
            if i==j && abs(mat[i,i]-1.0)>thr
                OK = false
            end
            if i!=j && abs(mat[i,j])>thr
                OK = false
            end
        end
    end
    return OK
end

"""
    get_rot_for_unit_vector(uvec)

returns a rotation matrix that transforms vector (0,0,1) into uvec
"""
function get_rot_for_unit_vector(uvec)

    # get the direction cosine and sines
    costh = min(1.,max(-1.,uvec[3]))
    sinth = sqrt(1. - costh^2)
    if abs(sinth)>1e-9
        cosphi = min(1.,max(-1.,uvec[1]/sinth))
        sinphi = sqrt(1. - cosphi^2)
        if uvec[2] < 0.
            sinphi *= -1
        end
    else
        cosphi=1.
        sinphi=0.
    end

    # rotation around y by theta
    yrot = [[costh,0.,sinth] [0.,1.,0.] [-sinth,0.,costh]]
    # rotation around z by -phi
    zrot = [[cosphi,-sinphi,0.] [sinphi,cosphi,0.] [0.,0.,1.]]

    # total rotation matrix
    Rot = zrot * yrot

    return Rot

end

"""
    make_pair_list(distance_coordinates_el_nucs,r_max_bath)

make a list of all pairs of bath spins that are considered, where r_max_bath is the
max. distance to be considered.
Return a list of unique pairs and a list which states to how many pairs a given spin 
contributes 
"""
function make_pair_list(distance_coordinates_el_nucs,r_max_bath)

    n_nuc = size(distance_coordinates_el_nucs,1)

    n_pairs = 0
    pair_list = []
    n_pair_contr = zeros(Int, n_nuc)

    for i = 1:n_nuc-1
        for j = i+1:n_nuc
            r12 = distance_coordinates_el_nucs[i] - distance_coordinates_el_nucs[j]
            if norm(r12) > r_max_bath
                continue
            end

            n_pairs += 1

            push!(pair_list,[i,j])

            n_pair_contr[i] += 1
            n_pair_contr[j] += 1

        end
    end

    return n_pairs,pair_list,n_pair_contr

end


function cce_exact(distance_coordinates_el_nucs,n_nuc,r_max_bath,
    s_nuc,gamma_n,s_el,gamma_el,mag_axes,B00,do_cce1,time_hahn_echo)

    if ! is_unit(mag_axes,1e-8)
        print("non unit magnetic axes not yet debugged!")
        error("non unit magentic axes not yet debugged!")
    end

    n_time_step = size(time_hahn_echo)

    intensity = ones(n_time_step)
    intensity_CCE1 = ones(n_time_step)
    intensity_CCE2 = ones(n_time_step)
    intensity_correction = ones(n_time_step)

    # Here, it is more convenient to keep the magnetic field at (0,0,1) and to
    # transform the nuclear coordinates
    
    # field strength
    Bstrength = norm(B00)

    Rmat = get_rot_for_unit_vector(B00/Bstrength)

    print("Found this field: ",B00,"\n")

    print("Transforming by:\n")
    display(Rmat)
    print("\n\n")

    n_pairs, pair_list, n_pair_contr = make_pair_list(distance_coordinates_el_nucs,r_max_bath)

    print("Total number of pairs in bath: ",n_nuc*(n_nuc-1)÷2,"\n")
    print("Screened number of CCE2 pairs: ",n_pairs,"\n")
    print("Screening distance was: ",r_max_bath/aacm," Å\n")

    B0 = [0.,0.,Bstrength]

    # transform coordinates
    distance_el_nuc_traf = []
    for coord in distance_coordinates_el_nucs
        coord_traf = transpose(Rmat) * coord
        push!(distance_el_nuc_traf,coord_traf)
    end

    mag_axes_t = transpose(Rmat) * mag_axes

    # "g tensor" from gamma and magnetic axes
    gamma_t = zeros(3,3)
    gamma_t[1,1] = gamma_el[1]; gamma_t[2,2] = gamma_el[2]; gamma_t[3,3] = gamma_el[3] 
    gamma_t = mag_axes_t * gamma_t * transpose(mag_axes_t)

    print("gamma_t: ",gamma_t,"\n")

    print("gamma_n: ",gamma_n,"\n")
    
    print("\n")
    print("computing hyperfine couplings ...\n")

    A_list = []
    for n in 1:n_nuc
        r1 = distance_el_nuc_traf[n]
        A_current = hyperfine(gamma_t,gamma_n,r1)
        push!(A_list,A_current)
    end

    dim_el = trunc(Int,(2*s_el+1))
    dim_nuc = trunc(Int,(2*s_nuc+1))

    Smat = []
    push!(Smat,sx_mat(s_el))
    push!(Smat,sy_mat(s_el))
    push!(Smat,sz_mat(s_el))
    S1 = unit_mat(dim_el)

    Imat = []
    push!(Imat,sx_mat(s_nuc))
    push!(Imat,sy_mat(s_nuc))
    push!(Imat,sz_mat(s_nuc))
    I1 = unit_mat(dim_nuc)

    # initial spin state |+><+| (max. comp along X)
    if s_el == 0.5
        rho0_el = ones(2,2).*0.5
    else
        error("S>0.5 not yet implemented")
        # TODO: diagonalize Sx and get eigenvector for -M state == |+> ; rho = |+><+|
    end
        
    # initial nuclear state: all M equal likely (kT >> delta E)
    rho0_nuc = I1.*(1. / trunc(Int,2*s_nuc+1) )
    
    if do_cce1
        # CCE1 contrubutions
        print("\n")
        print("Computing CCE1 contributions...\n")

        dimH1 = dim_el * dim_nuc

        # set magnetic interaction matrix
        Mmat1 = []
        for i = 1:3
            mat = zeros(dimH1,dimH1)
            # needs to be updated for non-isotropic g
            mat += kron(Smat[i],I1).*gamma_t[i,i]
            mat += kron(S1,Imat[i]).*gamma_n
            push!(Mmat1,mat)
        end

        # magnetic field contribution is universal
        Hmag1 = zeros(dimH1,dimH1)
        B0t = gamma_t * B0
        for i = 1:3
            Hmag1 += kron(Smat[i],I1).*B0t[i]
            Hmag1 += kron(S1,Imat[i]).*(B0[i]*gamma_n)
        end

        # set initial density matrix for 2-particle system    
        rho0 = kron(rho0_el,rho0_nuc)*(1. + 0*im)

        for n in 1:n_nuc
            r1 = distance_el_nuc_traf[n]
            A_n = A_list[n]
            ne_cont = n_e_contribution(dim_el,dim_nuc,gamma_n,r1,A_n,Hmag1,Mmat1,Smat,Imat,rho0,time_hahn_echo)

            mult = n_pair_contr[n]
            intensity_CCE1 = intensity_CCE1 .* ne_cont
            intensity_correction = intensity_correction .* ne_cont.^mult

        end
    end

    print("\n")
    print("Computing CCE2 contributions ...\n")

    dimH = dim_el * dim_nuc^2

    # Magnetic interaction matrix
    Mmat = []
    for i = 1:3
        mat = zeros(dimH,dimH)
        # needs to be updated for non-isotropic g
        mat += kron(Smat[i],kron(I1,I1)).*gamma_t[i,i]
        mat += kron(S1,kron(Imat[i],I1)).*gamma_n
        mat += kron(S1,kron(I1,Imat[i])).*gamma_n
        push!(Mmat,mat)
    end

    # magnetic field contribution is universal
    Hmag = zeros(dimH,dimH)
    B0t = gamma_t * B0
    for i = 1:3
        Hmag += kron(Smat[i],kron(I1,I1)).*B0t[i]
        Hmag += kron(S1,kron(Imat[i],I1)).*(B0[i]*gamma_n)
        Hmag += kron(S1,kron(I1,Imat[i])).*(B0[i]*gamma_n)
    end

    # set initial density matrix for 3-particle system    
    rho0 = kron(rho0_el,kron(rho0_nuc,rho0_nuc))*(1. + 0*im)

    if Threads.nthreads() == 1 

        for ipair = 1:n_pairs

            n = pair_list[ipair][1]
            m = pair_list[ipair][2]

            r1 = distance_el_nuc_traf[n]
            r2 = distance_el_nuc_traf[m]

            A_n = A_list[n]
            A_m = A_list[m]
            nne_cont = n_n_e_contribution(dim_el,dim_nuc,gamma_n,r1-r2,A_n,A_m,Hmag,Mmat,Smat,Imat,rho0,time_hahn_echo)

            intensity_CCE2 = intensity_CCE2 .* nne_cont
        end

    else

        print("entered threaded route")
        
        # block the pair list among the threads
        blocks = Iterators.partition(1:n_pairs, n_pairs ÷ Threads.nthreads())
        print(blocks,"\n")
        tasks = map(blocks) do block

            Threads.@spawn begin
                int_block = ones(n_time_step)
                for ipair in block
                    n = pair_list[ipair][1]
                    m = pair_list[ipair][2]

                    r1 = distance_el_nuc_traf[n]
                    r2 = distance_el_nuc_traf[m]
                    A_n = A_list[n]
                    A_m = A_list[m]
                    nne_cont = n_n_e_contribution(dim_el,dim_nuc,gamma_n,r1-r2,A_n,A_m,Hmag,Mmat,Smat,Imat,rho0,time_hahn_echo)

                    int_block = int_block .* nne_cont
                end
                int_block
            end
        end

        # get all the results ...
        int_blocks = fetch.(tasks)
        # ... and accumulate them
        for int in int_blocks
            intensity_CCE2 = intensity_CCE2 .* int
        end

    end

    intensity_CCE2 = intensity_CCE2 ./ intensity_correction

    intensity = intensity_CCE1 .* intensity_CCE2

    return intensity, intensity_CCE1, intensity_CCE2

end

end
