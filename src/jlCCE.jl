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
    # nuclei defining spin bath ("H", etc.)
    nuc_spin_bath::String
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
    set_custom_times::Bool
    times::Vector{Float64}
    report_pair_contrib::Bool
    pair_log_file::String
    # select a range of nuclei
    select::Bool
    rselect::Vector{Float64}
end

# convenient constructor with defaults for all but the first 3 parameters
SpinSystem(coord_file,spin_center,spin_center_index) = SpinSystem(
    coord_file,spin_center,spin_center_index,"H",
    "highfield_analytic",false,true,
    0.5,[2.0,2.0,2.0],[1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0],
    0.5,5.58569468,
    [0.,0.,1.],0.0,50.0,100.0,0.0,1e-3,25,false,[5e-6],false,"pair_log.txt",
    false,[0.0,50.0,0.0,100.0])

SpinSystem(coord_file,spin_center,spin_center_index,nuc_spin_bath) = SpinSystem(
    coord_file,spin_center,spin_center_index,nuc_spin_bath,
    "highfield_analytic",false,true,
    0.5,[2.0,2.0,2.0],[1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0],
    0.5,5.58569468,
    [0.,0.,1.],0.0,50.0,100.0,0.0,1e-3,25,false,[5e-6],false,"pair_log.txt",
    false,[0.0,50.0,0.0,100.0])


"""
    cce(system::SpinSystem)

input: Spinsystem - cif file of the spin system, name of the metal spin center,
        index of the metal atom, g factor of the central spin, magnetic axes of 
        the spin center, nuclei defining spin center, corresponding nuclear g 
        factor, magnetic field, maximun und manimum radius around the spin center

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

        atomic_number_metal = symbol_to_atomic_number(system.spin_center)
        print("\n")
        print("Atom symbol for selecting spin center: ",system.spin_center,"\n")

        # identify nuclear spin bath 
        if system.nuc_spin_bath == "H"
            atomic_number_nuclei = 1
        else 
            print("Error only proton bath \n")
            exit()
        end
        print("Atom symbols for selecting spin bath: ",system.nuc_spin_bath,"\n")

        # get list of spin bath nuclei using the module readCIF

        # call function get_coordinates: determine lattice of the spin system, the coordinates of the 
        # electron spin center (x,y,z) and coordinates of the nuclear spins of the unit cell     
        lattice,coord_electron_spin = 
            get_spin_center(system.coord_file,atomic_number_metal,system.spin_center_index)
        
        coords_nuclear_spins_unit_cell = 
            get_coordinates_nuclear_spins(system.coord_file,atomic_number_nuclei)
        
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
            get_bath_list(system.r_min,system.r_max,lattice,coords_nuclear_spins_unit_cell,coord_electron_spin)
            
	    
    else   # system.coords_file = test
    	println("Entered test mode")

    	if system.coord_file == "test"

            distance_coordinates_el_nucs = []

            # test - keep these (used for test suite)
            n_nuc = 2

	        I1 = [20.,0.,0.]
            I2 = [23.,0.,0.]
	
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
	
	        I1 = [20.,0.,0.]
            I2 = [23.,0.,0.]

	        I1_rot = rot_mat_z * (rot_mat_y*I1)
	        I2_rot = rot_mat_z * (rot_mat_y*I2)

            push!(distance_coordinates_el_nucs,I1_rot)
            push!(distance_coordinates_el_nucs,I2_rot)
	    end
    end

    println("Number of bath nuclei: ",n_nuc,"\n")   

    # rescale distance coordinates from AA to cm (cgs unit system)
    distance_coordinates_el_nucs = distance_coordinates_el_nucs .* aacm

    println("Applied magnetic field:")
    @printf " x  %20.6f Gauss\n" system.B0[1]
    @printf " y  %20.6f Gauss\n" system.B0[2]
    @printf " z  %20.6f Gauss\n\n" system.B0[3]

    
    println("Magnetic axes:")
    @printf " x [%10.6f %10.6f %10.6f] Å\n" system.magnetic_axes[1,1] system.magnetic_axes[2,1] system.magnetic_axes[3,1]
    @printf " y [%10.6f %10.6f %10.6f] Å\n" system.magnetic_axes[1,2] system.magnetic_axes[2,2] system.magnetic_axes[3,2]
    @printf " z [%10.6f %10.6f %10.6f] Å\n\n" system.magnetic_axes[1,3] system.magnetic_axes[2,3] system.magnetic_axes[3,3]

    println("g factor: ")
    @printf " x  %10.6f \n" system.g_factor[1]
    @printf " y  %10.6f \n" system.g_factor[2]
    @printf " z  %10.6f \n\n" system.g_factor[3]

    if system.simulation_type == "highfield_analytic"
        # considering the anisotropy of the g factor --> determine an effective g factor
        Bnorm = system.B0/norm(system.B0)
        g_eff = Bnorm' * system.magnetic_axes * diagm(system.g_factor) * system.magnetic_axes' * Bnorm 

        @printf "Determined effective g factor: %10.6f  (along magnetic field, relevant for high-field approximation)\n" g_eff
        # calculation of the gryomagnetic ratios of the central electron spin center and the nucle ar spins of the spin bath
        gamma_electron_sc = g_eff .* (mu_b / hbar)
    else
        gamma_electron = system.g_factor .*  (mu_b / hbar)
    end

    gamma_n = (system.gn_spin_bath * mu_n) / hbar

    @printf " Nuclear g factor:             %10.6f \n\n" system.gn_spin_bath 

    # set time for the simulation (corresponds to delay time tau)
    if ! system.set_custom_times
        @printf "Setting delay times (tau):\n"
        @printf "t_min = %12.3f µs\n" system.t_min*1e6
        @printf "t_max = %12.3f µs\n" system.t_max*1e6
        @printf "number of time steps: %6i\n" system.n_time_step
        time_hahn_echo = collect(range(system.t_min,system.t_max,system.n_time_step))
        if system.n_time_step < 4
            @printf "Times: ["
            for idx in 1:system.n_time_step-1
                @printf " %12.3f, " time_hahn_echo[idx]*1e6
            end
            @printf "%12.3f ] µs\n" time_hahn_echo[system.n_time_step]*1e6
        else
            @printf "Times: [ %12.3f, %12.3f, ..., %12.3f ] µs\n" time_hahn_echo[1]*1e6 time_hahn_echo[2]*1e6 time_hahn_echo[system.n_time_step]*1e6
        end
    else
        @printf "Setting custom delay times (tau):\n"
        time_hahn_echo = system.times
        n_time_step = size(time_hahn_echo)[1]
        @printf "Times: ["
        for idx in 1:n_time_step-1
            @printf " %12.3f, " time_hahn_echo[idx]*1e6
        end
        @printf "%12.3f ] µs\n" time_hahn_echo[n_time_step]*1e6
    end

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
        if err
            error("Unsuitable settings found (see above)")
        end

        intensity = cce_hf_analytic(distance_coordinates_el_nucs,n_nuc,
                                    gamma_n,gamma_electron_sc,time_hahn_echo,system)
        # dummies for iCCE1 and iCCE2
        iCCE0 = ones(size(time_hahn_echo))
        iCCE1 = ones(size(time_hahn_echo))
        iCCE2 = intensity
    elseif system.simulation_type == "exact"
        if system.report_pair_contrib
            error("pair contribution report only implemented for analytic model")
        end
        intensity,iCCE0,iCCE1,iCCE2 = cce_exact(distance_coordinates_el_nucs,n_nuc,
                gamma_n,gamma_electron,time_hahn_echo,system)
    else
        print("Unkonwn simulation type: ",system.simulation_type,"\n")
        intensity,iCCE0,iCCE1,iCCE2 = zeros(size(time_hahn_echo))
    end

    # return intensity and time
    return time_hahn_echo,intensity,iCCE0,iCCE1,iCCE2
end

"""
    cce_hf_analytic(distance_coordinates_el_nucs,n_nuc,r_max_bath,gamma_n,gamma_electron,B0,time_hahn_echo)

Use the analytic pair product approximation (high-field approximation)

should be called by cce (rather than directly); distance_coordinates_el_nucs are the distances 
of the bath nuclei from the spin center, gamma_n and gamma_electron the respective gyromagnetic 
rations required for computing the dipolar interaction, B0 is the magnetic field direction,
time_hahn_echo the time set for which the signal is to be computed

the function returns the intensity for the given time set on time_hahn_echo
note that time_hahn_echo corresponds to the HE delay time, so the full time from the initial
pi/2 pulse until echo is 2 * time_hahn_echo (2 tau)

References used:

Witzel, W. M.; Das Sarma, S. Quantum theory for electron spin decoherence in-
duced by nuclear spin dynamics in semiconductor quantum computer architectures:
Spectral diffusion of localized electron spins in the nuclear solid-state environment.
Phys. Rev. B 2006, 74, 035322 https://doi.org/10.1103/physrevb.74.035322

Jeschke, G. Nuclear Pair Electron Spin Echo Envelope Modulation. J. Magn. Reson. Open 2023, 14, 100094. 
https://doi.org/10.1016/j.jmro.2023.100094. Erratum: J. Magn. Reson. Open 2023, 14, 100094. 
https://doi.org/10.1016/j.jmro.2023.100115.

"""
function cce_hf_analytic(distance_coordinates_el_nucs,n_nuc,
                         gamma_n,gamma_electron,time_hahn_echo,system)
    
    r_max_bath = system.r_max_bath*aacm
    B0 = system.B0

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

    # in case, there are no nuclear spin baths within the given r_max, do nothing --> return intensity of 1 for all time steps
    if !isempty(A_n)
        minA = minimum(abs.(A_n))
        maxA = maximum(abs.(A_n))
        print("\n")
        @printf "Hyperfine coupling constants (min, max) in Hz: %12.4g %12.4g\n\n" minA maxA
    else 
        nothing
    end

    # set up list of pairs
    n_pairs, pair_list, n_pair_contr = make_pair_list(distance_coordinates_el_nucs,r_max_bath,
                                                      system.select,system.rselect)

    @printf "Total number of pairs in bath: %12i\n" n_nuc*(n_nuc-1)÷2
    @printf "Screened number of pairs:      %12i\n" n_pairs
    @printf "Radius for spin bath was:      %12.2f Å\n" system.r_max
    @printf "Pair screening distance was:   %12.2f Å\n\n" system.r_max_bath

    flush(stdout)

    if system.report_pair_contrib
        print("Writing out pair contributions to: ",system.pair_log_file,"\n\n")
        flush(stdout)
        plf = open(system.pair_log_file,"w")
    end

    # initialize Intensity to 1. for all times
    intensity = ones(n_time_step)
    
    # loop over pairs of bath nuclei 
    for ipair = 1:n_pairs
            n = pair_list[ipair][1]
            m = pair_list[ipair][2]
            
            vr_nm = (distance_coordinates_el_nucs[m] - distance_coordinates_el_nucs[n]) 
            r_nm = norm(vr_nm)

            # compute angle with magnetic field
            r_nm_x_B0 = cross(vr_nm, B0)
            r_nm_dot_B0 = dot(vr_nm, B0)
            theta_nm = atan(norm(r_nm_x_B0), r_nm_dot_B0)
            # nuclear dipole interaction
	        b_nm = - gamma_n^2 * hbar * (1 - 3 * cos(theta_nm)^2) / r_nm^3
            # modulation depth
	        c_nm = (A_n[n] - A_n[m]) / b_nm
            L_nm = 4.0*((c_nm^2) / (1 + c_nm^2)^2)
            # nuclear zero-quantum frequency
            #or: w_nm = 0.5 * b_nm * sqrt(1 + c_nm^2)
            w_nm = 0.5 * sqrt( b_nm^2  + (A_n[n] - A_n[m])^2 )

            if system.report_pair_contrib
                r_cent = norm(0.5*(distance_coordinates_el_nucs[m] + distance_coordinates_el_nucs[n]))
                @printf(plf,"%18.6f %18.6f %20.10g %20.10g",r_cent/aacm,r_nm/aacm,L_nm,abs(w_nm))
            end
            
            # for testing: use exponential (not recommended)
	        if system.use_exp
                for j in 1:size(time_hahn_echo)[1]
                    # note that time_hahn_echo contains tau (delay)
                    v_nm = - L_nm * (sin(0.5 * w_nm * time_hahn_echo[j]))^4
                    intensity[j] = intensity[j] * exp(v_nm)
                end
            else
                for j in 1:size(time_hahn_echo)[1]
                    # note that time_hahn_echo contains tau (delay)
                    v_nm = - L_nm * (sin(0.5 * w_nm * time_hahn_echo[j]))^4
                    intensity[j] = intensity[j] * (1 + v_nm)
                    if system.report_pair_contrib
                        @printf(plf," %20.10g",(1 + v_nm))
                    end
                end
            end

            if system.report_pair_contrib
                write(plf,"\n")
            end 
        
    end    # loop over pair_list    
     
    if system.report_pair_contrib
        close(plf)
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
    
    A = fact*(I(3) - 3.0 * rn*transpose(rn))

    # assuming this form of the HFC:  S^+ A I
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
    get_signal(mx,my)

wrapper function for consistent definition of detected signal
"""
function get_signal(mx,my)
     return sqrt(mx*mx + my*my)
     #return abs(mx)
end


"""
    e_contribution(various args)

core routine for "exact" propagation of free spin system (internal use only)
"""
function e_contribution(dim,Hmag,Mmat,PI,rho0,time_he)

    nt = size(time_he)[1]

    # only magnetic field contribution
    Hmat = Hmag

    mx = real(tr(Mmat[1]*rho0))
    my = real(tr(Mmat[2]*rho0))

    # initial signal0
    signal_0 = get_signal(mx,my)

    # propagate using the eigenbasis:
    Hev, UH = eigen(Hmat)
    # initialise kernel matrix
    K = zeros(ComplexF64,dim,dim)

    signal = zeros(nt)
    for idx in 1:nt

        # set up kernel
        for j in 1:dim
            K[j,j] = exp(-1im*time_he[idx]*Hev[j])
        end
        # get time propagator
        Ut = UH * K * adjoint(UH)

        U = Ut * PI * Ut

        rho_tau = U * rho0 * adjoint(U)

        mx = real(tr(Mmat[1]*rho_tau))
        my = real(tr(Mmat[2]*rho_tau))

        signal[idx] = get_signal(mx,my)/signal_0

    end

    return signal

end


"""
    n_e_contribution(various args)

core routine for "exact" propagation of 2-particle system (internal use only)
"""
function n_e_contribution(dim_el,dim_nuc,A_1,Hmag,Mmat,Smat,Imat,PI,rho0,time_he)

    dim = dim_nuc*dim_el
    nt = size(time_he)[1]

    # spin H for 2 spin system
    I1 = unit_mat(dim_nuc)

    # interaction contribution
    Hmat = zeros(dim,dim)
    for i in 1:3
        for j in 1:3
            Hmat += kron(Smat[i],Imat[j]).*A_1[i,j]
        end
    end

    # add magnetic contribution
    Hmat += Hmag

    # Expand π pulse operator to tensor basis
    PIexp = kron(PI,I1)

    mx = real(tr(Mmat[1]*rho0))
    my = real(tr(Mmat[2]*rho0))

    # initial signal0
    signal_0 = get_signal(mx,my)

    # propagate using the eigenbasis:
    Hev, UH = eigen(Hmat)
    # initialise kernel matrix
    K = zeros(ComplexF64,dim,dim)

    signal = zeros(nt)
    for idx in 1:nt

        # set up kernel
        for j in 1:dim
            K[j,j] = exp(-1im*time_he[idx]*Hev[j])
        end
        # get time propagator
        Ut = UH * K * adjoint(UH)

        U = Ut * PIexp * Ut

        rho_tau = U * rho0 * adjoint(U)

        mx = real(tr(Mmat[1]*rho_tau))
        my = real(tr(Mmat[2]*rho_tau))

        signal[idx] = get_signal(mx,my)/signal_0

    end

    return signal

end



"""
    n_n_e_contribution(various args)

core routine for "exact" propagation of 3-particle system (internal use only)
"""
function n_n_e_contribution(dim_el,dim_nuc,gamma_n,r12,A_1,A_2,Hmag,Mmat,Smat,Imat,PI,rho0,time_he)

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

    # interaction contribution
    Hmat = zeros(dim,dim)
    for i in 1:3
        for j in 1:3
            Hmat += kron(Smat[i],kron(Imat[j],I1)).*A_1[i,j]
            Hmat += kron(Smat[i],kron(I1,Imat[j])).*A_2[i,j]
            Hmat += kron(S1,kron(Imat[i],Imat[j])).*B[i,j]
        end
    end

    # add magnetic contribution
    Hmat += Hmag

    # expand π pulse operator to full tensor basis
    PIexp = kron(PI,kron(I1,I1))

    mx = real(tr(Mmat[1]*rho0))
    my = real(tr(Mmat[2]*rho0))

    # initial signal
    signal_0 = get_signal(mx,my)

    # propagate using the eigenbasis:
    Hev, UH = eigen(Hmat)
    # initialise kernel matrix
    K = zeros(ComplexF64,dim,dim)

    signal = zeros(nt)
    for idx in 1:nt

        # set up kernel
        for j in 1:dim
            K[j,j] = exp(-1im*time_he[idx]*Hev[j])
        end
        # get time propagator
        Ut = UH * K * adjoint(UH)

        U = Ut * PIexp * Ut

        rho_tau = U * rho0 * adjoint(U)

        mx = real(tr(Mmat[1]*rho_tau))
        my = real(tr(Mmat[2]*rho_tau))

        signal[idx] = get_signal(mx,my)/signal_0
        
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
    yrot = [ costh 0. sinth ;  0.  1.  0.; -sinth  0.  costh]

    # rotation around z by -phi
    zrot = [ cosphi  -sinphi  0. ; sinphi  cosphi  0.; 0.  0.  1. ]

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
function make_pair_list(distance_coordinates_el_nucs,r_max_bath,select,rselect)

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
            if norm(r12) < 0.5*aacm
                @printf "WARNING: Pair %4i %4i has r = %10.4f AA \n" i j norm(r12)/aacm
            end
            if norm(r12) < 0.1*aacm
                println("This must be wrong!")
                error("Nonsensical spin bath coordinates!")
            end
            if select
                if norm(r12) < rselect[3]*aacm || norm(r12) > rselect[4]*aacm
                    continue
                end
                ravg = norm(0.5*(distance_coordinates_el_nucs[i] + distance_coordinates_el_nucs[j]))
                if ravg < rselect[1]*aacm || ravg > rselect[2]*aacm
                    continue
                end

            end

            n_pairs += 1

            push!(pair_list,[i,j])

            n_pair_contr[i] += 1
            n_pair_contr[j] += 1

        end
    end

    return n_pairs,pair_list,n_pair_contr

end


function cce_exact(distance_coordinates_el_nucs,n_nuc,
    gamma_n,gamma_el,time_hahn_echo,system)

#    if ! is_unit(system.magnetic_axes,1e-8)
#        print("non unit magnetic axes not yet debugged!")
#        error("non unit magnetic axes not yet debugged!")
#    end

    r_max_bath = system.r_max_bath*aacm
    B00 = system.B0

    s_nuc = system.s_nuc
    s_el = system.s_el

    mag_axes = system.magnetic_axes

    n_time_step = size(time_hahn_echo)[1]

    intensity = ones(n_time_step)
    intensity_CCE0 = ones(n_time_step)
    intensity_CCE1 = ones(n_time_step)
    intensity_CCE2 = ones(n_time_step)
    intensity_correction = ones(n_time_step)

    # Here, it is more convenient to keep the magnetic field at (0,0,1) and to
    # transform the nuclear coordinates.
    # This is because the treatment employs the rotating wave approximation and
    # uses the magnetic field as the main quantization axis z
    
    # field strength
    Bstrength = norm(B00)

    Rmat = get_rot_for_unit_vector(B00/Bstrength)

    @printf "Found this field: [ %12.6g %12.6g %12.6g ] G\n" B00[1] B00[2] B00[3]
    @printf "Norm: %12.6g G\n\n" Bstrength

    print("Transforming by:\n")
    display(Rmat)
    print("\n\n")

    println("B field transformed:")
    display(Rmat' * B00)
    print("\n")

    n_pairs, pair_list, n_pair_contr = make_pair_list(distance_coordinates_el_nucs,r_max_bath,
                                                      system.select,system.rselect)

    @printf "Total number of pairs in bath: %12i\n" n_nuc*(n_nuc-1)÷2
    @printf "Screened number of CCE2 pairs: %12i\n" n_pairs
    @printf "Radius for spin bath was:      %12.2f Å\n" system.r_max
    @printf "Pair screening distance was:   %12.2f Å\n\n" system.r_max_bath

    B0 = [0.,0.,Bstrength]

    # transform coordinates
    distance_el_nuc_traf = []
    for coord in distance_coordinates_el_nucs
        coord_traf = Rmat' * coord
        push!(distance_el_nuc_traf,coord_traf)
    end

    # "g tensor" from gamma and magnetic axes
    gamma_t = zeros(3,3)
    gamma_t[1,1] = gamma_el[1]; gamma_t[2,2] = gamma_el[2]; gamma_t[3,3] = gamma_el[3] 
    mag_axes_t = Rmat' * mag_axes
    gamma_t = mag_axes_t * gamma_t * mag_axes_t'

    println("mag_axes:")
    display(mag_axes)
    print("\n")

    println("mag_axes_t:")
    display(mag_axes_t)
    print("\n")


    print("gamma_t: \n")
    display(gamma_t)
    print("\n")

    print("gamma_n: ",gamma_n,"\n")
    
    print("\n")
    print("computing hyperfine couplings ...\n")

    # compute list of all hyperfine interactions
    A_list = []
    for n in 1:n_nuc
        r1 = distance_el_nuc_traf[n]
        A_current = hyperfine(gamma_t,gamma_n,r1)
        push!(A_list,A_current)
    end

    dim_el = trunc(Int,(2*s_el+1))
    dim_nuc = trunc(Int,(2*s_nuc+1))

    # get all spin matrices
    Smat = []
    push!(Smat,sx_mat(s_el))
    push!(Smat,sy_mat(s_el))
    push!(Smat,sz_mat(s_el))
    # and a unit matrix of the same dimension
    S1 = unit_mat(dim_el)

    # get all spin matrices for nuclei
    Imat = []
    push!(Imat,sx_mat(s_nuc))
    push!(Imat,sy_mat(s_nuc))
    push!(Imat,sz_mat(s_nuc))
    # and again a unit matrix of the same dimension
    I1 = unit_mat(dim_nuc)

    # initial spin state |+><+| (max. comp along X)
    PI = zeros(dim_el,dim_el)
    # note: currently, we assume only interaction with electronic spin
    # this is likely a good approximation, althoung in practice there might be
    # cross talk to the nuclei
    if s_el == 0.5
        # represent π/2 and π pulses (in lab system!)

        # eigenvalues and vectors of Pauli matrix
        #
        evSX,vcSX = eigen(Smat[1]*2.0)
        
        # π pulse:
        KPI = zeros(ComplexF64,2,2)
        for i = 1:2
            KPI[i,i] = exp(-0.5im*evSX[i]*π)
        end
        PI = vcSX * KPI * adjoint(vcSX)

        # init electronic system to situation after ideal π/2 pulse:
        rho0_el = (vcSX[:,2]) * adjoint(vcSX[:,2])
    else
        error("S>0.5 not yet implemented")
        # TODO: diagonalize Sx and get eigenvector for -M state == |+> ; rho = |+><+|
    end
        
    # initial nuclear state: all M equal likely (kT >> delta E)
    rho0_nuc = I1.*(1. / trunc(Int,2*s_nuc+1) )

    # do lower order CCE (0 + 1); recommended, wrong results unless very high field
    if system.do_cce1
        # free evolution contribution
        # (relevant for anisotropic systems)
        Mmat0 = []
        for i = 1:3
            mat = zeros(dim_el,dim_el)
            for j = 1:3
                mat -= Smat[j].*gamma_t[i,j]
            end
            push!(Mmat0,mat)
        end
        Hmag0 = zeros(dim_el,dim_el)
        B0t = gamma_t' * B0
        for i = 1:3
            Hmag0 += Smat[i].*B0t[i]
        end

        intensity_CCE0 = e_contribution(dim_el,Hmag0,Mmat0,PI,rho0_el,time_hahn_echo)

        # intensity correction for CCE2
        intensity_correction = intensity_CCE0.^n_pairs

    end
    
    if system.do_cce1
        # CCE1 contrubutions
        print("\n")
        print("Computing CCE1 contributions...\n")

        dimH1 = dim_el * dim_nuc

        # set magnetic interaction matrix
        # gamma_t is by default positive for electrons, so we need a minus sign
        # for the nuclei, we assume that gamma_n carries the appropriate sign
        Mmat1 = []
        for i = 1:3
            mat = zeros(dimH1,dimH1)
            for j = 1:3
                mat -= kron(Smat[j],I1).*gamma_t[i,j]
            end
            mat += kron(S1,Imat[i]).*gamma_n
            push!(Mmat1,mat)
        end

        # magnetic field contribution is universal
        # signs: H = - M B ... therefore signs are reversed relative to Mmat1 setup
        Hmag1 = zeros(dimH1,dimH1)
        B0t = gamma_t' * B0 
        for i = 1:3
            Hmag1 += kron(Smat[i],I1).*B0t[i]
            Hmag1 -= kron(S1,Imat[i]).*(B0[i]*gamma_n)
        end

        # set initial density matrix for 2-particle system    
        rho0 = kron(rho0_el,rho0_nuc)*(1. + 0*im)

        for n in 1:n_nuc
            A_n = A_list[n]
            ne_cont = n_e_contribution(dim_el,dim_nuc,A_n,Hmag1,Mmat1,Smat,Imat,PI,rho0,time_hahn_echo)

            ne_cont ./= intensity_CCE0  # correct for CCE0 (free propagation)

            intensity_CCE1 = intensity_CCE1 .* ne_cont

            # number of pairs to which this CCE1 increment will contribute
            mult = n_pair_contr[n]
            # update correction for CCE2
            intensity_correction = intensity_correction .* ne_cont.^mult

        end

    end

    print("\n")
    print("Computing CCE2 contributions ...\n")

    dimH = dim_el * dim_nuc^2

    # Magnetic interaction matrix
    # gamma_t is by default positive for electrons, so we need a minus sign
    # for the nuclei, we assume that gamma_n carries the appropriate sign
    Mmat = []
    for i = 1:3
        mat = zeros(dimH,dimH)
        for j = 1:3
            mat -= kron(Smat[j],kron(I1,I1)).*gamma_t[i,j]
        end
        mat += kron(S1,kron(Imat[i],I1)).*gamma_n
        mat += kron(S1,kron(I1,Imat[i])).*gamma_n
        push!(Mmat,mat)
    end

    # magnetic field contribution is universal
    # signs: H = - M B ... therefore signs are reversed relative to Mmat1 setup
    Hmag = zeros(dimH,dimH)
    B0t = gamma_t' * B0
    for i = 1:3
        Hmag += kron(Smat[i],kron(I1,I1)).*B0t[i]
        Hmag -= kron(S1,kron(Imat[i],I1)).*(B0[i]*gamma_n)
        Hmag -= kron(S1,kron(I1,Imat[i])).*(B0[i]*gamma_n)
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
            nne_cont = n_n_e_contribution(dim_el,dim_nuc,gamma_n,r1-r2,A_n,A_m,Hmag,Mmat,Smat,Imat,PI,rho0,time_hahn_echo)

            intensity_CCE2 = intensity_CCE2 .* nne_cont
        end

    else
        
        # get number of threads and prepare result array for each thread
        nthreads = Threads.nthreads()
        int_thread = [ones(n_time_step) for _ in 1:nthreads]

        @printf "... Entered threaded route (nthreads = %i)\n" nthreads

        # unset BLAS multithreading
        n_BLAS_threads = BLAS.get_num_threads()
        BLAS.set_num_threads(1)

        # parallel loop
        Threads.@threads for ipair in eachindex(pair_list)

            thread_id = Threads.threadid()

            n = pair_list[ipair][1]
            m = pair_list[ipair][2]

            r1 = distance_el_nuc_traf[n]
            r2 = distance_el_nuc_traf[m]

            A_n = A_list[n]
            A_m = A_list[m]
            nne_cont = n_n_e_contribution(dim_el,dim_nuc,gamma_n,r1-r2,A_n,A_m,Hmag,Mmat,Smat,Imat,PI,rho0,time_hahn_echo)

            int_thread[thread_id] .*= nne_cont
        end
        # multiply contributions from threads:
        intensity_CCE2 = reduce(.*, int_thread)

        # reset BLAS
        BLAS.set_num_threads(n_BLAS_threads)

    end

    intensity_CCE2 = intensity_CCE2 ./ intensity_correction

    intensity = intensity_CCE0 .* intensity_CCE1 .* intensity_CCE2

    return intensity, intensity_CCE0, intensity_CCE1, intensity_CCE2

end

end
