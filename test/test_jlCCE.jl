push!(LOAD_PATH,"../src")

using jlCCE
using jlCCEtools
using SpinBase
using readCIF
using Test, LinearAlgebra


@testset "Spin functions" begin

    @test sz_mat(0.5) ≈ diagm([0.5,-0.5])

    @test sz_mat(1.0) ≈ diagm([1.,0.,-1.])

    @test sz_mat(1.5) ≈ diagm([1.5,0.5,-0.5,-1.5])

    @test sx_mat(0.5) ≈ [0. 0.5; 0.5 0.]

    r2h = sqrt(2)*0.5
    @test sx_mat(1.0) ≈ [0. r2h 0.; r2h 0. r2h; 0. r2h 0.]

    r3h = sqrt(3)*0.5
    @test sx_mat(1.5) ≈ [0. r3h 0. 0.; r3h 0. 1. 0.; 0. 1. 0. r3h; 0. 0. r3h 0.]

    @test sy_mat(0.5) ≈ [0. -0.5im; 0.5im 0.]

    @test sy_mat(1.0) ≈ [0. -r2h 0.; r2h 0. -r2h; 0. r2h 0.] * 1im

    @test sy_mat(1.5) ≈ [0. -r3h 0. 0.; r3h 0. -1. 0.; 0. 1. 0. -r3h; 0. 0. r3h 0.] * 1im

    @test sp_mat(1.5) ≈ [0. r3h 0. 0.; 0. 0. 1. 0.; 0. 0. 0. r3h; 0. 0. 0. 0.] * 2

    @test sm_mat(1.5) ≈ [0. 0. 0. 0.; r3h 0. 0. 0.; 0. 1. 0. 0.; 0. 0. r3h 0.] * 2

end


@testset "Read CIF" begin

    lattice,coord_e_spin,coords_n_spins_unit_cell = 
            get_coordinates("vodbm2.cif",23,1,1)
    @test lattice[1] ≈ [16.462, 0., 0.] atol=1e-3
    @test lattice[2] ≈ [0.,20.269,0.] atol=1e-3
    @test lattice[3] ≈ [0.,0.,14.95] atol=1e-3

    @test coord_e_spin ≈ [14.11814044,7.55729665,0.797134] atol=1e-6

    @test coords_n_spins_unit_cell[3] ≈ [9.901893, 2.3653923, 13.468455] atol=1e-6

    cell_list = set_supercell_list(20.,lattice)

    distance_coordinates_el_nucs,n_nuc = 
         get_bath_list(0.,20.,lattice,coords_n_spins_unit_cell,coord_e_spin)

    @test n_nuc == 1196

    @test distance_coordinates_el_nucs[3] ≈ [-14.49512024,-9.20922015,-3.347604] atol=1e-6
    @test distance_coordinates_el_nucs[1000] ≈ [ 5.57699636,11.10437165,10.089456] atol=1e-6

end

@testset "jlCCE integration test 1" begin

    spinsystem = SpinSystem("test","Pd",1)  # test mode, nucleus is dummy
    spinsystem.use_exp = false    
    spinsystem.g_factor = [1.9846,1.9846,1.9846]
    spinsystem.B0 = [0.,0.,100.]*10000 
    spinsystem.t_min = 0.
    spinsystem.t_max = 1.5e-4
    spinsystem.n_time_step = 101    

    times,intensity = cce(spinsystem)

    @test intensity[1:9]   ≈ [  1.0,     0.9999999725400402, 0.9999995607920832,
                              0.9999977777871571, 0.9999929823578974, 0.9999828847984324,
                              0.9999645547652702, 0.999934431396304,  0.9998883356186314] atol=1e-9

    @test intensity[62:70] ≈ [ 0.7543971635543107, 0.7417395059784252, 0.7287835484535932,
                               0.7155399813065866, 0.7020203427723999, 0.6882370044006308,
                               0.6742031541145092, 0.6599327769574507, 0.645440633568423] atol=1e-9

    @test times[40] ≈ 5.85e-5 atol=1e-15

    spinsystem.simulation_type="exact"

    times,intensityX = cce(spinsystem)

    @test intensity ≈ intensityX atol=1e-6

end

@testset "jlCCE integration test 2" begin

    spinsystem = SpinSystem("test","Pd",1)  # test mode, nucleus is dummy
    spinsystem.use_exp = false    
    spinsystem.g_factor = [1.9846,1.9846,1.9846]
    spinsystem.B0 = [100.,10.,100.]*10000    # tilted B field
    spinsystem.t_min = 0.
    spinsystem.t_max = 1.5e-4
    spinsystem.n_time_step = 101    

    times,intensityA = cce(spinsystem)

    spinsystem.simulation_type="exact"

    times,intensityX = cce(spinsystem)

    @test intensityA ≈ intensityX atol=1e-6

end

@testset "jlCCEtools tests" begin

    times  = [0.00000000E+00,3.06122449E-07,6.12244898E-07,9.18367347E-07,1.22448980E-06,
              1.53061224E-06,1.83673469E-06,2.14285714E-06,2.44897959E-06,2.75510204E-06,
              3.06122449E-06,3.36734694E-06,3.67346939E-06,3.97959184E-06,4.28571429E-06,
	      4.59183673E-06,4.89795918E-06,5.20408163E-06,5.51020408E-06,5.81632653E-06]
    signal = [1.00000000E+00,9.98974170E-01,9.90650813E-01,9.74138471E-01,9.50294218E-01,
              9.16723414E-01,8.66878057E-01,7.99398290E-01,7.20677728E-01,6.33920983E-01,
              5.41495348E-01,4.51089619E-01,3.69808807E-01,2.98071630E-01,2.35133270E-01,
              1.81208695E-01,1.36952725E-01,1.01017267E-01,7.25367807E-02,5.10004848E-02 ]

    T2 = get_decay_time(times,signal)

    @test T2 ≈ 3.68170253E-6 atol=1e-11

end
