push!(LOAD_PATH,"../src")

using jlCCE
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
            get_coordinates("vodbm2.cif",23,1)
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
