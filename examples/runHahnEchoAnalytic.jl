push!(LOAD_PATH,"../src")

using jlCCE
using BenchmarkTools

# set system - use the simple constructor
spinsystem = SpinSystem("vmnt3.cif","V",1)  # <--- note: put this file into the folder

# modify the values (SpinSystem creates a mutable object):
spinsystem.g_factor = [2.25,2.25,2.25]  # isotropic here
# pddbm2: g = 2.25
# vmnt3: 1.9846,1.9846,1.9846



# defaults, no need to define here
#voacac2.magnetic_axes = [1.0 0.0 0.0 ; 0.0 1.0 0.0 ; 0.0 0.0 1.0]
#voacac2.nuc_spin_bath = "H"
#voacac2.gn_spin_bath = 5.58569468
spinsystem.B0 = [0.,0.,1.]
spinsystem.r_min = 0.
spinsystem.r_max = 10. # AA
spinsystem.t_min = 0.
spinsystem.t_max = 1.5e-5 # s
spinsystem.n_time_step = 20


# run
@time time_hahn_echo,intensity = cce(spinsystem)
print("\n")
print("Our great result:\n")
print("time of the Hahn echo: ",time_hahn_echo,"\n")
print("Simulated intensity: ",intensity,"\n")
print("\n")

