push!(LOAD_PATH,"../src")

using jlCCE

# set system - use the simple constructor
vmnt3 = SpinSystem("vmnt3.cif","V",1)  # <--- note: put this file into the folder

# modify the values (SpinSystem creates a mutable object):
vmnt3.g_factor = [1.9846,1.9846,1.9846]  # isotropic here
# defaults, no need to define here
#voacac2.magnetic_axes = [1.0 0.0 0.0 ; 0.0 1.0 0.0 ; 0.0 0.0 1.0]
#voacac2.nuc_spin_bath = "H"
#voacac2.gn_spin_bath = 5.58569468
vmnt3.B0 = [0.,0.,1.]
vmnt3.r_min = 0.
vmnt3.r_max = 15. # AA
vmnt3.t_min = 0.
vmnt3.t_max = 1.5e-5 # s
vmnt3.n_time_step = 20

# run
time, intensity = cce(vmnt3)

print("Our great result:\n")
print(intensity)
print("\n")

