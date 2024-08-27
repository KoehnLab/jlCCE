push!(LOAD_PATH,"../src")

using jlCCE

# set system - use the simple constructor
voacac2 = SpinSystem("vmnt3.cif","V",1)  # <--- note: put this file into the folder

# modify the values (SpinSystem creates a mutable object):
voacac2.g_factor = [2.3,2.3,2.3]  # isotropic here
# defaults, no need to define here
#voacac2.magnetic_axes = [1.0 0.0 0.0 ; 0.0 1.0 0.0 ; 0.0 0.0 1.0]
#voacac2.nuc_spin_bath = "H"
#voacac2.gn_spin_bath = 5.58569468
voacac2.B0 = [0.,0.,1.]
voacac2.r_min = 0.
voacac2.r_max = 30. # AA
voacac2.t_min = 0.
voacac2.t_max = 1.5e-5 # s
voacac2.n_time_step = 30

# run
time, intensity = cce(voacac2)

print("Our great result:\n")
print(intensity)
print("\n")

