using jlCCE

# set system
voacac2 = SpinSystem("voacac2.cif","V",1)

voacac2.g_factor = [2.3,2.3,2.3]  # isotropic here
voacac2.magnetic_axes = I(3)
voacac2.nuc_spin_bath = "H"
voacac2.gn_spin_bath = 5.58569468
voacac2.B0 = [0.,0.,1.]
voacac2.rmin = 0.
voacac2.rmax = 30. # AA
voacac2.t_min = 0.
voacac2.t_max = 1.5e-5 # s
voacac2.n_time_step = 30

# run
time, intesity = cce(voacac2)
