module RotationMatrices

export rotation_matrices

push!(LOAD_PATH,"../../src")

using LinearAlgebra

function rotation_matrices(theta::Float64,phi::Float64) 
     print("\n")
     println("Start Rotation of the magnetic field")
     println("==================================== \n")
     #println("given B0: ", B0)
     #println("Angle around the y axis: ", rad2deg(theta))
     #println("Angle around the z axis: ", rad2deg(phi))

     # define rotation matrixes -- cartesian coordinate system 
     # rotation around y axis
     rot_mat_y_theta = [cos(theta) 0 sin(theta); 0 1 0; -sin(theta) 0 cos(theta)]
     # rotation around z axis
     rot_mat_z_phi = [cos(phi) -sin(phi) 0; sin(phi) cos(phi) 0; 0 0 1]

     # rotate B0
#     B0_rot = transpose(R_m) *(rot_mat_z * (rot_mat_y * B0))
#     println("Rotated magnetic field: ", B0_rot, "\n")

     return rot_mat_y_theta,rot_mat_z_phi
end


end
