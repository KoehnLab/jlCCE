module RotationMatrices

export rotation_matrices

push!(LOAD_PATH,"../../src")

using LinearAlgebra

function rotation_matrices(theta::Float64,phi::Float64) 
     # rotation around y axis
     rot_mat_y_theta = [cos(theta) 0 sin(theta); 0 1 0; -sin(theta) 0 cos(theta)]
     # rotation around z axis
     rot_mat_z_phi = [cos(phi) -sin(phi) 0; sin(phi) cos(phi) 0; 0 0 1]

     return rot_mat_y_theta,rot_mat_z_phi
end

end
