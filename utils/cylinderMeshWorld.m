function [msh] = cylinderMeshWorld(R, L, location, orientationVector)

N = 20;

theta = linspace(0,2*pi,N);

perpV1 = cross(orientationVector,[1 0 0]);
if (all(perpV1 == 0))
    perpV1 = cross(orientationVector,[0 1 0]);
end
perpV1 = perpV1 / norm(perpV1);
perpV2 = cross(perpV1, orientationVector);

x = R * [zeros(1,N); cos(theta); cos(theta); zeros(1,N)];
y = R * [zeros(1,N); sin(theta); sin(theta); zeros(1,N)];
z = L * repmat([0; 0; 1; 1],1,N) - L/2;

coords = repmat(permute(location,[1 3 2]), size(x)) + repmat(x,[1 1 3]).*repmat(permute(perpV1,[1 3 2]),[size(x) 1]) + repmat(y,[1 1 3]).*repmat(permute(perpV2,[1 3 2]),[size(x) 1]) + repmat(z,[1 1 3]).*repmat(permute(orientationVector,[1 3 2]),[size(x) 1]);
msh = surf2patch(coords(:,:,1),coords(:,:,2),coords(:,:,3),'triangles');

end
