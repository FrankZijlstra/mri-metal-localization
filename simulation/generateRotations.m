function [xyz] = generateRotations (subdivisions)
% Generates a near perfectly uniform distribution of vectors on the unit
% sphere by subdividing the triangles of an icosahedron.

phi = (1+sqrt(5))/2;  % Golden ratio

xyz = [0   1  phi 
           0  -1  phi 
           0   1 -phi 
           0  -1 -phi 
           1  phi  0  
          -1  phi  0  
           1 -phi  0  
          -1 -phi  0  
           phi 0   1  
          -phi 0   1  
           phi 0  -1  
          -phi 0  -1];

F = [1  2  9
     1  9  5
     1  5  6
     1  6  10
     1  10 2
     2  7  9
     9  7  11
     9  11 5
     5  11 3
     5  3  6
     6  3  12
     6  12 10
     10 12 8
     10 8  2
     2  8  7
     4  7  8
     4  8  12
     4  12 3
     4  3  11
     4  11 7];

Z = xyz(1,:)';   % Z passes through vertex 1
X = xyz(2,:)';   % Choose adjacent vertex as an approximate X
Y = cross(Z,X);  % Y is perpendicular to Z and this approx X
X = cross(Y,Z);  % Final X is perpendicular to Y and Z
X = X/norm(X); Y = Y/norm(Y); Z = Z/norm(Z);  % Ensure unit vectors
xyz = ([X Y Z]'*xyz')';  % Transform points;


normalize = @(x) x ./ norm(x);

for J=1:size(xyz,1)
    xyz(J,:) = normalize(xyz(J,:));
end

% Subdivide the mesh
for D=1:subdivisions
    newF = [];
    K = 1;
    
    for J=1:size(F,1)
        newP1 = normalize(xyz(F(J,1),:)+xyz(F(J,2),:));
        newP2 = normalize(xyz(F(J,2),:)+xyz(F(J,3),:));
        newP3 = normalize(xyz(F(J,1),:)+xyz(F(J,3),:));
        
        ind1 = size(xyz,1)+1;
        xyz(ind1,:) = newP1;
        ind2 = size(xyz,1)+1;
        xyz(ind2,:) = newP2;
        ind3 = size(xyz,1)+1;
        xyz(ind3,:) = newP3;
        

        newF(K,:) = [ind1 ind2 ind3];
        K = K+1;
        newF(K,:) = [F(J,1) ind1 ind3];
        K = K+1;
        newF(K,:) = [ind1 F(J,2) ind2];
        K = K+1;
        newF(K,:) = [ind2 F(J,3) ind3];
        K = K+1;
    end
    
    F = newF;
end

% Remove double vertices (each edge generates a midpoint twice)
xyz = unique(xyz,'rows');

% Determine mirrored vertices
N = size(xyz,1);
res = zeros(N,N);
for I=1:N
    for J=1:N
        if (xyz(I,1) == -xyz(J,1) && xyz(I,2) == -xyz(J,2) && xyz(I,3) == -xyz(J,3))
            res(I,J) = 1;
        end
    end
end

% Keep only one of each mirrored vertex
mask = zeros(1,N);
xyzNew = [];
for I=1:N
    if (~mask(I))
        xyzNew(size(xyzNew,1)+1,:) = xyz(I,:);
        mask = mask | res(I,:);
    end
end

xyz = xyzNew;

end

