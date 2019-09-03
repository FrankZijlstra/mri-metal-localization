function [] = orthoplotWorldMatrix (im, sliceY, sliceX, sliceZ, worldMatrix, showXYZ)
held = ishold;

if (nargin < 6 || isempty(showXYZ))
    showXYZ = [1 1 1];
end

cmap = gray(256);
hold on

if (showXYZ(1))
    C1 = worldMatrix * [sliceX-1; -0.5; -0.5; 1];
    C3 = worldMatrix * [sliceX-1; size(im,1)-0.5; -0.5; 1];
    C4 = worldMatrix * [sliceX-1; size(im,1)-0.5; size(im,3)-0.5; 1];
    C2 = worldMatrix * [sliceX-1; -0.5; size(im,3)-0.5; 1];

    X = [C1(1) C2(1); C3(1) C4(1)];
    Y = [C1(2) C2(2); C3(2) C4(2)];
    Z = [C1(3) C2(3); C3(3) C4(3)];

    surf(X,Y,Z,'FaceColor','texturemap','Cdata',squeeze(im(:,sliceX,:,1)) * size(cmap,1),'CdataMapping', 'direct', 'EdgeColor', 'none')
end

if (showXYZ(2))
    C1 = worldMatrix * [-0.5; sliceY-1; -0.5; 1];
    C3 = worldMatrix * [size(im,2)-0.5; sliceY-1; -0.5; 1];
    C4 = worldMatrix * [size(im,2)-0.5; sliceY-1; size(im,3)-0.5; 1];
    C2 = worldMatrix * [-0.5; sliceY-1; size(im,3)-0.5; 1];

    X = [C1(1) C2(1); C3(1) C4(1)];
    Y = [C1(2) C2(2); C3(2) C4(2)];
    Z = [C1(3) C2(3); C3(3) C4(3)];

    surf(X,Y,Z,'FaceColor','texturemap','Cdata',squeeze(im(sliceY,:,:,1)) * size(cmap,1),'CdataMapping', 'direct', 'EdgeColor', 'none')
end

if (showXYZ(3))
    C1 = worldMatrix * [-0.5; -0.5; sliceZ-1; 1];
    C2 = worldMatrix * [size(im,2)-0.5; -0.5; sliceZ-1; 1];
    C3 = worldMatrix * [-0.5; size(im,1)-0.5; sliceZ-1; 1];
    C4 = worldMatrix * [size(im,2)-0.5; size(im,1)-0.5; sliceZ-1; 1];

    X = [C1(1) C2(1); C3(1) C4(1)];
    Y = [C1(2) C2(2); C3(2) C4(2)];
    Z = [C1(3) C2(3); C3(3) C4(3)];

    surf(X,Y,Z,'FaceColor','texturemap','Cdata',im(:,:,sliceZ,1) * size(cmap,1),'CdataMapping', 'direct', 'EdgeColor', 'none')
end


axis equal
colormap(cmap)

if (~held)
    hold off
end


end

