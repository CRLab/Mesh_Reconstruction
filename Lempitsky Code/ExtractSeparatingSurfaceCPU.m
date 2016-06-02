function [faces,vertices] = ExtractSeparatingSurfaceCPU(V,niter)

%ExtractSeparatingSurface implements the procedure described in the submission
%INPUT: V - the binary volume (should be logical, 1 corresponding to the interior)
%       niter - number of iterations within the quadratic programming solver
%OUTPUT: the mesh components
%Running without the input argument  will apply the algorithm to the volume
%correspnding to a ball. Running without the output argument will display
%the original zero-isosurface and the resulting separating surface on
%screen.

%Victor Lempitsky, 2009
%LICENSE
%See license.txt for the terms of use.
%ACKNOWLEDGEMENT
%If you use the code for the publication, please acknowledge it by citing the following technical report or the subsequent publication:
%V. Lempitsky, "Surface Extraction from Binary Volumes with Higher-Order Smoothness", MSR-TR-2009-31 (also in submission).


    close all;

    if nargin < 2
        niter = 1000;
    end
    if nargin == 0
        V = CreateBall;
        %V = CreateCube;
    end


    perim = FastPerim(uint8(V)); % the B set - voxels adjacent to the other value

    %computing margin
    if exist('bwsqdist')
        margin = sqrt(double(bwsqdist(perim))); %use the Saito-Toriwaki algorithm, if available
    else
        margin = bwdist(perim,'euclidean'); %use the built-in MATLAB function (may be very slow)
    end

    clear perim;

    disp('Dist transform computed. Creating bands...');

    BAND_SIZE = 4;

    morphel = strel('arbitrary',cat(3,[0 0 0; 0 1 0; 0 0 0],...
        [0 1 0; 1 1 1; 0 1 0], [0 0 0; 0 1 0; 0 0 0]));
    

    %switch this on/off to enable the preservation of thin structures
    preserveThin = 0;
    if preserveThin
        thick = imdilate(imerode(V, morphel), ones([3 3 3]));
        thin = V & ~thick;
        margin(find(thin)) = 0.1;
        clear thick thin
    end      

    band = margin <= BAND_SIZE;
    tightBand = imerode(band, morphel);
    tightBand(1,:,:) = 0;
    tightBand(end,:,:) = 0;
    tightBand(:,1,:) = 0;
    tightBand(:,end,:) = 0;
    tightBand(:,:,1) = 0;
    tightBand(:,:,end) = 0;

    band = imdilate(tightBand, morphel);
    %tightBand now contains voxels for which we evaluate smoothness
    %band now contains voxels from the tightBand and the 6-ajacent ones

    index = zeros(size(band));
    where = find(band);
    %size(where)
    nband = size(where,1);
    index(where) = 1:nband;

    disp('Bands set. Creating matrices...');

    %turning margin into the lower and upper bounds on F
    lbnd = margin;
    lbnd(find(V == 0)) = -1000;%inf;
    ubnd = -margin;
    ubnd(find(V == 1)) = +1000;%inf;

    lb = lbnd(where);
    ub = ubnd(where);

    clear lbnd ubnd margin 

    %creating the sparse matrix, implementing smoothness
    mid = index(find(tightBand));
    ntight = size(mid,1);
    [x y z] = ind2sub(size(tightBand),find(tightBand));
    left = index(sub2ind(size(index),x-1,y,z));
    right = index(sub2ind(size(index),x+1,y,z));
    top = index(sub2ind(size(index),x,y-1,z));
    bottom = index(sub2ind(size(index),x,y+1,z));
    front = index(sub2ind(size(index),x,y,z-1));
    back = index(sub2ind(size(index),x,y,z+1));

    Hi = [1:ntight 1:ntight 1:ntight ...
        ntight+1:2*ntight ntight+1:2*ntight ntight+1:2*ntight...
        2*ntight+1:3*ntight 2*ntight+1:3*ntight 2*ntight+1:3*ntight]';

    Hj = [mid; left; right; mid; top; bottom; mid; front; back];

    clear mid left right mid top bottom mid front back tightBand index

    Hs = repmat([-2*ones(ntight,1); ones(ntight,1); ones(ntight,1)],3,1);

    H = sparse(Hi,Hj,Hs); %smoothness matrix created
    clear Hi Hj Hs
    
    H = H'*H;
 
    %solving convex quadratic program
    disp('Program prepared. Starting QP...');
    tic

    x = zeros(size(H,1),1);

    w = find(x < lb);
    x(w) = lb(w);
    w = find(x > ub);
    x(w) = ub(w);
    
    omega = 0.5;

    dg = diag(H);
    R = H - sparse(diag(dg));
    clear H
    invdg = 1.0./dg;
    clear dg

    x = RunQP(R, invdg, lb, ub, x, 1000.0);
    toc
    
    clear R dg x_ lb ub

    disp('Done, now extracting smooth surface...');
    V = 2*V-1;
    F = V*(BAND_SIZE+1);
    F(where) = x; 

    if(any(F(:).*V(:) < 0)) %check that the embedding function is consistent with the original volume
        error('Error in the algorithm');
    end

    clear V x
    %extracting the resulting surface
    if nargout > 0
        [faces,vertices] = isosurface(F,0); 
    else 
        figure;
        isosurface(F,0); 
        daspect([1 1 1]);
        title('Smoothed surface');
    end
    disp('Done');

end


function V = CreateCube
%creating the cube test example
    sz = 128;
    [x y z] = ndgrid(1:sz,1:sz,1:sz);
    x = x/max(x(:))-0.5;
    y = y/max(y(:))-0.5;
    z = z/max(z(:))-0.5;

    phi = 0.3;
    ksi = 0.2;

    mat1 = [cos(phi) sin(phi) 0; -sin(phi) cos(phi) 0; 0 0 1];
    mat2 = [cos(ksi) 0 sin(ksi); 0 1 0; -sin(ksi)  0 cos(ksi)];

    mat = mat1*mat2;

    x1 = x*mat(1)+y*mat(2)+z*mat(3);
    y1 = x*mat(4)+y*mat(5)+z*mat(6);
    z1 = x*mat(7)+y*mat(8)+z*mat(9);

    V = (x1 < 0.3 & x1 > -0.3 & y1 < 0.3 & y1 > -0.3 & z1 < 0.3 & z1 > -0.3);
    
end


function V = CreateBall
%creating the ball test example
    sz = 64;
    [x y z] = ndgrid(1:sz,1:sz,1:sz);
    x = x/max(x(:))-0.5;
    y = y/max(y(:))-0.5;
    z = z/max(z(:))-0.5;
    % 
    p = floor(sz*0.33);
    V = single((func(x,y,z) < func(x(p,p,p),y(p,p,p),z(p,p,p))));

    function r = func(x,y,z)
    r=x.*x+y.*y+z.*z;
    end
end
