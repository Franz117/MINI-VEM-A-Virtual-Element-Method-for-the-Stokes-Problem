function [P,PX,PY] = get_base_unvem(glob_coords,glob_xq,glob_yq)
%superexpensive function  but unvirtualize VEM element

    glob_xv = glob_coords(:,1);
    glob_yv = glob_coords(:,2);
    all_points = [glob_xq,glob_yq; glob_xv , glob_yv];

    dt = delaunayTriangulation(all_points);
    ne = dt.size(1);
    dim = size(dt.Points,1);

    A = sparse(dim,dim);
    F = zeros(dim,1);

    CL = dt.ConnectivityList;
    PTS = dt.Points;

    for ie=1:ne 

        dof = CL(ie,:);
        coords = PTS(dof,:);
        
        % area and h^2
        xv = coords(:,1);
        yv = coords(:,2);
        dx = xv([2 3 1]) - xv([1 2 3]);
        dy = yv([2 3 1]) - yv([1 2 3]);
        h2 = max(dx.^2 + dy.^2);
        area = 0.5 * abs( ...
        (xv(2)-xv(1))*(yv(3)-yv(1)) - ...
        (xv(3)-xv(1))*(yv(2)-yv(1)) );

        % gradients
        P = [yv([2 3 1]) - yv([3 1 2])] / (2*area);
        Q = [xv([3 1 2]) - xv([2 3 1])] / (2*area);

        % assembly
        A(dof,dof) = A(dof,dof) +(P*P.' + Q*Q.') * area;
        F(dof)=F(dof) + area/(3*h2);
    end


    % border condition
    nvi = size(glob_xq,1);
    int_dof = 1:nvi;
    bd_dof = nvi+1:dim;
    nv = numel(bd_dof);
    BD = sparse(dim, nv);
    for k = 1:nv
        BD(bd_dof(k), k) = 1;
    end

    % ball
    Ub = zeros(dim,1);
    Fmod = F;                    
    Ub(int_dof) = A(int_dof,int_dof) \ Fmod(int_dof);

    %vertices
    Fmod = -A * BD;
    Uint = A(int_dof,int_dof) \ Fmod(int_dof,:);
    Uvert = BD;
    Uvert(int_dof,:) = Uint;
    P = [Uvert, Ub];

%     % plot
%     for k=1:nv+1
%         figure
%         trisurf(dt.ConnectivityList, ...
%                 all_points(:,1), ...
%                 all_points(:,2), ...
%                 full(P(:,k)));
%         shading interp
%     end
    

    % GRADIENT with FDM: PX and PY using P
    tri = dt.ConnectivityList;   % M x 3 triangles
    pts = dt.Points;             % N x 2 coordinates
    numPts = dim;
    numFuncs = nv+1;
    
    gradP = zeros(numPts, 2, numFuncs);  % N x 2 x 6
    counts = zeros(numPts,1);            % triangle counts per vertex
    
    for t = 1:size(tri,1)
        verts = tri(t,:);
        A = [pts(verts,:), ones(3,1)];   % 3x3 for linear plane
        coeff = A\P(verts,:);            % 3x6: each col = [a;b;c] for a function
        gradTri = coeff(1:2,:)';          % 6 x 2 gradient
        
        % accumulate to triangle vertices
        for k = 1:numFuncs
            gradP(verts,:,k) = gradP(verts,:,k) + repmat(gradTri(k,:),3,1);
        end
        counts(verts) = counts(verts) + 1;
    end
    % average per vertex
    gradP = gradP ./ counts(:,[1 1],:);
    
    % select internal vertices only
    gradInternal = gradP(int_dof,:,:);   % size: numInternal x 2 x 6
    PX = squeeze(gradInternal(:,1,:));   % 36 x 7
    PY = squeeze(gradInternal(:,2,:));   % 36 x 7

    P = P(int_dof,:);


end


