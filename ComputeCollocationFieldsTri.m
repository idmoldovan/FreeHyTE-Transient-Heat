function [Loops] = ComputeCollocationFieldsTri(Nodes,Loops,X,Xp,deltaT,theta)
% COMPUTECOLLOCATIONFIELDSTRI computes the temperature and its derivative 
% in the collocation points, to be used as initial conditions in the next
% time step. The results are stored in the Loops structure.
%
% COMPUTECOLLOCATIONFIELDSTRI is called by MAINTRI. 
% Input:
% * It receives structures Nodes and Loops structures, the solution vectors 
% X (complementary solution) and Xp (particular solution), the time step
% deltaT, and the time discretization parameter theta.  
% Ouput:
% The Loops structures updated with the temperatures and flux values in
% the collocation points.
%
% BIBLIOGRAPHY
% 1. FreeHyTE Page - https://sites.google.com/site/ionutdmoldovan/freehyte
% 2. Moldovan ID, Cismasiu I - FreeHyTE: theoretical bases and developer’s 
% manual, https://drive.google.com/file/d/0BxuR3pKS2hNHTzB2N2Q4cXZKcGc/view
% 3. FreeHyTE Heat HTTE User's Manual - 
%    https://drive.google.com/drive/folders/0BxuR3pKS2hNHaFhiSjZHOE9TMzg
% 4. Moldovan ID, Coutinho AP, Cismasiu I - Hybrid-Trefftz finite elements
% for non-homogeneous parabolic problems using a novel dual reciprocity
% variant, Engineering Analysis with Boundary Elements, under review, 2018.
%
% The estimates of the temperature field are obtained by
% substituting the solutions X and Xp of the solving systems in the domain
% approximations of the respective fields. This is done element by element.
%
% The computation of the field estimates is further covered in Section 7.1
% of reference [2]. 
%

%% Computation of the final fields, to store in the collocation pts
% Iteration on the elements
for ii=1:length(Loops.area)
    %% Initialization
    % LocLoop is a local structure where the features of the current
    % element which are directly useful for the calculation of the output
    % fields are stored.
    LocLoop = struct('id',ii,'nodes',Loops.nodes(ii,:), ...
        'edges',Loops.edges(ii,:),'center',Loops.center(ii,:),...
        'orderP',Loops.orderP(ii,:),'orderC',Loops.orderC(ii),...
        'insertP',Loops.insertP(ii,2),'insertC',Loops.insertC(ii),...
        'dimP',Loops.dimP(ii,:),'dimC',Loops.dimC(ii),...
        'lambda', Loops.lambda(ii,:),'T0',cell2mat(Loops.T0(ii)),...
        'V0',cell2mat(Loops.V0(ii)),'material',Loops.material(ii,:),...
        'ohm',Loops.ohm(ii));
    
    % orders of the complementary and particular solution bases
    nd = -LocLoop.orderC:LocLoop.orderC;
    nsB1 = -LocLoop.orderP(1):LocLoop.orderP(1);
    nsB2 = -LocLoop.orderP(2):LocLoop.orderP(2);
    nsB3 = -LocLoop.orderP(3):LocLoop.orderP(3);
    
    %% Generating the geometric data
    % Getting coordinates of the nodes of the element (global)
    LocNodes = Nodes(LocLoop.nodes(:),:);
    
    % Computing the number of collocation points
    ncolpts = sqrt(LocLoop.dimP(1));
    % Getting the collocation points in global coordinates
    [GlobalX,GlobalY,~,~]=triquad(ncolpts,LocNodes);
    % Getting the collocation points in local coordinates
    x = GlobalX - LocLoop.center(1);
    y = GlobalY - LocLoop.center(2);
    % Getting the collocation points in polar coordinates
    r = sqrt(x.^2 + y.^2);
    t = atan2(y,x);
    
    
    %% Computing the basis functions - particular solution
    % Arguments of the basis functions
    z1 = LocLoop.lambda(1)*r;
    z2 = LocLoop.lambda(2)*r;
    z3 = LocLoop.lambda(3)*r;
    
    % Computing temperature shape functions in the collocation points.
    % For a full description of the basis, please refer to Section 4.3.2 of
    % reference [4].
    UPB1 = 1/(LocLoop.ohm^2-LocLoop.lambda(1)^2)*...
        bsxfun(@besselj,nsB1,z1(:)) .* exp(1i*bsxfun(@times, t(:), nsB1));
    UPB2 = 1/(LocLoop.ohm^2-LocLoop.lambda(2)^2)*...
        bsxfun(@besselj,nsB2,z2(:)) .* exp(1i*bsxfun(@times, t(:), nsB2));
    UPB3 = 1/(LocLoop.ohm^2-LocLoop.lambda(3)^2)*...
        bsxfun(@besselj,nsB3,z3(:)) .* exp(1i*bsxfun(@times, t(:), nsB3));
    % Concatenating the basis
    UP = cat(2,UPB1,UPB2,UPB3);
    % Extracting the weights corresponding to the current element,
    % according to the insertion point and dimensions allocated for the
    % current element in ASSIGNPARTSPART.
    XP = Xp(LocLoop.insertP:LocLoop.insertP+LocLoop.dimP(2)-1);
    
    %% Computing the basis functions - complementary solution
    % Arguments of the basis functions
    z = LocLoop.ohm*r;
    
    % Computing temperature shape functions (complementary solution)
    UC = bsxfun(@besselj,nd,z(:)) .* exp(1i*bsxfun(@times,nd,t(:)));
    % Extracting the weights corresponding to the current element,
    % according to the insertion point and dimensions allocated for the
    % current element in ASSIGNPARTSCOMPL.
    XC = X(LocLoop.insertC:LocLoop.insertC+LocLoop.dimC-1);
    
    %% Computing the temperature and its time derivative
    T = UP*XP + UC*XC;
    V = 1/(theta*deltaT) .* (T - LocLoop.T0(:)) - (1-theta)/ theta .* LocLoop.V0(:);
    
    % Storing the fields in the Loops data structure.
    Loops.T0(ii) = {T};
    Loops.V0(ii) = {V};
    
    
end
end