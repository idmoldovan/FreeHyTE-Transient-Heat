function [RHS] =Gen_RHS_VectorReg(Nodes,Loops,InternHeat,Dim,...
    time,deltaT,theta)
% GEN_RHS_VECTORREG sweeps through the elements and calls the functions 
% that generate the free vector of the solving system for the
% particular solution.
%
% GEN_RHS_VECTORREG is called by MAINREG. 
% Input data: 
% the Nodes and Loops structures, vector Dim that stores the total number
% of lines and columns of the matrix of coefficients (as determined in
% ASSIGNPARTSPART), InternHeat, a handle to the function that computes the
% heat source, the current time 'time', the mid-point scheme calibration
% parameter 'theta', and the time step 'deltaT'.
% Output/returns to MAINREG:
% the matrix of coefficients of the particular solution solving system, RHS
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

% The particular solution is obtained by solving the collocation system 
% Fp(x) * Xp = f0(x)
% where Fp(x) is a matrix that collects the particular solution trial basis
% at the collocation points, Xp are the unknown weights of the functions
% collected in the basis, and f0(x) is the free vector, listing the values 
% of the source function in the collocation points. 
% Vector f0(x) is computed by GEN_RHS_VECTORREG.
% The procedure is detailed in Reference [2], Section 6.1, and Reference
% [4], Section 4.3. The expression of the source function is derived in
% Section 3 of Reference [4]. Its definition for a general parabolic
% equation is given in Section 5.2.1 of Reference [2].

%% Initialization
RHS = zeros(Dim(1),1);

%% Sweeping the elements
for ii=1:length(Loops.area)
    % for some issue related to the translation from cell to matrix, we had
    % to separate the single-element case from the rest.
    if length(Loops.area) == 1
        % LocLoop is a structure where the features of the current element
        % which are directly useful for the calculation of f0(x)are stored.
        LocLoop = struct('id',ii,'nodes',Loops.nodes(ii,:),'center',...
            Loops.center(ii,:),'order',Loops.orderP(ii,:),...
            'ohm',Loops.ohm(ii),...
            'insert',Loops.insertP(ii,:),'dim',Loops.dimP(ii,:),...
            'material',Loops.material(ii,:),...
            'T0',Loops.T0,'V0',Loops.V0);
    else
        LocLoop = struct('id',ii,'nodes',Loops.nodes(ii,:),'center',...
            Loops.center(ii,:),'order',Loops.orderP(ii,:),...
            'ohm',Loops.ohm(ii),...
            'insert',Loops.insertP(ii,:),'dim',Loops.dimP(ii,:),...
            'material',Loops.material(ii,:),...
            'T0',cell2mat(Loops.T0(ii)),'V0',cell2mat(Loops.V0(ii)));
    end
    
    % Computing the f0(x) vector of element ii
    RHSi = RHS_Vector(Nodes,LocLoop,InternHeat,time,deltaT,theta);
    
    % Inserting the vector in the global RHS vector at the loci determined
    % in ASSIGNPARTSPART
    RHS(LocLoop.insert(1):LocLoop.insert(1)+LocLoop.dim(1)-1) = RHSi;
    
end

end

function RHSi = RHS_Vector(Nodes,LocLoop,InternHeat,time,deltaT,theta)
% RHS_VECTOR is a local function that computes the values of the source
% function f0 in the collocation points.

%% Computing the collocation points in the global Cartesian referential
% Computing the length of the sides of the element in x and y direction.
% sze simply collects the distances between the two most far apart points
% of the element in x and y directions. ASSUMES THAT THE ELEMENT IS
% RECTANGULAR!!
sze = max(Nodes(LocLoop.nodes(:),:))-min(Nodes(LocLoop.nodes(:),:));

% The number of collocation points in each Cartesian direction
ncolpts = sqrt(LocLoop.dim(1));
% Generate Gauss-Legendre points, mapped on a [-1,1] interval
abscissa = gauleg(ncolpts, -1, 1); 
% Uncomment the code below to use Gauss-Chebyshev collocation points
% abscissa = (-1)*cos(pi * (2*(1:ncolpts)-1)/2/ncolpts); % Chebyshev points

% Compute the local x and y coordinates of the Gauss collocation points
x=1/2*sze(1)*abscissa;
y=1/2*sze(2)*abscissa;
[X,Y] = ndgrid(x,y); % matrices with the collocation points
X = X(:); % transform matrix to vector
Y = Y(:);

% Compute the global x and y coordinates of the Gauss collocation points
GlobalX = X + LocLoop.center(1);
GlobalY = Y + LocLoop.center(2);


%% Computing the source function f0 in the collocation points
% The source function f0 is derived in Section 3 of Reference [4].

% This is the part of the source function controlled by the initial
% conditions
t0 =  -1/LocLoop.material(4)*(1/(theta*deltaT).*LocLoop.T0(:) + ((1-theta)/theta) .* LocLoop.V0(:)) ; 

% This is the part of the source function controlled by the internal heat
% generation
Q = InternHeat(GlobalX,GlobalY,time);

% f0
f0 = -(1/LocLoop.material(1))*Q + t0;
RHSi = f0;

end

