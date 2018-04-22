function [LHS] = Gen_LHS_MatrixReg(Nodes,Loops,Dim,RegCode)
% GEN_LHS_MATRIXREG sweeps through the elements and calls the functions 
% that generate the matrix of coefficients of the solving system for the
% particular solution.
%
% GEN_LHS_MATRIXREG is called by MAINREG. 
% Input data: 
% the Nodes and Loops structures, vector Dim that stores the total number
% of lines and columns of the matrix of coefficients (as determined in
% ASSIGNPARTSPART), and RegCode, which is 1 if all elements are identical
% in all regards (geometry and apprximation basis), and 0 otherwise.
% Output/returns to MAINREG:
% the matrix of coefficients of the particular solution solving system, LHS
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
% collected in the basis, and f0(x) is a vector listing the values of the
% source functions in the collocation points. 
% Basis Fp(x) collects 3*(2*Np+1) trial functions of type,
% Fp = Jn(lambda*r).* exp(1i*Th.*n)
% where Np is the order of the basis and lambda are the algorithmic wave
% numbers supplied by the user in the GUI(1). Therefore, Fp has as many
% lines as collocation points and 3*(2*Np+1) columns. Therefore, it is 
% generally not a square matrix. The procedure is detailed in
% Reference [2], Section 6.1, and Reference [4], Section 4.3.

%% Initialization
LHS = zeros(Dim);

%% Sweeping the elements
for ii=1:length(Loops.area)
    % LocLoop is a structure where the features of the current element 
    % which are directly useful for the calculation of Fp(x)are stored.    
    LocLoop = struct('id',ii,'nodes',Loops.nodes(ii,:),'center',...
        Loops.center(ii,:),'order',Loops.orderP(ii,:),'lambda',...
        Loops.lambda(ii,:),'insert',Loops.insertP(ii,:),...
        'dim',Loops.dimP(ii,:));
    
    % If RegCode == 1 AND it is not the first element to compute the
    % coefficient bolck for, it reuses the same block over and over again.
    % Otherwise, it computes the blocks for every element.
    if ~RegCode || ~exist('LHSi','var')
        % Computing the Fp(x) matrix of element ii
        LHSi = LHS_Matrix(Nodes,LocLoop);
    end
    
    % Inserting the matrix in the global Fp matrix at the loci determined
    % in ASSIGNPARTSPART
    LHS(LocLoop.insert(1):LocLoop.insert(1)+LocLoop.dim(1)-1,...
        LocLoop.insert(2):LocLoop.insert(2)+LocLoop.dim(2)-1) = LHSi;
    
end

end


function LHSi = LHS_Matrix(Nodes,LocLoop)
% LHS_MATRIX is a local function that computes the values of the trial
% functions Fp in the collocation points.

%% Computing the collocation points in the polar referential
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
% abscissa = (-1)*cos(pi * (2*(1:ncolpts)-1)/2/ncolpts); 

% Compute the local x and y coordinates of the Gauss collocation points
x=1/2*sze(1)*abscissa;
y=1/2*sze(2)*abscissa;
[X,Y] = ndgrid(x,y); % matrices with the collocation points
X = X(:); % transform matrix to vector
Y = Y(:);

% Compute the local r and theta coordinates of the Gauss collocation points
R = sqrt(X.^2 + Y.^2);
T = atan2(Y,X);

%% Computing the shape functions
% For each lambda, there are 2*nsBi+1 shape functions
nsB1 = -LocLoop.order(1):LocLoop.order(1);
nsB2 = -LocLoop.order(2):LocLoop.order(2);
nsB3 = -LocLoop.order(3):LocLoop.order(3);
Z1 = LocLoop.lambda(1)*R;
Z2 = LocLoop.lambda(2)*R;
Z3 = LocLoop.lambda(3)*R;
LHSB1 = bsxfun(@besselj,nsB1,Z1) .* exp(1i*bsxfun(@times, T, nsB1));
LHSB2 = bsxfun(@besselj,nsB2,Z2) .* exp(1i*bsxfun(@times, T, nsB2));
LHSB3 = bsxfun(@besselj,nsB3,Z3) .* exp(1i*bsxfun(@times, T, nsB3));
% Concatenating the shape functions for the three lambdas
LHSi = cat(2,LHSB1,LHSB2,LHSB3);

end

