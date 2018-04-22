function LHS = Gen_H_Matrix(Edges, LHS, BConds)
% GEN_H_MATRIX iterates through the edges and calls the functions that 
% generate the H block of the convection matrix in the LHS. 
%
% BGEN_H_MATRIX is called by MAIN***. It receives as input data the Edges 
% structure, the LHS matrix (that is, the matrix of coefficients of the
% solving system), BConds structure. It returns to MAIN*** the LHS matrix
% with the H blocks of all Robin boundaries, inserted at the correct
% positions (as determined in ASSIGNPARTSCOMPL). 
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
% GEN_H_MATRIX computes the internal product between the heat flux basis, 
%      Z  = cos(n*ArcCos(abscissa))
% with itself, on every Robin boundary of the mesh. All integrals are
% computed analytically.
%

%% Start iteration on edges
for ii=1:length(Edges.type)
    
    % Convection blocks are constructed for Robin boundaries only. 
    if (strcmpi(Edges.type(ii),'R'))
        % LocEdge is a local structure where the features of the current
        % edge which are directly useful for the calculation of the
        % convection block are stored.
        LocEdge = struct('id',ii,'nini',Edges.nini(ii),'nfin',Edges.nfin(ii),...
            'parametric',Edges.parametric(ii,:), 'lleft',Edges.lleft(ii),...
            'lright',Edges.lright(ii), 'order',Edges.order(ii),...
            'insert',Edges.insert(ii), 'dim',Edges.dim(ii));
        
        % Generating the convection block corresponding to the current edge
        Hi = H_Matrix_i(LocEdge, BConds);
        
        % Inserting the matrix on the diagonal of the global LHS matrix
        LHS(LocEdge.insert(1):LocEdge.insert(1)+LocEdge.dim(1)-1,...
            LocEdge.insert(1):LocEdge.insert(1)+LocEdge.dim(1)-1) = Hi; 
    end
     
end

end

function Hi = H_Matrix_i(LocEdge, BConds)
% H_MATRIX_I is a local function that computes the H block for the LocEdge
% edge. The integration is analytical.

%% Initialization 
% Retrieving the heat transfer coefficient 
if (isnan(BConds.Robin{LocEdge.id}))
    error('local:consistencyChk',...
        'No Robin boundary conditions are defined on edge %d. \n',...
        LocEdge.id);
else
    hc = BConds.Robin{LocEdge.id};
end

% Initializing the line & column indices
n = 0:LocEdge.order(1); 
m = 0:LocEdge.order(1);
[N,M] = ndgrid(n,m);

%% Generating the geometric data
% Computing the length of the current edge
L = sqrt(LocEdge.parametric(3)^2 + LocEdge.parametric(4)^2); 

%% Computing the integral
% Analytic expression (valid for M~=N)
Hi = (L/2/hc * (1./(1-(M+N).^2) + 1./(1-(M-N).^2))) .* (rem(M + N, 2)==0);
% Setting the entries to zero for N=M (the expression above yielded NaN)
Hi(isnan(Hi)) = 0;

end