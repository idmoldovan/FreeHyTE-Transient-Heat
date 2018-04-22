function RHS = Gen_Tf_Vector(Edges, BConds, RHS, time)
% GEN_TF_VECTOR iterates through the edges and calls the functions that 
% generate Tf block of the free vector (RHS) of the complementary solution
% solving system . The Tf blocks are only computed on the Robin
% (convection) boundaries. 
%
% GEN_TF_VECTOR is called by MAIN***.
% Input data:
%  * the Edges and BConds structures, the RHS vector (that is, the free
%  vector of the solving system), and the current time.
% Output/Returns to MAIN*** : 
% * the RHS vector with the Tf blocks of all elements inserted at the
% correct positions (as determined in ASSIGNPARTSCOMPL).
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
% GEN_TF_VECTOR computes the product between the (constant) temperature of
% the fluid surrounding the current boundary, and the boundary integral of
% the flux basis Z, defined by Chebyshev polynomials,
%        Z  = cos(n*ArcCos(abscissa))
% where n is the line of the current entry in vector Tf. The boundary
% integral is evaluated analytically.
% The fluid temperature is an arbitrary function of time, defined by any
% valid Matlab command in the GUI or as an input file (see Section 4.6 of
% reference [3]).  
%
% Further details on the structure of the solving system are presented in
% reference [4] (Section 4.4).

%% Sweeping the edges and selecting the exterior Robin boundaries
for ii=1:length(Edges.type)
    % Looking for external Robin boundaries
    if ( strcmpi(Edges.type(ii),'R') && Edges.lright(ii) == 0 )
        % LocEdge is a local structure where the features of the current
        % edge which are directly useful for the calculation of the
        % Tf block are stored.
        LocEdge =struct('id',ii,'nini',Edges.nini(ii),'nfin',Edges.nfin(ii),...
            'parametric',Edges.parametric(ii,:), 'lleft',Edges.lleft(ii),...
            'lright',Edges.lright(ii), 'order',Edges.order(ii),...
            'insert',Edges.insert(ii), 'dim',Edges.dim(ii));
                 
        % Computing the Tf vector of element ii. Function Tf_Vector_i
        % is a local function defined below.
        Tfi = Tf_Vector_i(LocEdge, BConds, time); 
        
        % Inserting the Tfi vector in the global RHS vector
        RHS(LocEdge.insert:LocEdge.insert+LocEdge.dim-1) = -Tfi;
    end

end

end

function Tfi = Tf_Vector_i(LocEdge, BConds, time)
% TF_VECTOR_I local function computes the Tf vector of the LocEdge
% boundary. The edge is mapped to a [-1,1] interval to perform the
% integration.

% n+1 denotes the current line of the block
n = 0:LocEdge.order;

%% Generating the geometric data
% Computing the length of the current edge
L = sqrt(LocEdge.parametric(3)^2 + LocEdge.parametric(4)^2); 

%% Computing the (constant) value of the fluid temperature
% assessing if the time variation of the temperature was defined in
% an input file or as an analytic expression
if isnumeric(BConds.RobinTime{LocEdge.id})   % from an input file
    Tt = BConds.RobinTime{LocEdge.id};
    % interpolates the input for the current time
    Tinterp = interpn(Tt(:,1),Tt(:,2),time);
    tg = Tinterp;
else   % analytic definition
    tg = eval(BConds.RobinTime{LocEdge.id}{1});
end

%% Computing the integral, analytically
% The integral of a Chebyshev function on a [-1,1] interval is,
Tfi = L/2*((-1).^n+1)./(1-n.^2)*tg;
% for n~=1, and zero if n==1. Setting the integral to zero for n==1:
Tfi(isnan(Tfi)) = 0;

end