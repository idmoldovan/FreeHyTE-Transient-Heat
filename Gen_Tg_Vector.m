function RHS = Gen_Tg_Vector(Edges, BConds, RHS, abscissa, weight, time)
% GEN_TG_VECTOR iterates through the edges and calls the functions that 
% generate Tgamma block of the free vector (RHS) of the complementary
% solving system (The Tgamma vector is abbreviated in the code as tg). The
% Tg "blocks" are only computed on the exterior Dirichlet boundaries. For
% all other elements, the Tg blocks are filled with zeros.
%
% GEN_TG_VECTOR is called by MAIN***.
% Input data:
%  * the Edges, Loops and BConds structures, the RHS vector
% (that is, the free vector of the solving system), and the Gauss-Legendre
% integration parameters abscissa and weights, and the current time.
% Output/Returns to MAIN*** :
% * the RHS vector with the Tg blocks of all elements inserted at the
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
% * the boundary temperature function is defined as, 
% Tg(x,t) = tg(x) * Tt(time)
% tg(x) is a polynomial function, defined here through its values in an
% arbitrary number of equally spaced points on the boundary. A polynomial
% interpolation is performed between these values to obtain the analytic
% expression of the applied temperature. The degree of the polynomial is
% equal to the number of temperature values, minus one. Tt(time) is an
% arbitrary function of time, defined by any valid Matlab command in the
% GUI or as an input file (see Section 4.6 of reference [3]). 
%
% Further details on the structure of the solving system are presented in
% reference [4] (Section 4.4).


%% Sweeping the edges and selecting the exterior Dirichlet boundaries
for ii=1:length(Edges.type)
    
    % Exterior Dirichlet boundaries have no right element
    if ( strcmpi(Edges.type(ii),'D') && Edges.lright(ii) == 0 )
        % LocEdge is a local structure where the features of the current
        % edge which are directly useful for the calculation of the
        % tg block are stored.
        LocEdge =struct('id',ii,'nini',Edges.nini(ii),'nfin',Edges.nfin(ii),...
            'parametric',Edges.parametric(ii,:), 'lleft',Edges.lleft(ii),...
            'lright',Edges.lright(ii), 'order',Edges.order(ii),...
            'insert',Edges.insert(ii), 'dim',Edges.dim(ii));
                 
        % Computing the Tg vector of element ii. Function Tg_Vector_i
        % is a local function defined below.
        Tgi = Tg_Vector_i(LocEdge, BConds, abscissa, weight, time); 
        
        % Inserting the Tgi vector in the global RHS vector
        RHS(LocEdge.insert:LocEdge.insert+LocEdge.dim-1) = -Tgi;
    end

end

end

function Tgi = Tg_Vector_i(LocEdge, BConds, abscissa, weight, time)
% TG_VECTOR_I local function computes the Tg vector of the LocEdge exterior
% Dirichlet boundary. The edge is mapped to a [-1,1] interval to perform 
% the integration.

% n+1 denotes the current line of the block
n = 0:LocEdge.order;

%% Generating the geometric data
% Computing the length of the current edge
L = sqrt(LocEdge.parametric(3)^2 + LocEdge.parametric(4)^2); % length

%% Computing the integrands at the integration points
% Computing the values of the normal flux basis. Z* -> the order is 'n'
Zstar = conj(cos(bsxfun(@times,n,acos(abscissa))));
Zstar = Zstar.';

% Computing the values of the enforced temperatures at the abscissas:
% obtaining the equally spaced points on [-1,1] interval where the
% temperatures are defined and stored in BConds.Dirichlet
a = linspace(-1,1,length(BConds.Dirichlet{LocEdge.id}));

% obtaining the polynomial that has the values given in BConds.Dirichlet
% at the points a
if (isnan(BConds.Dirichlet{LocEdge.id}))
    error('local:consistencyChk',...
        'No Dirichlet boundary conditions are defined on edge %d. \n',...
        LocEdge.id);
else
    pol = polyfit(a,BConds.Dirichlet{LocEdge.id},...
        length(BConds.Dirichlet{LocEdge.id})-1);
end

% computing the values of the "pol"  (interpolation polynomials) at the
% abscissas 
tg = polyval(pol,abscissa);

% Assessing if the time variation of the temperature was defined in
% an input file or as an analytic expression

if isnumeric(BConds.DirTime{LocEdge.id})   % from an input file
    Tt = BConds.DirTime{LocEdge.id};
    % interpolates the input for the current time
    Tinterp = interpn(Tt(:,1),Tt(:,2),time);
    tg = Tinterp*tg.';
else   % analytic definition
    tg = eval(BConds.DirTime{LocEdge.id}{1})*tg.';
end

%% Computing the integral on the side
% The integral is the internal product between the flux basis and the 
% applied temperature 
Tgi2D = bsxfun(@times, Zstar, tg); 

% computes the integral
Tgi = L/2 * sum(bsxfun(@times,Tgi2D,weight.'),2); 


end
