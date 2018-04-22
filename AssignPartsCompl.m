function [Edges,Loops,Dim] = AssignPartsCompl(Edges, Loops)
% ASSIGNPARTSCOMPL maps the (complementary solution) solving system 
% and assigns each  block an entry point and a dimension.
%
%
% ASSIGNPARTSCOMPL is called by MAINREG and MAINTRI. 
% Input: 
%  receives structures Edges and Loops as input arguments, 
% Return: 
%  Edges and Loops structures, updated with two fields containing the 
%  insertion points of the respective blocks in the solving system, 
%  and their dimensions. It also returns Dim, which is the dimension 
%  of the solving system.
%
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
% The layout of a solving system with a single finite element a single
% essential boundary is presented below.
%
% The matrix of coefficients of the solving sistem  described in section
% 4.4 of reference [4], is constant in time. Therefore, it is computed
% only once, before the time stepping process is started. The general
% layout of the (complementary solution) solving system is given below. The
% blocks with index 'c' are specific to the complementary solution. 
%   __Matrix of coefficients__
%   |       |       |       | 
%   |   Kc  |  -Bc  |  -Bc  | 
%   |_______|_______|_______| 
%   |       |       |       | 
%   |  -Bc  |   H   |   0   | 
%   |_______|_______|_______| 
%   |       |       |       | 
%   |  -Bc  |   0   |   0   | 
%   |_______|_______|_______| 
%   _________________
%
% The Kc matrix belonging to element i is inserted on the main 
% diagonal of the solving system at position Loops.insert(i). 
% Its dimension is equal to the order of the domain basis of the respective 
% element, Loops.order(i) ( 2*order of the approximation basis +1). 
%
% Each boundary 'j' has Bc block. Interior boundaries always have approximations.
% Exterior Dirichlet boundaries have approximations where temperatures and Robin boundaries 
% where temperature and normal heat flux relation are enforced.  
% The H is the convection matrix of edge j where where Robin boundary conditions are 
% imposed. The H is inserted in diagonal, at entry line and column corresponding to boundary j 
% For more insight on how to model different types of 
% boundary conditions, please consult Section 4.6 of reference [3].
% The Bcj blocks are inserted on the lines that correspond to the entries of
% the complementary conductivity block of the element(s) the boundary belongs to. 
% The entry line and column of the Bc block corresponding to boundary j and
% neighbouring element i are Loops.insert(i) and Edges.insert(j).
% The Bcj blocks are not square. They have as many lines as the shape
% functions of the neighbouring elements and Edges.dim(i) = Edges.orderD(i)+1 
% ( or Edges.orderR(i)+1) lines.
%
% The layout and storage of the solving system is discussed at length in
% reference [2] (Section 6.2).
% 

%% Initialization
% Initializes the current line indicator, entry.
entry = 1;

% Initialize the insertion points for the complementary solution parts
Loops.insertC = zeros(length(Loops.area),1); 
% initialize the dimensions of the complementary solution parts
Loops.dimC = zeros(length(Loops.area),1);

%% Mapping of the Kc (complementary) conductivity entry blocks of the system
% The insertion point of the Kc block and its dimension are
% computed for each element and stored in the Loops structure.
for i = 1:length(Loops.area)
    Loops.insertC(i) = entry;
    Loops.dimC(i) = 2*Loops.orderC(i)+1;
    entry = entry + Loops.dimC(i);
end

% Initialize the insertion points for the edges
Edges.insert = zeros(length(Edges.type),1); 
% Initialize the dimensions of the static and dynamic parts
Edges.dim = zeros(length(Edges.type),1); 

%% Mapping of the boundary entry blocks (for Dirichlet and Robin) of the system
% The insertion point of the boundary block and its dimension are
% computed for each edge and stored in the Edges structure.
for i = 1:length(Edges.insert)
    if strcmpi(Edges.type(i),'D') || strcmpi(Edges.type(i), 'R')
        Edges.insert(i) = entry;
        Edges.dim(i) = Edges.order(i)+1;
        entry = entry + Edges.dim(i);
    end
end  
% note that Neumann edges correspond to zero dim, zero insert point

% Computing the total dimension of the solving system
Dim = entry-1;

end
