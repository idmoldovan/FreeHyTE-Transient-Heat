function [Loops,Dim] = AssignPartsPart(Loops)
% ASSIGNPARTSPART maps the (particular solution) solving system and assigns
% each element an entry point and a dimension.
%
% ASSIGNPARTSCOMPL is called by MAINREG and MAINTRI. 
% Input: 
%  receives structure Loops as input argument, 
% Return: 
%  Loops structure, updated with two fields containing the insertion points
%  of the respective blocks in the solving system, and their dimensions. It
%  also returns Dim, which is the dimension of the solving system.
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
% The particular solution blocks corresponding to each finite element are
% completely disjunct (which is why each block is solved separately). Each
% block has as many lines as collocation points and as many columns as
% functions is the particular solution basis. 
%
% Further details on the construction of the particular solution system can
% be found in Section 6.1 of Reference [2] and Section 4.3 of Reference
% [4].

%% Initialization
% Initializes the current line and column indices.
entryline = 1;
entrycol = 1;

% Initialize the insertion points for the particular solution blocks
Loops.insertP = zeros(length(Loops.area),2); % line & column insertion points
% Initialize the dimensions for the particular solution blocks
Loops.dimP = zeros(length(Loops.area),2);   % number of lines is the number of collocation points;
                                            % number of columns is the number of functions in the collocation basis

%% Mapping of the entry points and dimensions
% The insertion points and dimensions are computed for each element and
% stored in the Loops structure. 
for i = 1:length(Loops.area)
    % insertion points
    Loops.insertP(i,:) = [entryline entrycol];
    % the total number of functions in the particular solution basis is the
    % sum of the orders of the basis for each generalized wave number
    % lambda, plus 3.
    Loops.dimP(i,2) = 2*Loops.orderP(i,1)+2*Loops.orderP(i,2)+...
        2*Loops.orderP(i,3)+3; 
    % Number of collocation points
    Loops.dimP(i,1) = Loops.gc(i)^2; 
    % Incrementing the entry line and column
    entryline = entryline + Loops.dimP(i,1);
    entrycol = entrycol + Loops.dimP(i,2);
end

%% Dimension of the particular solution system
Dim = [entryline-1 entrycol-1];

end

