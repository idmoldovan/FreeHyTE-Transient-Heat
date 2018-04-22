function [Loops] = ComputeGaussInitialFieldsTri(Nodes,Loops,Ini,NGP)
% COMPUTEGAUSSINITIALFIELDSTRI computes the values of the initial
% fields in a grid of NGPxNGP Gauss points and stores them in the Loops
% structure. Used for plotting purposes.
%
% COMPUTEGAUSSINITIALFIELDSTRI is called by MAINTRI.
% Input:
%  receives structures Nodes and Loops, handle INI to the the function
%  that calculates the initial fields, and NGP, the number of Gauss points.
% Return:
%  Loops structure, updated with the initial field values in the plotting
%  points.
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
% There are two grids of Gauss points in the domain of each element:
% * the collocation grid, formed by Loops.gc(i)^2 Gauss-Legendre points,
% and used for the computation of the particular solution;
% * the solution plotting grid, formed by NGP^2 Gauss-Legendre points, and
% used for plotting purposes.
% COMPUTEGAUSSINITIALFIELDSTRI computes the initial fields in the
% SECOND grid.

%% Start iteration on elements
for ii=1:length(Loops.area)

    % LocLoop is a local structure where the features of the current
    % element which are directly useful for the calculation of the
    % plotting points are stored.
    LocLoop = struct('id',ii,'nodes',Loops.nodes(ii,:), 'center',...
        Loops.center(ii,:));
    
    % Initialize the initial temperature and 'velocity' fields
    Loops.GT0(ii) = {zeros(NGP)};
    Loops.GV0(ii) = {zeros(NGP)};
 
    %% Generating the geometric data
    % The following code computes the global coordinates of the plotting
    % points in the current element.
    
    % Getting coordinates of the nodes of the element (global)
    LocNodes = Nodes(LocLoop.nodes(:),:);
    % Generating the Gauss points (global)
    [GlobalX,GlobalY,~,~]=triquad(NGP,LocNodes);
    
    %% Computing the initial fields
    % Calling INI
    [T0,V0] = Ini(GlobalX,GlobalY);
    
    % Storing the fields in the Loops structure
    Loops.GT0(ii) = {T0};
    Loops.GV0(ii) = {V0};
    
end

end