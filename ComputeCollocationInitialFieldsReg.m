function [Loops] = ComputeCollocationInitialFieldsReg(Nodes,Loops,Ini)
% COMPUTECOLLOCATIONINITIALFIELDSREG computes the values of the initial
% fields in the collocation points and stores them in the Loops structure.
%
% COMPUTECOLLOCATIONINITIALFIELDSREG is called by MAINREG. 
% Input: 
%  receives structures Nodes and Loops, and handle INI to the the function
%  that calculates the initial fields.
% Return: 
%  Loops structure, updated with the initial field values in the
%  collocation points.
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
% COMPUTECOLLOCATIONINITIALFIELDSREG computes the initial fields in the
% FIRST grid.


%% Start iteration on elements
for ii=1:length(Loops.area)
    
    % LocLoop is a local structure where the features of the current
    % element which are directly useful for the calculation of the
    % collocation points are stored.
    LocLoop = struct('id',ii,'nodes',Loops.nodes(ii,:), 'center',...
        Loops.center(ii,:), 'dim', Loops.dimP(ii,:));
    
    % Initialize the initial temperature and 'velocity' fields
    Loops.T0(ii) = {zeros(LocLoop.dim(1),1)};
    Loops.V0(ii) = {zeros(LocLoop.dim(1),1)};

    %% Generating the geometric data
    % The following code computes the global coordinates of the collocation
    % points in the current element.
    
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
        
    % Compute the collocation points in the global referential
    GlobalX = X + LocLoop.center(1); 
    GlobalY = Y + LocLoop.center(2);
    
    %% Computing the initial fields
    % Calling INI
    [T0,V0] = Ini(GlobalX(:),GlobalY(:));
    
    % Storing the fields in the Loops structure
    Loops.T0(ii) = {T0};
    Loops.V0(ii) = {V0};
 
end

end