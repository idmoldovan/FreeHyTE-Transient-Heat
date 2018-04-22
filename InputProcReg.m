function [NGP, Nodes, Edges, Loops, BConds, RegCode, SpecOutput] = InputProcReg
% INPUTPROCREG is the input processing function for regular meshes.
%
% INPUTPROCREG is called by MAINREG. It reads the input data inserted in 
% the GUIs and organizes it in the data structures Edges, Loops and BConds.
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
% INPUT DATA (read from the *.mat files created by the GUIs):
% * from TrHeatStructDef: Nx, Ny (number of divisions in x and y) for
% regular mesh generation, SpecOutput (number of steps to plot the 
% results), k (thermal conductivity), c (specific heat), rho (density), 
% NumberGaussPoints (the number of Gauss-Legendre points for the side
% integration), EdgesOrder and LoopsOrderC/LoopsOrderP (the orders of the
% approximation bases on the essential edges of the mesh and
% complementary/particular solution in the finite elements), lambdai
% (i=1,2,3), the generalized wave numbers for the definition of the
% particular solution basis;
% * from TrHeatRegBC1: nodes, edges_nodes, edges_loops, loops_nodes and 
% loops_edges. These variables contain geometrical and topological 
% information regarding the mesh. For a full description of these
% variables, please consult Section 5.3 of reference [2];
% * from TrHeatRegBC2: edgesDirichlet, edgesNeumann, edgesRobin, dataDir, 
% dataNeu, dataRobin;
% * edgesDirichlet (edgesNeumann/edgesRobin): list with the Dirichlet
% (Neumann/Robin) edges;
% * dataDir (dataNeu/dataRobin): cell array with three columns and as many
% lines as Dirichlet (Neumann/Robin) edges. For each edge, it stores the
% edge's index, and the Dirichlet (Neumann/Robin) boundary conditions. 
% The boundary conditions are stored as strings. For further details
% regarding the definition of the boundary conditions, please refer to
% Section 4.6 of reference [3].
%
%
% OUTPUT DATA (to function MAINREG):
% * NGP is the number of Gauss points for the line integration;
% * Nodes is a (NNODE x 2) matrix, where NNODE is the number of nodes in 
% mesh. It stores the coordinates of each node;
% * Edges, Loops and BConds are data structures storing information
% on the edges, finite elements (loops) and boundary conditions,
% respectively. They are documented in Section 5.3 of reference [2];
% * RegCode, is equal to 1 if all finite elements are identical in terms of
% both shape and domain refinement and 0 otherwise. If RegCode=1, time is
% saved by only computing the left-hand side of the particular solution
% collocation system once;
% * SpecOutput, specifies at how many steps should the results be plotted as
% colour maps.

%% ------------------ User Input Zone ------------------

% Loading mesh information from STRUCTDEF GUI
load('TrHeatStructDef','Nx','Ny','SpecOutput');

% RegCode controls the regular/irregular meshes 
% RegCode =1 if ALL elements are identical, otherwise RegCode =0 
% The faulty definition of this coefficient may have serious consequences. 
% It is therefore not included in the options available in the GUI.
%
RegCode = 0; 
%
%% --------- Mesh & struct definition ---------
% Loading mesh information from STRUCTDEF GUI
load('TrHeatRegBC1','nodes','edges_nodes','edges_loops','loops_nodes',...
     'loops_edges');

% Definition of the mesh-related data Nodes, and data structures Edges and
% Loops. For a full description of these data structure, please refer to
% Section 5.3 of reference [2].
% 
Nodes = nodes;
Edges=struct('nini',edges_nodes(:,1),'nfin',edges_nodes(:,2),...
    'parametric',createLine(Nodes(edges_nodes(:,1),:),...
    Nodes(edges_nodes(:,2),:)),...
    'lleft',edges_loops(:,1),'lright',edges_loops(:,2),...
    'type',char(zeros(length(edges_nodes(:,1)),1)),'order',...
    zeros(length(edges_nodes(:,1)),1));

Edges.type(:) = 'D';           % all edges are predefined as Dirichlet 
Edges.order(:) = NaN;          % all degrees are predefined as NaN

Loops=struct('nodes',loops_nodes,'edges',loops_edges,...
    'center',zeros(length(loops_nodes(:,1)),2),...
    'area',zeros(length(loops_nodes(:,1)),1),...
    'orderP',zeros(length(loops_nodes(:,1)),3),...  % orders for the particular solution
    'lambda',zeros(length(loops_nodes(:,1)),3),...
    'orderC',zeros(length(loops_nodes(:,1)),1),...  % orders for the complementary solution
    'ohm',zeros(length(loops_nodes(:,1)),1),...
    'gc',zeros(length(loops_nodes(:,1)),1),...      % number of Gauss collocation points in each direction
    'material',zeros(length(loops_nodes(:,1)),4));

% Computation of the barycenter and area of each finite element:     
% It uses the POLYGONCENTROID function by David Legland.
for i=1:length(loops_nodes(:,1))
    [Loops.center(i,:),Loops.area(i)] = ...
        polygonCentroid(Nodes(loops_nodes(i,:),:));
    Loops.area(i)=abs(Loops.area(i)); 
end

%% MATERIAL DATA
% Loads the material data and stores them in the Loop structure. The
% material data is assumed uniform for all elements. However, users may
% overwrite the data loaded from the GUI to define elements with distinct
% material properties. An example is given below.

% Loading the material data
load('TrHeatStructDef','k','c','rho');

% Allocate the material characteristics defined in the GUI to all elements
Loops.material(:,1) = k; %thermal condutivity
Loops.material(:,2) = c; %specific heat
Loops.material(:,3) = rho; %density
% Computing the thermal diffusivity based on the input data
Loops.material(:,4) = k/c/rho; %thermal diffusivity (alpha)
% ... or overwrite the material characteristics manually
% Loops.material(5,1) = 1;      <------- Examples
% Loops.material(5,2) = 2;      <------- Examples
% Loops.material(5,3) = 3;      <------- Examples
% Loops.material(5,4) = Loops.material(5,1)/Loops.material(5,2)/Loops.material(5,3)

%% ------------------- USER-DEFINED AREA ---------------------
%% RUN CONTROL DATA
% Loading algorithmic, refinement and edge data
load('TrHeatStructDef','NumberGaussPoints','EdgesOrder','LoopsOrderC',...
    'LoopsOrderP','lambda1','lambda2','lambda3');
load('TrHeatRegBC2','edgesDirichlet','edgesNeumann','edgesRobin',...
    'dataDir','dataNeu','dataRobin');
NGP = NumberGaussPoints;   %Number of Gauss integration points per interval

%% EDGE TYPE DATA
% Registration of the Neumann edges, using edgesNeumann vector from the GUI
% It is recalled that all edges were predefined as Dirichlet.
% Users may overwrite the data loaded from the GUI to change the boundary
% types. 
if exist('edgesNeumann')
    for i=1:length(edgesNeumann)
        Edges.type(edgesNeumann(i),1) = 'N'; 
    end
end

% Registration of the Robin edges, using edgesRobin vector from the GUI
% It is recalled that all edges were predefined as Dirichlet.
% Users may overwrite the data loaded from the GUI to change the boundary
% types. 
if exist('edgesRobin')
    for i=1:length(edgesRobin)
        Edges.type(edgesRobin(i),1) = 'R'; 
    end
end

%% EDGE REFINEMENT DATA
% Allocate the refinement order defined in the GUI to all Dirichlet and Robin
% boundaries ...
Edges.order(Edges.type=='D') = EdgesOrder;
Edges.order(Edges.type=='R') = EdgesOrder;
% ... or overwrite orders manually, if you need to have different orders
% for different essential boundaries
% Edges.order(:) = 4;      % <------- Examples
% Edges.order(3:4) = 0;    % <------- Examples
% Edges.order(8:3:11) = 0; % <------- Examples


%% ELEMENT DATA
% Allocate the refinement orders defined in the GUI to all elements
% 1. Particular solution basis
Loops.orderP(:,1) = LoopsOrderP; % order of the first Bessel basis
Loops.orderP(:,2) = LoopsOrderP; % order of the second Bessel basis
Loops.orderP(:,3) = LoopsOrderP; % order of the third Bessel basis
Loops.lambda(:,1) = lambda1; % lambda for the first Bessel functions
Loops.lambda(:,2) = lambda2; % lambda for the second Bessel functions
Loops.lambda(:,3) = lambda3; % lambda for the third Bessel functions
% 2. Complementary solution basis
% order of the complementary solution basis
Loops.orderC(:) = LoopsOrderC; 
% ... or overwrite orders manually, if you need to have different orders
% for different elements
% Loops.orderP(1,:) = 3;  % <------- Examples
% Loops.orderC(1) = 9;   % <------- Examples

% Computing the number of Gauss collocation points in each direction
Loops.gc(:) = ceil(sqrt(3*(2*LoopsOrderP+1))); 

%% TEMPERATURE AND FLUX DATA 
% Boundary conditions can be described by polynomials of any order. The 
% definition of a boundary condition is made by specifying its values in 
% as many equally spaced points along the boundary as needed to define its
% polynomial variation. For further details regarding the definition of the 
% boundary conditions, please refer to Section 4.6 of reference [3].
%
% BConds data structure collects information regarding the boundary 
% conditions. Its members are cell arrays with as many lines as the 
% external boundaries of the structure. The values of the fluxes (or/and
% temperatures) enforced on the boundaries are stored in the Neumann (or
% Dirichlet or Robin) fields of the structure. NaN is stored in the Dirichlet 
% and in the Robin field of a Neumann boundary, and vice-versa.
%
% Initialization of the BConds structure
BConds=struct('Neumann',{cell(length(edges_nodes),1)},...
    'NeuTime',{cell(length(edges_nodes),1)},...
    'Dirichlet',{cell(length(edges_nodes),1)},...
    'DirTime',{cell(length(edges_nodes),1)},...
    'Robin',{cell(length(edges_nodes),1)},...
    'RobinTime',{cell(length(edges_nodes),1)});
BConds.Neumann(:) = {NaN};
BConds.NeuTime(:) = {NaN};   
BConds.Dirichlet(:) = {NaN};
BConds.DirTime(:) = {NaN};   
BConds.Robin(:) = {NaN};
BConds.RobinTime(:) = {NaN}; 

% Dirichlet boundary conditions are assumed to be of the form 
% U(x,t) = u(x) * T(time)
% u(x) is a polynomial function, defined here through its values in an
% arbitrary number of equally spaced points on the boundary. T(time) is an
% arbitrary function of time, defined by any valid Matlab command. 

% Dirichlet boundary conditions, as imported from the GUI and stored in the
% Dirichlet field of the structure, in the normal and tangential
% directions. Users may overwrite the data loaded from the GUI to change 
% the boundary conditions. 
if exist('dataDir')
    for i=1:size(dataDir,1)
        % loading the space variation data
        BConds.Dirichlet{dataDir{i,1},1}=str2num(dataDir{i,2});
        
        try    % loading the time variation data from a file
            BConds.DirTime{dataDir{i,1},1} = load(dataDir{i,3});
        catch  % otherwise assuming an analytic expression was input
            BConds.DirTime{dataDir{i,1},1} = {strrep(dataDir{i,3},'t','time')};
        end
    end
end

% overwrite boundary conditions manually
% BConds.Dirichlet(1:4) = {0};

% Neumann boundary conditions are assumed to be of the form 
% F(x,t) = f(x) * T(time)
% f(x) is a polynomial function, defined here through its values in an
% arbitrary number of equally spaced points on the boundary. T(time) is an
% arbitrary function of time, defined by any valid Matlab command. 

% Neumann boundary conditions and stored in the
% Neumann field of the BConds structure, in the normal and tangential
% directions. Users may overwrite the data loaded from the GUI to change 
% the boundary conditions. 
if exist('dataNeu')
    for i=1:size(dataNeu,1)
        % loading the space variation data
        BConds.Neumann{dataNeu{i,1},1}=str2num(dataNeu{i,2});
        
        try    % loading the time variation data from a file
            BConds.NeuTime{dataNeu{i,1},1} = load(dataNeu{i,3});
        catch  % otherwise assuming an analytic expression was input
            BConds.NeuTime{dataNeu{i,1},1} = {strrep(dataNeu{i,3},'t','time')};
        end
    end  
end
% overwrite boundary conditions manually
% BConds.Neumann(1) = {[0.5 0]};  

% Robin boundary conditions are assumed to be of the form Tf(time). 
% T(time) is an arbitrary function of time, defined by any valid Matlab 
% command or loaded from a file. The heat transfer coefficient is also 
% important in the definition of the Robin boundary condition.

% Robin boundary conditions stored in the
% Robin field of the BConds structure, in the normal and tangential
% directions. Users may overwrite the data loaded from the GUI to change 
% the boundary conditions. 
if exist('dataRobin')
    for i=1:size(dataRobin,1)
        % loading the heat transfer coefficient
        BConds.Robin{dataRobin{i,1},1}=str2double(dataRobin{i,2});
        
        try    % loading the time variation data from a file
            BConds.RobinTime{dataRobin{i,1},1} = load(dataRobin{i,3});
        catch  % otherwise assuming an analytic expression was input
            BConds.RobinTime{dataRobin{i,1},1} = {strrep(dataRobin{i,3},'t','time')};
        end
    end
end
    
end 
