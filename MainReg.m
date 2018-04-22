function MainReg
% Main processor program
% MainReg is the main routine of FreeHyTE - Transient Heat Conduction 
% module with regular rectangular meshes.
% MAINREG is called upon exiting the data input (GUI) phase of the module.
% It is used to launch all functions required for the solution of the 
% heat transfer problem and centralize all data they provide. 
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

%% ---------- User input ------------
% Launch pre-processing routine: INPUTPROCREG and TIMEGENREG
% Processes the input data from GUI, generates the mesh and topological
% data structures
% * NGP is the number of Gauss points for the line integration;
% * Nodes is a (NNODE x 2) matrix, where NNODE is the number of nodes in 
% mesh. It stores the coordinates of each node;
% * Edges, Loops and BConds are data structures storing information
% on the edges, finite elements (loops) and boundary conditions,
% respectively. They are documented in Section 5.3 of reference [2];
% * RegCode is equal to 1 if all finite elements are identical in terms of
% both shape and domain refinement and 0 otherwise. If RegCode=1, time is
% saved by only computing the left-hand side of the particular solution
% collocation system once;
% * SpecOutput specifies at how many steps should the results be plotted as
% colour maps.
[NGP,Nodes,Edges,Loops,BConds,RegCode,SpecOutput] = InputProcReg; 

% Launch time domain data processor. TIMEGENREG computes the generalized
% frequency for each element, stores it in the Loops structure, and returns
% handles to the functions that compute the initial fields (Ini) and the
% heat source (InternHeat). The other output of TIMEGENREG is listed below:
% * TotalTime is the total duration of the analysis;
% * deltaT is the time step;
% * theta is the time integration coefficient.
[Loops, TotalTime, deltaT, Ini, InternHeat, theta] = TimeGen(Loops);
% ---------- End of user input zone ------------

%% --------- Pre-flight initialization -----------

% Initialization of Gauss weights & abscissas (on a -1:1 interval).
% gauleg is a external routine, written by Greg von Winckel.
%
[abscissa,weight] = gauleg(NGP, -1, 1);
%
% ASSIGNPARTS**** functions map the solving systems for the particular and
% complementary solutions and assigns entries and dimensions to elements
% and sides. The information is used by the functions that generate the
% blocks of the solving systems to insert them at the right positions.
% * Dim**** is the total dimension of the finite element solving systems;
% * the computation of the particular solution in FreeHyTE is documented in
% Section 4.3 of reference [4];
% * the computation of the complementary solution is covered in Section 4.4
% of reference [4]. 

% Mapping the solving system for the computation of the particular solution
%
[Loops,DimPart] = AssignPartsPart(Loops);  
% Mapping the solving system for the computation of the complementary
% solution 
[Edges,Loops,DimCompl] = AssignPartsCompl(Edges, Loops);  

% Computing the initial fields in the collocation points
Loops = ComputeCollocationInitialFieldsReg(Nodes,Loops,Ini); 
% Computing the initial fields in the Gauss points, for plotting
Loops = ComputeGaussInitialFieldsReg(Nodes,Loops,Ini,abscissa); 
%
%% --------- Calculation of the constant matrices -----------
% The matrix of coefficients of the solving sistem  described in section
% 4.4.4 of reference [4], is constant in time. Therefore, it is computed
% only once, before the time stepping process is started. The general
% layout of the (complementary solution) solving system is given below. The
% blocks with index 'c' are specific to the complementary solution. The
% blocks with index 'p' refer to the particular solution. Terms qg, tf and
% tg correspond to the enforced boundary conditions.

%   __Matrix of coefficients______________free vector_
%   |       |       |       |  <------>   |         |
%   |   Kc  |  -Bc  |  -Bc  |             | qg-KpXp |
%   |_______|_______|_______|             |_________|
%   |       |       |       |             |         |
%   |  -Bc  |   H   |   0   |             | -tf+BpXp|
%   |_______|_______|_______|             |_________|
%   |       |       |       |             |         |
%   |  -Bc  |   0   |   0   |             | -tg+BpXp|
%   |_______|_______|_______|             |_________|

% Initialization of the matrix of coefficients
LHS = zeros(DimCompl);

% Generating & allocating the LHS matrices
% The following functions generate the coefficients blocks for each finite
% element and essential boundary and insert them in LHS at the right place,
% according to the mapping information generated in ASSIGNPARTSCOMPL. 
% The explicit expressions of the conductivity and boundary matrices are 
% given in Section 4.4 of reference [4].
LHS = Gen_Kc_Matrix(Edges, Loops, LHS, abscissa, weight);
LHS = Gen_Bc_Matrix(Edges, Loops, LHS, abscissa, weight); 
LHS = Gen_H_Matrix(Edges, LHS, BConds);                     

% Preconditioning and solving the system
% System scaling procedure:
% Generating the scaling matrix, Sc. Sc is a diagonal matrix, whose terms
% are defined as the inverse of the square roots of the diagonal terms of
% the coefficient matrix.
Sc = (diag(LHS)).^(-1/2);
% If a diagonal term is null, the corresponding line and column are not
% scaled.
Sc(isinf(Sc)) = 1.0;
Sc = diag(Sc);
% Scaling LHS. ScLHS is the scaled version of the matrix of coefficients.
ScLHS = Sc' * LHS * Sc;
% The reciprocal condition number CndNo of the scaled LHS:
CndNo=rcond(ScLHS);
% Deciding if the conditioning of the LHS matrices permits storing them 
% in the LU form. If the reciprocal condition number is larger than the
% precision of the machine, storing the ststem in the LU form should speed
% up the solution of the system. Conversely, if the Moore-Penrose
% pseudoinverse "pinv" needs to be used to solve the system, one gains
% nothing by storing it in the LU form.
if (CndNo>1*eps)
    [L,U,p] = lu(ScLHS,'vector');
else
    PinvLHS = pinv(ScLHS);
end

%% **********************************************************************
% TIME LOOPING ZONE STARTS
% ***********************************************************************
SubDir = '';
for time=deltaT:deltaT:TotalTime
    
    fprintf('Time step %d \n', round(time/deltaT));
    
    %% --------- Computation of the particular solution ------------
    % The a new form of the Dual Reciprocity Method (DRM) is used for the
    % computation of the particular solution. As opposed to the mainstream
    % DRM variants, the approach used here presents a simple particular 
    % solution basis, which need not be singular. The method is described
    % in detail in Sections 2.2.3 and 2.3.3 of reference [2] and in Chapter
    % 4.3 of reference [4].    
    % The particular solution of the governing differential equation is 
    % computed by solving the following collocation system, one element at
    % a time,
    %
    %   ____________PartLHS____________             _PartRHS_
    %   |          |       |          |  <------>   |        |
    %   | F(1,-np) |  ...  | F(1, np) |             |  f0(1) |
    %   |__________|_______|__________|             |________|
    %   |          |       |          |             |        |
    %   |   ...    |  ...  |   ...    |             |  ...   |
    %   |__________|_______|__________|             |________|
    %   |          |       |          |             |        |
    %   |F(NCP,-np)|  ...  |F(NCP, np)|             | f0(NCP)|
    %   |__________|_______|__________|             |________|
    %
    % where F is the aproximation basis for the source term f0,
    % NCP is the number of collocation points, and np is the order of the
    % particular basis.
    
    % Generating & allocating the LHS blocks for the particular solution
    % solving system.
    PartLHS = Gen_LHS_MatrixReg(Nodes,Loops,DimPart,RegCode); 
    
    % Generating the particular RHS vector with the source term
    PartRHS = Gen_RHS_VectorReg(Nodes,Loops,InternHeat,DimPart,...
        time,deltaT,theta);                                    
    
    % Solving the particular solution systems, one element at a time
    Xp = zeros(size(LHS,2),1);
    for ii=1:length(Loops.area)
        % extracting the blocks of the PartLHS and PartRHS that correspond
        % to the current element, according to the mapping performed in
        % ASSIGNPARTSPART.
        LocLHS = PartLHS(Loops.insertP(ii,1) : Loops.insertP(ii,1) + ...
            Loops.dimP(ii,1)-1, Loops.insertP(ii,2) : Loops.insertP(ii,2)+...
            Loops.dimP(ii,2)-1);
        LocRHS = PartRHS(Loops.insertP(ii,1) : Loops.insertP(ii,1) + ...
            Loops.dimP(ii,1)-1);
        % Solving the particular solution system using the Moore-Penrose
        % pseudoinverse
        LocXp = pinv(LocLHS)*LocRHS;
        Xp(Loops.insertP(ii,2) : Loops.insertP(ii,2) + ...
            Loops.dimP(ii,2)-1)=LocXp;
    end 
    
    %% --------- Computation of the complementary solution ------------
    % After the computation of the particular solution, the free vector of
    % the complementary solution solving system (above) can be constructed.
    %
    % Generating & allocating the RHS vectors
    %
    % The following functions generate the free vectors for each finite
    % element and essential boundary and insert them at the right place,
    % according to the mapping information generated in ASSIGNPARTSCOMPL.
    %  ________
    % |        |                - qg is the enforced heat flow vector,
    % |qg-KpXp |                - Bp is the particular boundary matrix,
    % |________|                - Kp is the particular conductivity matrix,
    % |        |                - tg is the enforced temperature vector,
    % |-tf+BpXp|                - tf is the convection vector,
    % |________|                - Xp is the particular solution.
    % |        |
    % |-tg+BpXp|
    % |________|
    %
    % 
    % The explicit expressions of the RHS vectors are given in the section 
    % 4.4 of reference [4].
    %
    RHS = zeros(DimCompl,1);
    RHS = Gen_qg_Vector(Edges,Loops,BConds,RHS,abscissa,weight,time);
    RHS = Gen_KpXp_Vector(Edges,Loops,Xp,RHS,abscissa,weight);   
    RHS = Gen_Tg_Vector(Edges, BConds, RHS, abscissa, weight, time); 
    RHS = Gen_Tf_Vector(Edges, BConds, RHS, time);              
    RHS = Gen_BpXp_Vector(Edges,Loops,Xp,RHS,abscissa,weight);   
    
    % Solving the FE system
    % According to the pre-conditioning, the LHS vector has already been
    % scaled into ScLHS.
    if (CndNo<=1*eps)
        % Scaling RHS. ScRHS is the scaled versions of the free vector.
        ScRHS = Sc' * RHS;
        % Solving the scaled system by using the Moore-Penrose 
        % pseudoinverse-based solution procedure
        ScX=PinvLHS*ScRHS;
        % Scaling back the solution
        X = Sc * ScX;
    else
        % Scaling RHS. ScRHS is the scaled versions of the free vector.
        ScRHS = Sc' * RHS;
        % Solving the scaled system by using default solution procedure
        ScX = U\(L\ScRHS(p,:));
        % Scaling back the solution
        X = Sc * ScX;
    end
 
    
    % Constructing the displacement field
    [Loops] = ComputeCollocationFieldsReg(Nodes,Loops,X,Xp,deltaT,theta); 
    
    % Preparing the fields for plotting and/or storing
    [Loops,SubDir] = ComputeGaussFieldsReg(Nodes,Loops,X,...
            Xp,time,deltaT,SubDir,SpecOutput,abscissa,theta);   
    
    % Plotting the displacement field
    PlotFieldsReg(Nodes,Edges,Loops,time,deltaT,SpecOutput,abscissa); 
    
end
%% **********************************************************************
% TIME LOOPING ZONE ENDS
% **********************************************************************

fclose('all');

end
