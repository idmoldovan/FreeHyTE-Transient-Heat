function [Loops,SubDir] = ComputeGaussFieldsTri(Nodes,Loops,X,...
            Xp,time,deltaT,SubDir,SpecOutput,NGP,theta)
% COMPUTEGAUSSFIELDSTRI  computes the final temperature and flux fields, 
% in a grid of NGP x NGP points at the requested instants. The results are 
% stored in the Loops structure and in an output file, if requested by 
% the user (SpecOutput).
%
% COMPUTEGAUSSFIELDSTRI is called by MAINTRI. 
% Input:
% * It receives structures Nodes and Loops structures, the solution vectors 
% X (complementary solution) and Xp (particular solution), the current
% time, the time step deltaT, SubDir the name of the directory to store the
% results if required by the user, SpecOutput the number time steps the
% results should be stored, the number of Gauss point, NGP, and the time
% discretization parameter theta.
% Besides the input arguments,  COMPUTEGAUSSFIELDSTRI loads from
% TrHeatStructDef.mat information regarding the file where the results
% should be stored. 
% Ouput:
% The Loops structures updated with the temperatures and flux values in
% the (NGP+1)^2 plotting points (NGP is the number of Gauss-Legendre 
 % integration points), and SubDir, the name of the directory where the
 % solution is stored.
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
% The estimates of the temperature and flux fields are obtained by
% substituting the solutions X and Xp of the solving systems in the domain
% approximations of the respective fields. This is done element by element.
% The estimates are used for:
% * storing the temperature and flux estimates, calculated in the 
% (NGP+1)^2 plotting points, in the Loops structure returned to MAINTRI, 
%  for posterior plotting;
% * storing the temperature and flux estimates, calculated in the 
% (NGP+1)^2 plotting points in a result file which can be accessed by 
% third party post-processing software (it is pre-formatted for TecPlot). 
% To avoid overloading the disk with (fairly large) undesired result files,
% the result file is only created if the input was saved in the first GUI
% (see Section 5.2.2 of reference [3]).
%
% The computation of the field estimates is further covered in Section 7.1
% of reference [2]. The output options are explained in detail in Section
% 7.2 of reference [2]. 
%

%% Initialization and impoting information
% Computing the current time step number
tstep = round(time/deltaT);
% Loading the save file data from TrHeatStructDef.
% * DirName is the folder to save the result file;
% * FileName is the name of the file where the analysis data was saved.
load('TrHeatStructDef','DirName','FileName'); 
% if there's no DirName AND there's no SpecOutput (no plotting) or no 
% plotting is required at this time step...
if isempty(DirName) && (~SpecOutput || rem(tstep,SpecOutput))  
    return; % ...do nothing
end

%% Prepares the file to plot
% if the user left DirName blank, does not generate written output
if ~isempty(DirName)    
    % in the first step, it creates the subdirectory needed to save the
    % results
    if tstep == 1
        % creates a sub-folder with the current date and time
        status =  mkdir(DirName,FileName(1:end-4));
        if ~status
            error('local:Path',...
                'Cannot create a sub-folder in the working folder defined in the GUI. Exiting. \n');
        end
        % getting the sub-folder path
        SubDir = fullfile(DirName,FileName(1:end-4));
    end

    % Generating the name of the file corresponding to the current time
    % step. The file is called TS_numberofthetimestep
    SubFileName = sprintf('TS_%d.dat',tstep);
    % link the path to the filename
    % the UFilename string, which allows the creation of the result file
    UFilename = fullfile(SubDir,SubFileName);
    % open file for writing
    FileU = fopen(UFilename,'w');   
    
    % The following lines are the header needed to prepare the result file
    % to be readable with TecPlot. If you wish to use another software for
    % post-processing, you probably need to change the header.
    % However, the format of the data should be fairly general. The results
    % are written as a matrix of NEL*(NGP+1)^2 lines (where NEL is the 
    % total number of finite elements and (NGP+1)^2 the total number of 
    % result points per element), and 8 columns. For each result point, the
    % columns list the x and y coordinates, in the global referential, 
    % followed by the values of the temperature and flux fields (complex numbers).
    fprintf(FileU,'TITLE="%s"\n',UFilename);
    fprintf(FileU,'VARIABLES="X", "Y", "Re.U", "Im.U", "Re.Qx", "Im.Qx", "Re.Qy", "Im.Qy" \n');
    fprintf(FileU,'ZONE T="TS = %d"\n',tstep);

end

%% Computation of the final fields, to store in the grid pts
% Iteration on the elements
for ii=1:length(Loops.area)
    %% Initialization
    % LocLoop is a local structure where the features of the current
    % element which are directly useful for the calculation of the output
    % fields are stored.
    LocLoop = struct('id',ii,'nodes',Loops.nodes(ii,:), ...
        'edges',Loops.edges(ii,:),'center',Loops.center(ii,:),...
        'orderP',Loops.orderP(ii,:),'orderC',Loops.orderC(ii),...
        'insertP',Loops.insertP(ii,2),'insertC',Loops.insertC(ii),...
        'dimP',Loops.dimP(ii,2),'dimC',Loops.dimC(ii),...
        'lambda', Loops.lambda(ii,:),'GT0',cell2mat(Loops.GT0(ii)),...
        'GV0',cell2mat(Loops.GV0(ii)),'material',Loops.material(ii,:),...
        'ohm',Loops.ohm(ii));
    
    k = LocLoop.material(1);
    
    % orders of the complementary and particular solution bases
    nc = -LocLoop.orderC:LocLoop.orderC;
    npB1 = -LocLoop.orderP(1):LocLoop.orderP(1);
    npB2 = -LocLoop.orderP(2):LocLoop.orderP(2);
    npB3 = -LocLoop.orderP(3):LocLoop.orderP(3);
    
    % Getting coordinates of the nodes of the element (global)
    LocNodes = Nodes(LocLoop.nodes(:),:);
    
    % Generating the output points in the global Cartesian referential.
    % The output points belong to the Gauss-Legendre quadrature of the
    % triangular element.
    [globalx,globaly,~,~]=triquad(NGP,LocNodes);
    % Generating the output grid in local coordinates.
    x = globalx - LocLoop.center(1);
    y = globaly - LocLoop.center(2);
    
    % Transforming the local Cartesian coordinates of the Gauss collocation 
    % points into r and th polar coordinates
    r = sqrt(x.^2 + y.^2);
    t = atan2(y,x);
    
    %% Computing the basis functions - particular solution
    % Arguments of the basis functions
    z1 = LocLoop.lambda(1)*r;
    z2 = LocLoop.lambda(2)*r;
    z3 = LocLoop.lambda(3)*r;

    % Computing temperature shape functions in the grid points.
    % For a full description of the basis, please refer to Sections 4.2 and
    % 4.3 of reference [4].
    UPB1 = 1/(LocLoop.ohm^2-LocLoop.lambda(1)^2)*...
        bsxfun(@besselj,npB1,z1(:)) .* exp(1i*bsxfun(@times, t(:), npB1));
    UPB2 = 1/(LocLoop.ohm^2-LocLoop.lambda(2)^2)*...
        bsxfun(@besselj,npB2,z2(:)) .* exp(1i*bsxfun(@times, t(:), npB2));
    UPB3 = 1/(LocLoop.ohm^2-LocLoop.lambda(3)^2)*...
        bsxfun(@besselj,npB3,z3(:)) .* exp(1i*bsxfun(@times, t(:), npB3));
    % Concatenating the basis
    UP = cat(2,UPB1,UPB2,UPB3);
    % Extracting the weights corresponding to the current element,
    % according to the insertion point and dimensions allocated for the
    % current element in ASSIGNPARTSPART.
    XP = Xp(LocLoop.insertP:LocLoop.insertP+LocLoop.dimP-1);
    
    %% Computing the basis functions - complementary solution
    % Arguments of the basis functions
    z = LocLoop.ohm*r;
    
    % Computing temperature shape functions (complementary solution)
    UC = bsxfun(@besselj,nc,z(:)) .* exp(1i*bsxfun(@times,nc,t(:)));
    % Extracting the weights corresponding to the current element,
    % according to the insertion point and dimensions allocated for the
    % current element in ASSIGNPARTSCOMPL.
    XC = X(LocLoop.insertC:LocLoop.insertC+LocLoop.dimC-1);
    
    %% Computing the temperature and its time derivative
    T = UP*XP + UC*XC;
    V = (1/(theta*deltaT)) .* (T - LocLoop.GT0(:)) - ((1-theta)/ theta) .* LocLoop.GV0(:);
    T = reshape(T,[NGP,NGP]);
    V = reshape(V,[NGP,NGP]);
    
    %% Computing the heat fluxes, in radial and tangential directions
    % Computing the flux shape functions, in the R direction. For details
    % on these bases, please refer to Sections 4.2 and 4.3 of reference [4]
    Scr = -0.5 *k * LocLoop.ohm * (bsxfun(@besselj,nc-1,z(:))...
        - bsxfun(@besselj,nc+1,z(:))).* exp(1i*t(:)*nc);
    Spr1 = -k/(LocLoop.ohm^2-LocLoop.lambda(1)^2)*0.5 *...
        LocLoop.lambda(1) * (bsxfun(@besselj,npB1-1,z1(:)) -...
        bsxfun(@besselj,npB1+1,z1(:))) .* exp(1i*t(:)*npB1);
    Spr2 = -k/(LocLoop.ohm^2-LocLoop.lambda(2)^2)*0.5 *...
        LocLoop.lambda(2) * (bsxfun(@besselj,npB2-1,z2(:)) -...
        bsxfun(@besselj,npB2+1,z2(:))) .* exp(1i*t(:)*npB2);
    Spr3 = -k/(LocLoop.ohm^2-LocLoop.lambda(3)^2)*0.5 *...
        LocLoop.lambda(3) * (bsxfun(@besselj,npB3-1,z3(:)) -...
        bsxfun(@besselj,npB3+1,z3(:))) .* exp(1i*t(:)*npB3);
    Spr = cat(2,Spr1,Spr2,Spr3);
    
    % Computing the flux shape functions, in the Th direction. For details
    % on these bases, please refer to Sections 4.2 and 4.3 of reference [4]
    Scth = -0.5 *k * LocLoop.ohm * 1i * (bsxfun(@besselj,nc-1,z(:)) + ...
        bsxfun(@besselj,nc+1,z(:))).* exp(1i*t(:)*nc);
    Spth1 = -k/(LocLoop.ohm^2-LocLoop.lambda(1)^2)*0.5 * ...
        LocLoop.lambda(1) * 1i * (bsxfun(@besselj,npB1-1,z1(:)) +...
        bsxfun(@besselj,npB1+1,z1(:))) .* exp(1i*t(:)*npB1);
    Spth2 = -k/(LocLoop.ohm^2-LocLoop.lambda(2)^2)*0.5 * ...
        LocLoop.lambda(2) * 1i * (bsxfun(@besselj,npB2-1,z2(:)) +...
        bsxfun(@besselj,npB2+1,z2(:))) .* exp(1i*t(:)*npB2);
    Spth3 = -k/(LocLoop.ohm^2-LocLoop.lambda(3)^2)*0.5 * ...
        LocLoop.lambda(3) * 1i * (bsxfun(@besselj,npB3-1,z3(:)) +...
        bsxfun(@besselj,npB3+1,z3(:))) .* exp(1i*t(:)*npB3);
    Spth = cat(2,Spth1,Spth2,Spth3);
    
    % Computing the flux values qr and qth in the grid points, in the polar 
    % referential.
    qr = Scr*XC + Spr*XP;
    qth = Scth*XC + Spth*XP;
    
    % Transforming the flux field qr and qth from the polar to the Cartesian
    % referential. All pages in Th are equal, so the first one is selected 
    % to compute the  normal cosines.
    Loops.Qx(ii) = {zeros(NGP)};
    Loops.Qy(ii) = {zeros(NGP)};
    Qx = cos(t(:)).*qr - sin(t(:)).*qth;
    Qy = sin(t(:)).*qr + cos(t(:)).*qth;
    Qx = reshape(Qx,[NGP,NGP]);
    Qy = reshape(Qy,[NGP,NGP]);
    
    % Writing the fields in the TecPlot compatible file, if requested by
    % the user. The results are written as a matrix of NEL*(NGP+1)^2 lines
    % (where NEL is the total number of finite elements and (NGP+1)^2 the
    % total number of result points per element), and 8 columns. For each
    % result point, the columns list the x and y coordinates, in the global
    % referential, followed by the values of the temperature and flux
    % fields (complex numbers).
    if ~isempty(DirName)
        for jj = 1:NGP
            for kk = 1:NGP
                fprintf(FileU,'%0.6e %0.6e %0.6e %0.6e %0.6e %0.6e %0.6e %0.6e \n',...
                    globalx(jj,kk), globaly(jj,kk), real(T(jj,kk)), ...
                    imag(T(jj,kk)), real(Qx(jj,kk)), imag(Qx(jj,kk)), ...
                    real(Qy(jj,kk)), imag(Qy(jj,kk)));
            end
        end
    end
    
    % Storing the temperature and heat flux values in the Loops data
    % structure.
    Loops.GT0(ii) = {T};
    Loops.GV0(ii) = {V};
    Loops.Qx(ii) = {Qx};
    Loops.Qy(ii) = {Qy};
    
end
end
