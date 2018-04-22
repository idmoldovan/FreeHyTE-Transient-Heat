function PlotFieldsTri(Nodes,Edges,Loops,time,deltaT,SpecOutput,NGP)
% PLOTFIELDSTRI plots the final temperature and flux fields in NGP*NGP
% Gauss points of each element if the plot is requested for the current
% time step.
%
% PLOTFIELDSTRI is called by MAINTRI
%
% Input:
%  * structures Edges and Loops where the final temperature and flux fields
%  are stored, the node position list Nodes, the current time, deltaT (time
%  step), SpecOutput the number time steps at which the results should be
%  stored, and the Gauss abscissas
% Output:
%  * 4 figures - three color maps of the temperature and heat flux fields,
%  and one surf plot of the temperature field.
%
% This is a plotting-only function. The values of the final temperature and
% heat flux values were computed and stored in the Loops structure in
% COMPUTEGAUSSFIELDSTRI.
%

%% Pre-flight
tstep = round(time/deltaT);

%% Mesh drawing zone
if (SpecOutput && ~rem(tstep,SpecOutput)) % if the current step number is a
    % multiple of SpecOutput, a plot is required for this time step
    
    % get the coordinates of the mesh points
    xmesh = [Nodes(Edges.nini(:),1) Nodes(Edges.nfin(:),1)];
    ymesh = [Nodes(Edges.nini(:),2) Nodes(Edges.nfin(:),2)];
    
    %% Drawing the mesh, to prepare the plotting of the fields
    Fig=figure;
    Figname = sprintf('T field, TS %d',tstep);
    set(Fig,'name',Figname,'numbertitle','off','color','w') ;
    
    %% Temperature colour map plot
    % Preparing the mesh for the temperature plot
    Tre = subplot(2,2,1);
    hold on; title('T');
    axis([0. max(max(xmesh)) 0. max(max(ymesh))]);
    % Drawing the mesh, one edge at a time
    for ii=1:length(Edges.type)
        line(xmesh(ii,:),ymesh(ii,:),'color','b');
    end
    daspect([1 1 1]);
    
    %% Temperature surf plot
    % Preparing the mesh for the temperature plot
    Tsrf = subplot(2,2,2);
    hold on; title('T (surf)');
    axis([0. max(max(xmesh)) 0. max(max(ymesh))]);
    for ii=1:length(Edges.type)
        line(xmesh(ii,:),ymesh(ii,:),'color','b');
    end
    daspect([1 1 1]);
    
    %% Heat flux plot, in X
    % Preparing the mesh for the Qx plot
    Qxre = subplot(2,2,3);
    hold on; title('Qx');
    axis([0. max(max(xmesh)) 0. max(max(ymesh))]);
    % Drawing the mesh, one edge at a time
    for ii=1:length(Edges.type)
        line(xmesh(ii,:),ymesh(ii,:),'color','b');
    end
    daspect([1 1 1]);
    
    %% Heat flux plot, in Y
    % Preparing the mesh for the Qy plot
    Qyre = subplot(2,2,4);
    hold on; title('Qy');
    axis([0. max(max(xmesh)) 0. max(max(ymesh))]);
    % Drawing the mesh, one edge at a time
    for ii=1:length(Edges.type)
        line(xmesh(ii,:),ymesh(ii,:),'color','b');
    end
    daspect([1 1 1]);
    
    %% Computation of the temperature and heat flux fields
    % Sweeping the elements to draw the colour maps
    for ii=1:length(Loops.area)
        % LocLoop is a local structure where the features of the current
        % element which are directly useful for the plotting of the output
        % fields are stored.
        LocLoop = struct('id',ii,'nodes',Loops.nodes(ii,:), ...
            'edges',Loops.edges(ii,:),'center',Loops.center(ii,:),...
            'T',cell2mat(Loops.GT0(ii)),'V',cell2mat(Loops.GV0(ii)),...
            'Qx',cell2mat(Loops.Qx(ii)),'Qy',cell2mat(Loops.Qy(ii)));
        
        % Getting coordinates of the nodes of the element (global)
        LocNodes = Nodes(LocLoop.nodes(:),:);
        
        % Getting the Gauus points of the element
        [GlobalX,GlobalY,~,~]=triquad(NGP,LocNodes);
        
        % PLOTTING AREA
        subplot(Tre); contourf(GlobalX, GlobalY, real(LocLoop.T),...
            20,'edgecolor','none');
        
        subplot(Tsrf); surf(GlobalX, GlobalY, real(LocLoop.T)); axis square;
        
        subplot(Qxre); contourf(GlobalX, GlobalY, real(LocLoop.Qx),...
            20,'edgecolor','none');
        
        subplot(Qyre); contourf(GlobalX, GlobalY, real(LocLoop.Qy),...
            20,'edgecolor','none');
        
    end
    
    % Plotting the legends
    subplot(Tre); caxis(caxis); colorbar; colormap(jet);
    subplot(Qxre); caxis(caxis); colorbar; colormap(jet);
    subplot(Qyre); caxis(caxis); colorbar; colormap(jet);
    
end
