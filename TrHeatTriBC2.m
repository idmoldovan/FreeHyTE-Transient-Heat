function varargout = TrHeatTriBC2(varargin)
% TRHEATTRIBC2 MATLAB code for TrHeatTriBC2.fig
%      TRHEATTRIBC2, by itself, creates a new TRHEATTRIBC2 or raises the existing
%      singleton*.
%
%      H = TRHEATTRIBC2 returns the handle to a new TRHEATTRIBC2 or the handle to
%      the existing singleton*.
%
%      TRHEATTRIBC2('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TRHEATTRIBC2.M with the given input arguments.
%
%      TRHEATTRIBC2('Property','Value',...) creates a new TRHEATTRIBC2 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before TrHeatTriBC2_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to TrHeatTriBC2_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to next (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help TrHeatTriBC2

% Last Modified by GUIDE v2.5 22-Nov-2017 14:37:51

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @TrHeatTriBC2_OpeningFcn, ...
    'gui_OutputFcn',  @TrHeatTriBC2_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Outputs from this function are returned to the command line.
function varargout = TrHeatTriBC2_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


%% OPENING FUNCTION
% * Executes just before |TrHeatTriBC2| is made visible;
% * Reads the mesh data from the |mat| file of the previous GUI and
% constructs the lists of Dirichlet and Neumann boundaries;
% * If the |mat| file of the current GUI exists (meaning that the 
% previous GUI was not changed), it loads the boundary information,
% otherwise it sets all boundaries to Dirichlet.
function TrHeatTriBC2_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to TrHeatTriBC2 (see VARARGIN)

%%
% Creates the |handles| structure
% Choose default command line output for |TrHeatTriBC2|
handles.output = hObject;

%%
% Loads the exterior edges and their types from the previous GUI
load('TrHeatTriBC1','edgesArray','data','edgesType');

%%
% Creates three arrays listing the exterior Dirichlet, Neumann and Robin
% edges (|edgesDirichlet|, |edgesNeumann| and |edgesRobin|)
for i=1:length(edgesArray)
    if strcmp(data{i,2},edgesType{1}) % if the edge type is Dirichlet
        handles.edgesDirichlet(i,1)=edgesArray(i); % list with the Dirichlet edges
    elseif strcmp(data{i,2},edgesType{2}) % if the edge type is Neumann
        handles.edgesNeumann(i,1)=edgesArray(i);  % list with the Neumann edges
    else    % if the edge type is Robin
        handles.edgesRobin(i,1)=edgesArray(i);  % list with the Robin edges
    end
end

%%
% *Dirichlet boundary conditions*

%%
% Generates the |dataDir| matrix to store the id of the exterior Dirichlet 
% edges and their enforced boundary conditions. 
% The |dataDir| matrix is created from from scratch or imported from the
% |TrHeatTriBC2| file, if such file exists (meaning that the previous GUI
% was left unchanged).

if isfield(handles,'edgesDirichlet')==1
    
    handles.edgesDirichlet(handles.edgesDirichlet==0)=[];
    
    %%
    % Writes all |edgesDirichlet| into the |listboxDir|
    set(handles.listboxDir,'string',handles.edgesDirichlet);
    
    %%
    % Creates the |dataDir| matrix.
    handles.dataDir=cell(length(handles.edgesDirichlet),3);
    
    %%
    % If there exists a local |mat| file, it loads it ...
    if exist('TrHeatTriBC2.mat','file') % reuse the previous data
        load('TrHeatTriBC2','dataDir');
        handles.dataDir=dataDir;
        %%
    % ... otherwise it just sets all Dirichlet boundary conditions to NaN
    else
        for i=1:length(handles.edgesDirichlet)
            handles.dataDir{i,1}=handles.edgesDirichlet(i);
            handles.dataDir{i,2}='NaN';
            handles.dataDir{i,3}='NaN';
        end
    end
    
    %%
    % Creates the Dirichlet bc table in the interface
    column1={'Bnd ID','Variation in space','Variation in time'};
    uitable('units','Normalized','Position',[0.23,0.01,0.25,0.44],...
        'ColumnWidth',{'auto' 'auto' 200},'Data',handles.dataDir,...
        'ColumnName',column1,'RowName',[]);
end


%%
% *Neumann boundary conditions*

%%
% Generates the |dataNeu| matrix to store the id of the exterior Neumann 
% edges and their enforced boundary conditions. The |dataNeu| matrix is 
% created from from scratch or imported from the |TrHeatTriBC2| file, if such 
% file exists (meaning that the previous GUI was left unchanged).

if isfield(handles,'edgesNeumann')==1

    handles.edgesNeumann(handles.edgesNeumann==0)=[];
    
    %%
    % Writes all |edgesNeumann| into the |listboxNeu|
    set(handles.listboxNeu,'string',handles.edgesNeumann);
    
    %%
    % Creates the |dataNeu| matrix.
    handles.dataNeu=cell(length(handles.edgesNeumann),3);
    
    %%
    % If there exists a local |mat| file, it loads it...
    if exist('TrHeatTriBC2.mat','file') % reuse the previous data
        load('TrHeatTriBC2','dataNeu','NeuTime');
        handles.dataNeu=dataNeu;
        %%
    % ... otherwise it just sets all Neumann boundary conditions to NaN
    else
        for i=1:length(handles.edgesNeumann)
            handles.dataNeu{i,1}=handles.edgesNeumann(i);
            handles.dataNeu{i,2}='NaN';
            handles.dataNeu{i,3}='NaN';
        end
    end
    
    %%
    % Creates the Neumann bc table in the interface
    column2={'Bnd ID','Variation in space','Variation in time'};
    uitable('units','Normalized','Position',[0.50,0.01,0.25,0.44],...
        'ColumnWidth',{'auto' 'auto' 200},'Data',handles.dataNeu,...
        'ColumnName',column2,'RowName',[]);
end


%%
% *Robin boundary conditions*

%%
% Generates the |dataRobin| matrix to store the id of the exterior Robin 
% edges and their enforced boundary conditions. The |dataRobin| matrix is 
% created from from scratch or imported from the |TrHeatTriBC2| file, if such 
% file exists (meaning that the previous GUI was left unchanged).

if isfield(handles,'edgesRobin')==1

    handles.edgesRobin(handles.edgesRobin==0)=[];
    
    %%
    % Writes all |edgesRobin| into the |listboxRobin|
    set(handles.listboxRobin,'string',handles.edgesRobin);
    
    %%
    % Creates the |dataRobin| matrix.
    handles.dataRobin=cell(length(handles.edgesRobin),3);
    
    %%
    % If there exists a local |mat| file, it loads it...
    if exist('TrHeatTriBC2.mat','file') % reuse the previous data
        load('TrHeatTriBC2','dataRobin','RobinTime');
        handles.dataRobin=dataRobin;
        %%
    % ... otherwise it just sets all Robin boundary conditions to NaN
    else
        for i=1:length(handles.edgesRobin)
            handles.dataRobin{i,1}=handles.edgesRobin(i);
            handles.dataRobin{i,2}='NaN';
            handles.dataRobin{i,3}='NaN';
        end
    end
    
    %%
    % Creates the Robin bc table in the interface
    column2={'Bnd ID','Heat coeff','Variation in time'};
    uitable('units','Normalized','Position',[0.78,0.12,0.20,0.33],...
        'ColumnWidth',{'auto' 'auto' 200},'Data',handles.dataRobin,...
        'ColumnName',column2,'RowName',[]);
end


%%
% Generates the code for drawing the mesh, along with the mesh information
% buttons

%load the mesh data
load('TrHeatTriBC1','nodes','edges_nodes','edges_loops','loops_nodes',...
    'loops_edges');

nel = size(loops_nodes,1);        % Total Number of Elements in the Mesh
nnel = size(loops_nodes,2);           % Number of nodes per Element
nnode = size(nodes,1);      % Total Number of Nodes in the Mesh
nedge = size(edges_nodes,1); % Total Number of Edges in the Mesh

% For drawing purposes
limxmin = min(nodes(:,1));
limxmax = max(nodes(:,1));
limymin =  min(nodes(:,2));
limymax =  max(nodes(:,2));

%
% Plotting the Finite Element Mesh
% Initialization of the required matrices
X = zeros(nnel,nel) ;
Y = zeros(nnel,nel) ;
% Extract X,Y coordinates for the (iel)-th element
for iel = 1:nel
    X(:,iel) = nodes(loops_nodes(iel,:),1) ;
    Y(:,iel) = nodes(loops_nodes(iel,:),2) ;
end

patch(X,Y,'w');
axis([limxmin-0.01*abs(limxmin) limxmax+0.01*abs(limxmax) limymin-0.01*abs(limymin) limymax+0.01*abs(limymax)]);
axis equal;
axis off ;

% To display Node Numbers % Element Numbers
axpos = getpixelposition(handles.axes3); % position & dimension of the axes object
% Define button's weight and height
bweight = 55;
bheight = 20;
pos = [((axpos(1)+axpos(3))/2) (axpos(2)-1.5*bheight) bweight bheight]; % align the second button with the center of the axes obj limit
ShowNodes = uicontrol('style','toggle','string','Nodes',....
    'position',[(pos(1)-2.5*bweight) pos(2) pos(3) pos(4)],'background',...
    'white');

ShowEdges = uicontrol('style','toggle','string','Edges',....
    'position',[pos(1)-0.5*bweight pos(2) pos(3) pos(4)],'background','white');

ShowElements = uicontrol('style','toggle','string','Elements',....
    'position',[(pos(1)+1.5*bweight) pos(2) pos(3) pos(4)],'background',...
    'white');

set(ShowNodes,'callback',...
    {@SHOWNODES,ShowEdges,ShowElements,nodes,edges_nodes,loops_nodes,X,Y,nnode,nedge,nel});
set(ShowEdges,'callback',...
    {@SHOWEDGES,ShowNodes,ShowElements,nodes,edges_nodes,loops_nodes,X,Y,nnode,nedge,nel});
set(ShowElements,'callback',....
    {@SHOWELEMENTS,ShowNodes,ShowEdges,nodes,edges_nodes,loops_nodes,X,Y,nnode,nedge,nel});


% Update handles structure
guidata(hObject, handles);



%% LOAD DIRICHLET FUNCTION
% * Executes on button press in |loadDirichlet|;
% * Gets the path of the data file for the temperature variation;
% * Checks the file for consistency.
function loadDirichlet_Callback(hObject, eventdata, handles)
% hObject    handle to loadDirichlet (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%%
% Reads the |handles| structure
handles=guidata(hObject);

%%
% Loads the |mat| file indicated by the user
[FileName,DirName] = uigetfile('*.*','File to load');

%%
% Writes the |FileNameDir| and |DirNameDir| variables
handles.FileNameDir = FileName;
handles.DirNameDir = DirName;

%%
% Checking if the file is readable and if it is, writing the file name to
% the |editDirT| window
try
    load(fullfile(DirName,FileName));
    set(handles.editDirT,'string',fullfile(DirName,FileName));
catch err
    errordlg(err.message,'Invalid input','modal');
    uicontrol(hObject);
end

    
%% LOAD NEUMANN FUNCTION
% * Executes on button press in |loadNeumann|;
% * Gets the path of the data file for the temperature variation;
% * Checks the file for consistency.
function loadNeumann_Callback(hObject, eventdata, handles)
% hObject    handle to loadNeumann (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%%
% Reads the |handles| structure
handles=guidata(hObject);

%%
% Loads the |mat| file indicated by the user
[FileName,DirName] = uigetfile('*.*','File to load');

%%
% Writes the |FileNameNeu| and |DirNameNeu| variables
handles.FileNameNeu = FileName;
handles.DirNameNeu = DirName;

%%
% Checking if the file is readable and if it is, writing the file name to
% the |editNeuT| window
try
    load(fullfile(DirName,FileName));
    set(handles.editNeuT,'string',fullfile(DirName,FileName));
catch err
    errordlg(err.message,'Invalid input','modal');
    uicontrol(hObject);
end


%% LOAD ROBIN FUNCTION
% * Executes on button press in |loadRobin|;
% * Gets the path of the data file for the temperature variation;
% * Checks the file for consistency.
function loadRobin_Callback(hObject, eventdata, handles)
% hObject    handle to loadRobin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%%
% Reads the |handles| structure
handles=guidata(hObject);

%%
% Loads the |mat| file indicated by the user
[FileName,DirName] = uigetfile('*.*','File to load');

%%
% Writes the |FileNameRobin| and |DirNameRobin| variables
handles.FileNameRobin = FileName;
handles.DirNameRobin = DirName;

%%
% Checking if the file is readable and if it is, writing the file name to
% the |editRobinT| window
try
    load(fullfile(DirName,FileName));
    set(handles.editRobinT,'string',fullfile(DirName,FileName));
catch err
    errordlg(err.message,'Invalid input','modal');
    uicontrol(hObject);
end


%% ASSIGN DIRICHLET FUNCTION
% * Executes on button press in |assignDirichlet|;
% * Fills in the enforced flux table for the boundaries selected in
% |listboxDir| with the strings defined in |editState| and |editDirT|.
function assignDirichlet_Callback(hObject, eventdata, handles)
% hObject    handle to assignDirichlet (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%%
% Gets the selected boundaries in |listboxDir| 
itemlist = get(handles.listboxDir,'Value'); % list of selected items in the listbox
nitems = length(itemlist); % number of selected items in the listbox

%%
% Writes the flux definition in the |dataDir| table.
for ii = 1:nitems
    crtitem = itemlist(ii);
    dataDir2 = get(handles.editState,'string');
    if ~isempty(dataDir2)
        handles.dataDir{crtitem,2}=dataDir2;
    end
    dataDir3 = get(handles.editDirT,'string');
    if ~isempty(dataDir3)
        handles.dataDir{crtitem,3}=dataDir3;
    end
end

%%
% Updates the handles structure
guidata(hObject,handles);

%%
% Redraws the table where the Dirichlet boundary conditions are listed
column1={'Bnd ID','Variation in space','Variation in time'};
uitable('units','Normalized','Position',[0.23,0.01,0.25,0.44],...
    'ColumnWidth',{'auto' 'auto' 200},'Data',handles.dataDir,...
    'ColumnName',column1,'RowName',[]);



%% ASSIGN NEUMANN FUNCTION
% * Executes on button press in |assignNeumann|;
% * Fills in the enforced flux table for the boundaries selected in
% |listboxNeu| with the strings defined in |editFlux| and |editNeuT|.
function assignNeumann_Callback(hObject, eventdata, handles)
% hObject    handle to assignNeumann (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%%
% Gets the selected boundaries in |listboxNeu| 
itemlist = get(handles.listboxNeu,'Value'); % list of selected items in the listbox
nitems = length(itemlist); % number of selected items in the listbox

%%
% Writes the flux definition in the |dataNeu| table.
for ii = 1:nitems
    crtitem = itemlist(ii);
    dataNeu2 = get(handles.editFlux,'string');
    if ~isempty(dataNeu2)
        handles.dataNeu{crtitem,2}=dataNeu2;
    end
    dataNeu3 = get(handles.editNeuT,'string');
    if ~isempty(dataNeu3)
        handles.dataNeu{crtitem,3}=dataNeu3;
    end
end

%%
% Updates the handles structure
guidata(hObject,handles);

%%
% Redraws the table where the Neumann boundary conditions are listed
column2={'Bnd ID','Variation in space','Variation in time'};
uitable('units','Normalized','Position',[0.50,0.01,0.25,0.44],...
    'ColumnWidth',{'auto' 'auto' 200},'Data',handles.dataNeu,...
    'ColumnName',column2,'RowName',[]);



%% ASSIGN ROBIN FUNCTION
% --- Executes on button press in assignRobin.
% * Fills in the enforced flux table for the boundaries selected in
% |listboxNeu| with the strings defined in |editHeatTransCoef| and |editRobinT|.
function assignRobin_Callback(hObject, eventdata, handles)
% hObject    handle to assignRobin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%%
% Gets the selected boundaries in |listboxNeu| 
itemlist = get(handles.listboxRobin,'Value'); % list of selected items in the listbox
nitems = length(itemlist); % number of selected items in the listbox

%%
% Writes the flux definition in the |dataNeu| table.
for ii = 1:nitems
    crtitem = itemlist(ii);
    dataRobin2 = get(handles.editHeatTransCoef,'string');
    if ~isempty(dataRobin2)
        handles.dataRobin{crtitem,2}=dataRobin2;
    end
    dataRobin3 = get(handles.editRobinT,'string');
    if ~isempty(dataRobin3)
        handles.dataRobin{crtitem,3}=dataRobin3;
    end    
end

%%
% Updates the handles structure
guidata(hObject,handles);

%%
% Redraws the table where the Neumann boundary conditions are listed
column2={'Bnd ID','Heat coeff','Variation in time'};
uitable('units','Normalized','Position',[0.78,0.12,0.20,0.33],...
    'ColumnWidth',{'auto' 'auto' 200},'Data',handles.dataRobin,...
    'ColumnName',column2,'RowName',[]);



%% RESET DIRICHLET/NEUMANN/ROBIN BOUNDARY CONDITION FUNCTIONS


%%
% *Reset Dirichlet*
%
% * Executes on button press in |resetDirichlet|;
% * Substitutes all previous definitions in |editState| and |editDirT| by |NaN|;
% * Redraws the state boundary condition table.
function resetDirichlet_Callback(hObject, eventdata, handles)
% hObject    handle to resetDirichlet (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%%
% Deletes the time variation field |dataDir| and |editDirT|
set(handles.editState,'String','');
set(handles.editDirT,'String','');

%%
% Substitutes all previous definitions in |dataDir| by |NaN|
for i=1:length(handles.edgesDirichlet)
    handles.dataDir{i,2}='NaN';
    handles.dataDir{i,3}='NaN';
end

%%
% Updates the handles structure
guidata(hObject,handles);

%%
% Redraws the table where the Dirichlet boundary conditions are listed
column1={'Bnd ID','Variation in space','Variation in time'};
uitable('units','Normalized','Position',[0.23,0.01,0.25,0.44],...
    'ColumnWidth',{'auto' 'auto' 200},'Data',handles.dataDir,...
    'ColumnName',column1,'RowName',[]);


%%
% *Reset Neumann*
%
% * Executes on button press in |resetNeumann|;
% * Substitutes all previous definitions in |editFlux| and |editNeuT| by |NaN|;
% * Redraws the flux boundary condition table.
function resetNeumann_Callback(hObject, eventdata, handles)
% hObject    handle to resetNeumann (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%%
% Deletes the time variation fields |dataNeu| and |editNeuT|
set(handles.editFlux,'String','');
set(handles.editNeuT,'String','');

%%
% Substitutes all previous definitions in |dataNeu| by |NaN|
for i=1:length(handles.edgesNeumann)
    handles.dataNeu{i,2}='NaN';
    handles.dataNeu{i,3}='NaN';
end

%%
% Updates the handles structure
guidata(hObject,handles);

%%
% Redraws the table where the Neumann boundary conditions are listed
column2={'Bnd ID','Variation in space','Variation in time'};
uitable('units','Normalized','Position',[0.50,0.01,0.25,0.44],...
    'ColumnWidth',{'auto' 'auto' 200},'Data',handles.dataNeu,...
    'ColumnName',column2,'RowName',[]);


%%
% *Reset Robin*
%
% * Executes on button press in |resetRobin|;
% * Substitutes all previous definitions in |editHeatTransCoef| and |editRobinT| by |NaN|;
% * Redraws the flux boundary condition table.
function resetRobin_Callback(hObject, eventdata, handles)
% hObject    handle to resetRobin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%%
% Deletes the definition fields |editHeatTransCoef| and |editRobinT|
set(handles.editHeatTransCoef,'String','');
set(handles.editRobinT,'String','');

%%
% Substitutes all previous definitions in |datRobin| by |NaN|
for i=1:length(handles.edgesRobin)
    handles.dataRobin{i,2}='NaN';
    handles.dataRobin{i,3}='NaN';
end

%%
% Updates the handles structure
guidata(hObject,handles);

%%
% Redraws the table where the Neumann boundary conditions are listed
column2={'Bnd ID','Heat coeff','Variation in time'};
uitable('units','Normalized','Position',[0.78,0.12,0.20,0.33],...
    'ColumnWidth',{'auto' 'auto' 200},'Data',handles.dataRobin,...
    'ColumnName',column2,'RowName',[]);


%% NEXT FUNCTION
% * Executes on button press in |next|;
% * Reads the boundary condition data, stores it in the local |mat| file
% and starts the checking GUI.

function next_Callback(hObject, eventdata, handles)
% hObject    handle to next (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%%
% Creating the save file
Dummy = 0;
save('TrHeatTriBC2','Dummy');

% Loading the save file data
load('TrHeatStructDef','DirName','FileName');

%%
% Recovering the user-defined data 
if isfield(handles,'edgesDirichlet')
    edgesDirichlet=handles.edgesDirichlet;
    dataDir=handles.dataDir;
    % Appending data to the local |mat| file
    save('TrHeatTriBC2','-append','edgesDirichlet','dataDir');
    
    % Saving to the save file if DirName is not empty (meaning that save
    % was requested by the user)
    if ~isempty(DirName)
        % saving the data
        save(fullfile(DirName,FileName),'-append','edgesDirichlet',...
            'dataDir');
    end
    
    % Converting strings in data arrays to matrices
    D = cellfun(@str2num,dataDir(:,2:end),'UniformOutput',0);
    % Looking for NaN in D 
    NaND = cellfun(@any,cellfun(@isnan,D,'UniformOutput',0));
    %%
    % Checking the data
    if any(NaND) % if the Dirichlet conditions on an edge are NaN
        errordlg('Dirichlet boundary conditions are set to NaN at least for one boundary.','Invalid input','modal');
        return;
    end
end
if isfield(handles,'edgesNeumann')
    edgesNeumann=handles.edgesNeumann;
    dataNeu=handles.dataNeu;
    save('TrHeatTriBC2','-append','edgesNeumann','dataNeu');
    
    % Saving to the save file if DirName is not empty (meaning that save
    % was requested by the user)
    if ~isempty(DirName)
        % saving the data
        save(fullfile(DirName,FileName),'-append','edgesNeumann',...
            'dataNeu');
    end
    
    % Converting strings in data arrays to matrices
    N = cellfun(@str2num,dataNeu(:,2:end),'UniformOutput',0);
    % Looking for NaN in N
    NaNN = cellfun(@any,cellfun(@isnan,N,'UniformOutput',0));
    %%
    % Checking the data
    if any(NaNN) % if the Neumann conditions on an edge are NaN
        errordlg('Neumann boundary conditions are set to NaN at least for one boundary.','Invalid input','modal');
        return;
    end
end
if isfield(handles,'edgesRobin')
    edgesRobin=handles.edgesRobin;
    dataRobin=handles.dataRobin;
    save('TrHeatTriBC2','-append','edgesRobin','dataRobin');
    
    % Saving to the save file if DirName is not empty (meaning that save
    % was requested by the user)
    if ~isempty(DirName)
        % saving the data
        save(fullfile(DirName,FileName),'-append','edgesRobin',...
            'dataRobin');
    end
    
    % Converting strings in data arrays to matrices
    R = cellfun(@str2num,dataRobin(:,2:end),'UniformOutput',0);
    % Looking for NaN in R
    NaNR = cellfun(@any,cellfun(@isnan,R,'UniformOutput',0));
    %%
    % Checking the data
    if any(NaNR) % if the Robin conditions on an edge are NaN
        errordlg('Heat transfer coefficient or fluid temperature are set to NaN at least for one boundary.','Invalid input','modal');
        return;
    end
end

%%
% Closes everything and launches the ckecking GUI
close(handles.figure1); % closing the GUI window

% Trying to close Visualize if it was opened
try
    close('Visualize');
catch
end

CheckStructTri;


%% PREVIOUS FUNCTION
% * Executes on button press in |previous|;
% * Just closes the current GUI and launches the previous one. All changes
% made in the current GUI are lost.

% --- Executes on button press in previous.
function previous_Callback(hObject, eventdata, handles)
% hObject    handle to previous (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close(TrHeatTriBC2);

% Trying to close Visualize if it was opened
try
    close('Visualize');
catch
end

TrHeatTriBC1;


%% ENLARGE FUNCTION
% * Executes on button press in |enlarge|;
% * Simply calls VisualizeTri
function enlarge_Callback(hObject, eventdata, handles)
% hObject    handle to enlarge (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

VisualizeTri;



%% GUI ENTITIES GENERATION CODE
% * Automatically generated code for the buttons and menus;
% * Some (not-so-sound) checks are performed on the data inserted by the
% user
% --- Executes on selection change in listboxNeu.
function listboxNeu_Callback(hObject, eventdata, handles)
% hObject    handle to listboxNeu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listboxNeu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listboxNeu

% --- Executes during object creation, after setting all properties.
function listboxNeu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listboxNeu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in listboxDir.
function listboxDir_Callback(hObject, eventdata, handles)
% hObject    handle to listboxDir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listboxDir contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listboxDir


% --- Executes during object creation, after setting all properties.
function listboxDir_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listboxDir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in listboxRobin.
function listboxRobin_Callback(hObject, eventdata, handles)
% hObject    handle to listboxRobin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listboxRobin contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listboxRobin


% --- Executes during object creation, after setting all properties.
function listboxRobin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listboxRobin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function editFlux_Callback(hObject, eventdata, handles)
% hObject    handle to editFlux (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editFlux as text
%        str2double(get(hObject,'String')) returns contents of editFlux as a double

[Flux, status] = str2num(get(hObject,'string'));
if any(isnan(Flux)) || ~status  % if the input is something else than
                                     % a vector of reals
    errordlg('Flux field must have real values','Invalid input','modal');
    uicontrol(hObject);
    return;
end


% --- Executes during object creation, after setting all properties.
function editFlux_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editFlux (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function editHeatTransCoef_Callback(hObject, eventdata, handles)
% hObject    handle to editHeatTransCoef (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editHeatTransCoef as text
%        str2double(get(hObject,'String')) returns contents of editHeatTransCoef as a double

HTC = str2double(get(handles.editHeatTransCoef,'string'));
if isnan(HTC) || ~isreal(HTC) || HTC <= 0  
    set(hObject,'String','');
    errordlg('You must enter a positive integer value','Invalid input','modal');
    uicontrol(hObject);
    return;
end

% --- Executes during object creation, after setting all properties.
function editHeatTransCoef_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editHeatTransCoef (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function editState_Callback(hObject, eventdata, handles)
% hObject    handle to editState (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editState as text
%        str2double(get(hObject,'String')) returns contents of editState as a double
[U, status] = str2num(get(hObject,'string'));
if any(isnan(U)) || ~status  % if the input is something else than
                                     % a vector of reals
    errordlg('Sate field must have real value','Invalid input','modal');
    uicontrol(hObject);
    return;
end

% --- Executes during object creation, after setting all properties.
function editState_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editState (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editDirT_Callback(hObject, eventdata, handles)
% hObject    handle to editDirT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editDirT as text
%        str2double(get(hObject,'String')) returns contents of editDirT as a double


% --- Executes during object creation, after setting all properties.
function editDirT_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editDirT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function editNeuT_Callback(hObject, eventdata, handles)
% hObject    handle to editNeuT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editNeuT as text
%        str2double(get(hObject,'String')) returns contents of editNeuT as a double


% --- Executes during object creation, after setting all properties.
function editNeuT_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editNeuT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function editRobinT_Callback(hObject, eventdata, handles)
% hObject    handle to editRobinT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editRobinT as text
%        str2double(get(hObject,'String')) returns contents of editRobinT as a double


% --- Executes during object creation, after setting all properties.
function editRobinT_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editRobinT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
