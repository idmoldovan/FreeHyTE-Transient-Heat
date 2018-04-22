function varargout = TrHeatStructDef(varargin)
% TRHEATSTRUCTDEF MATLAB code for TrHeatStructDef.fig
%      TRHEATSTRUCTDEF, by itself, creates a new TRHEATSTRUCTDEF or raises the existing
%      singleton*.
%
%      H = TRHEATSTRUCTDEF returns the handle to a new TRHEATSTRUCTDEF or the handle to
%      the existing singleton*.
%
%      TRHEATSTRUCTDEF('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TRHEATSTRUCTDEF.M with the given input arguments.
%
%      TRHEATSTRUCTDEF('Property','Value',...) creates a new TRHEATSTRUCTDEF or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before TrHeatStructDef_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to TrHeatStructDef_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help TrHeatStructDef

% Last Modified by GUIDE v2.5 02-Feb-2018 17:06:46

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @TrHeatStructDef_OpeningFcn, ...
                   'gui_OutputFcn',  @TrHeatStructDef_OutputFcn, ...
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

% UIWAIT makes TrHeatStructDef wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = TrHeatStructDef_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes just before TrHeatStructDef is made visible.
function TrHeatStructDef_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to TrHeatStructDef (see VARARGIN)

%%
% If no varargin was defined, it is set to zero and saved to handles
if isempty(varargin)
    handles.varargin = 0;
else
    handles.varargin = varargin{1};
end

%%
% Setting warnings to off. These warnings are caused by missing files and
% variables before the local |mat| files are written for the first time and
% by the possibility that the problem is purely Neumann or Dirichlet. The
% warnings are re-activated after the successful execution, at the end of
% the |main.m| function.
warning('off','MATLAB:DELETE:FileNotFound');
warning('off','MATLAB:load:variableNotFound');

%%
% Creates the |handles| structure
handles.output = hObject;

%% 
% Getting the current working folder (in R2012 may be different from the
% folder you started the app from!)
handles.WorkingFolder = pwd;

%% Annotation in the LaTeX format for the problem description
mystr = '$$k \cdot \nabla^2 T \left( x,y,t \right) + Q \left( x,y,t \right) = \rho c \, \dot{T} \left( x,y,t \right)$$'; 

% The annotation command works differently in R2012(3) and R2014(5). The
% try/catch block functions essentially as a conditional compiler.
try
    annotation(handles.uipanel2,'textbox', [0.02,0.81,0.1,0.1],'string',mystr,...
        'interpreter','latex', 'FontSize',16,'LineStyle','none',...
        'FontWeight','bold');
catch errors
    annotation(gcf,'textbox', [0.02,0.81,0.1,0.1],'string',mystr,...
        'interpreter','latex', 'FontSize',16,'LineStyle','none',...
        'FontWeight','bold');
end

%%
% If there exists a local |mat| file, it loads its key elements to fill in
% the fields with data taken from the previous run.
if exist('./TrHeatStructDef.mat','file')
    load('TrHeatStructDef','L','B','Nx','Ny','EdgesOrder','LoopsOrderC',...
        'LoopsOrderP','NumberGaussPoints','MeshOption','Theta','k','rho',...
        'c','Q','deltaT','TotalTime','lambda1','lambda2','lambda3','T0',...
        'v0','SpecOutput','DirName');
    
    %%
    % In rare situations, no values are stored for L, B, Nx, Ny. If this is
    % the case, it reads NaN and thus ChangesQ always results 1. To avoid
    % this, when NaN is red, it is substituted by zero.
    if isnan(L)
        L=0; B=0; Nx=0; Ny=0;
    end
    
    %%
    % Filling in the fields with the data from the previous iteration
    set(handles.editL,'String',sprintf('%d',L));
    set(handles.editB,'String',sprintf('%d',B));
    set(handles.editNx,'String',sprintf('%d',Nx));
    set(handles.editNy,'String',sprintf('%d',Ny));
    set(handles.editEO,'String',sprintf('%d',EdgesOrder));
    set(handles.editLOC,'String',sprintf('%d',LoopsOrderC));
    set(handles.editLOP,'String',sprintf('%d',LoopsOrderP));
    set(handles.edit_k,'String',sprintf('%g',k));
    set(handles.edit_rho,'String',sprintf('%g',rho));
    set(handles.edit_c,'String',sprintf('%g',c));
    set(handles.editNGP,'String',sprintf('%d',NumberGaussPoints));
    set(handles.edittheta,'String',sprintf('%g',Theta));
    set(handles.editdeltaT,'String',sprintf('%g',deltaT));
    set(handles.edittotalT,'String',sprintf('%g',TotalTime));
    set(handles.edit_lmd1,'String',sprintf('%g',lambda1));
    set(handles.edit_lmd2,'String',sprintf('%g',lambda2));
    set(handles.edit_lmd3,'String',sprintf('%g',lambda3));
    set(handles.edit_Q,'String',Q);
    set(handles.edit_T0,'String',T0);
    set(handles.editv0,'String',v0);
    set(handles.editSpecOut,'String',SpecOutput);
    set(handles.popupmenu_mesh,'Value',MeshOption);
     
    %%
    % If |MeshOption = 2|, that is, the mesh generation is automatic, it
    % makes no sense to edit the fields associated to the regular mesh
    % generator. They become inactive.
    if MeshOption == 1
        set(handles.editL, 'enable', 'on');
        set(handles.editB, 'enable', 'on');
        set(handles.editNx, 'enable', 'on');
        set(handles.editNy, 'enable', 'on');
    else
        set(handles.editL, 'enable', 'off');
        set(handles.editB, 'enable', 'off');
        set(handles.editNx, 'enable', 'off');
        set(handles.editNy, 'enable', 'off');
    end
end

%% 
% Returns control to the user.
% Update handles structure
guidata(hObject, handles);



%% NEXT BUTTON
% * Executes on button press in |pushbutton_next|;
% * It recovers all data provided by the user in the GUI fields;
% * If the mesh is regular, it computes the |edgesArray| list (consisting of
% the exterior edges of the srtucture);
% * If the mesh is formed by triangular elements, it launches the |pdetool|
% to generate the mesh;
% * In either case, after the mesh generation is complete, it checks for
% relevant (i.e. mesh) changes as compared to the previous |mat| file. If
% such changes exist, it deletes the |mat| file corresponding to the next
% GUIs to force their definition from scratch. If no mesh changes were
% detected, the next GUI will load its previous version;
% * If the user asked for the model to be saved, it saves the information
% regarding this GUI in the specified file and folder.
function pushbutton_next_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_next (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%%
% Reading the information provided by the user in the GUI fields.
EdgesOrder = str2double(get(handles.editEO,'String'));
LoopsOrderC = str2double(get(handles.editLOC,'String'));
LoopsOrderP = str2double(get(handles.editLOP,'String'));
NumberGaussPoints = str2double(get(handles.editNGP,'String'));
MeshOption = get(handles.popupmenu_mesh,'Value');
Theta = str2double(get(handles.edittheta,'String'));
deltaT = str2double(get(handles.editdeltaT,'String'));
TotalTime = str2double(get(handles.edittotalT,'String'));
lambda1 = str2double(get(handles.edit_lmd1,'String'));
lambda2 = str2double(get(handles.edit_lmd2,'String'));
lambda3 = str2double(get(handles.edit_lmd3,'String'));
k = str2double(get(handles.edit_k,'String'));
rho = str2double(get(handles.edit_rho,'String'));
c = str2double(get(handles.edit_c,'String'));
Q = get(handles.edit_Q,'String');
T0 = get(handles.edit_T0,'String');
v0 = get(handles.editv0,'String');
SpecOutput = str2double(get(handles.editSpecOut,'String'));

%%
% If the user asked for the model to be saved, it stores the path and the
% file name.
if isfield(handles,'DirName')
    DirName = handles.DirName;
    FileName = handles.FileName;
else
    DirName = '';
    FileName = '';
end


%% 
% *Procedure for the regular rectangular mesh*
if MeshOption == 1 
    %%
    % Reading mesh information
    L = str2double(get(handles.editL,'String'));
    B = str2double(get(handles.editB,'String'));
    Nx = str2double(get(handles.editNx,'String'));
    Ny = str2double(get(handles.editNy,'String'));
    %%
    % Computing the |edgesArray| list
    L1 = linspace(1,Ny,Ny);
    L2 = linspace((Nx*Ny)+1,((Nx*Ny)+1)+(Ny-1),Ny);
    L3 = linspace(((Nx*Ny)+1)+Ny,((Nx*Ny)+1)+Ny+((Ny+1)*(Nx-1)),Nx);
    L4 = linspace(((Nx*Ny)+1)+2*Ny,((Nx*Ny)+1)+Ny+((Ny+1)*(Nx-1))+Ny,Nx);
    edgesArray = sort(cat(1,L1',L2',L3',L4'));
    %%
    % |p| and |t| variables are allocated dummy values to avoid errors
    % related to the comparison between |mat| files corresponding to
    % regular (rectangular) and triangular meshes.
    p = 0; t = 0;
    %%
    % Check if critical changes were made in the current GUI session. Note
    % that critical changes only refer to the mesh, not necessarily to the
    % geometry of the structure.
    if exist('./TrHeatStructDef.mat','file')
        save('ToDetectChanges','MeshOption','Nx','Ny','p','t'); % temporary file to detect critical changes
        % Fields whose change triggers the change of the mesh
        CriticalFields = {'MeshOption','Nx','Ny','p','t'};
        % Checks the old |TrHeatStructDef.m| and the new |ToDetectChanges| file
        % for changes in the critical fields. ChangesQ = 1 if there are 
        % changes, 0 otherwise.
        ChangesQ = any(~isfield(comp_struct(load('TrHeatStructDef.mat'),...
            load('ToDetectChanges.mat')),CriticalFields));
        % deleting the auxiliary file
        delete('ToDetectChanges.mat');
    else
        % If there is no |TrHeatStructDef| file in the first place, it sets
        % ChangesQ to 1.
        ChangesQ = 1;
    end    
    
    %%
    % Saving the workspace to the local |mat| file. If requested by the 
    % user, the GUI information is also saved to the file specified by the
    % user.   
    save('TrHeatStructDef','L','B','Nx','Ny','EdgesOrder','LoopsOrderC',...
        'LoopsOrderP','NumberGaussPoints','MeshOption','Theta','deltaT',...
        'TotalTime','lambda1','lambda2','lambda3','k','rho','c','Q',...
        'T0','v0','SpecOutput','DirName','FileName','edgesArray',...
        'p','t','ChangesQ');
        
    % If the folder where it looks for pre-load files is different from the
    % current folder, it creates a copy of TrHeatStructDef in the former
    if  ~strcmp(pwd,handles.WorkingFolder)
        save(fullfile(handles.WorkingFolder,'TrHeatStructDef'),...
        'L','B','Nx','Ny','EdgesOrder','LoopsOrderC',...
        'LoopsOrderP','NumberGaussPoints','MeshOption','Theta','deltaT',...
        'TotalTime','lambda1','lambda2','lambda3','k','rho','c','Q',...
        'T0','v0','SpecOutput','DirName','FileName','edgesArray',...
        'p','t','ChangesQ');
    end
    
    % If a savefile is requested, it creates it and stores the run data
    if ~isempty(DirName)
        save(fullfile(DirName,FileName),...
            'L','B','Nx','Ny','EdgesOrder','LoopsOrderC',...
        'LoopsOrderP','NumberGaussPoints','MeshOption','Theta','deltaT',...
        'TotalTime','lambda1','lambda2','lambda3','k','rho','c','Q',...
        'T0','v0','SpecOutput','DirName','FileName','edgesArray',...
        'p','t','ChangesQ');
    end
    
    %%
    % Preparing to exit.
    close(handles.figure1);
    
    %%
    % If there are relevant changes, it physically deletes the |mat| files
    % associated to the following GUIs.
    if ChangesQ
        delete('TrHeatRegBC1.mat','TrHeatRegBC2.mat');
    end
    TrHeatRegBC1;

    %%
    % *Procedure for the automatic triangular mesh*
else 
    %%
    % Reading mesh information. This data is useless for the automatic mesh
    % generation. It is only stored to fill in the corresponding fields in
    % a future run.
    L = str2double(get(handles.editL,'String'));
    B = str2double(get(handles.editB,'String'));
    Nx = str2double(get(handles.editNx,'String'));
    Ny = str2double(get(handles.editNy,'String'));
    %%
    % |edgesArray| variable is allocated a dummy value to avoid errors
    % related to the comparison between |mat| files corresponding to
    % regular (rectangular) and triangular meshes.
    edgesArray = 0; 
    
    %%
    % Launching |pdetool| to define the mesh
    pdewindow = pdeinit;
    %% 
    % Stopping the execution until the |pdetool| window is closed 
    waitfor(pdewindow);
    %%
    % This is a basic check to confirm that the user saved the mesh 
    % information. A new mesh need not be created if the corresponding
    % information already exist in |TrHeatStructDef.mat|, so the program checks
    % if a new mesh was defined or if |TrHeatStructDef.mat| exists. Of course,
    % the check fails if |TrHeatStructDef.mat| exists, but does not contain
    % relevant mesh information (the program exists with an error).
    BaseVars = evalin('base','whos'); % collects variable info from base
    % if the mesh variables are not found in base and an old definition
    % does not exist in TrHeatStructDef.mat
    if (~ismember('p',[BaseVars(:).name]) || ~ismember('t',[BaseVars(:).name])) && ...
            (~exist('./TrHeatStructDef.mat','file'))
        errordlg('No mesh information was exported or variable names were changed. Please press "Next" again and export the mesh info as p, e and t.','Invalid input','modal');
        uicontrol(hObject);
        return;
    end
    
    %%
    % if the mesh variables are found in base, it loads them...
    if ismember('p',[BaseVars(:).name]) && ismember('t',[BaseVars(:).name])
        p = evalin('base','p');
        t = evalin('base','t');
        %%
        % ... otherwise, it loads them from the TrHeatStructDef.mat
    else
        load('TrHeatStructDef','p','t');
    end
    
    %%
    % Check if critical changes were made in the current GUI session. Note
    % that critical changes only refer to the mesh, not necessarily to the
    % geometry of the structure.
    if exist('./TrHeatStructDef.mat','file')
        save('ToDetectChanges','MeshOption','Nx','Ny','p','t'); % temporary file to detect critical changes
        % Fields whose change triggers the change of the mesh
        CriticalFields = {'MeshOption','Nx','Ny','p','t'};
        % Checks the old TrHeatStructDef and the new ToDetectChanges for changes in
        % the critical fields. ChangesQ = 1 if there are changes, 0 otherwise
        ChangesQ = any(~isfield(comp_struct(load('TrHeatStructDef.mat'),...
            load('ToDetectChanges.mat')),CriticalFields));
        % deleting the auxiliary file
        delete('ToDetectChanges.mat');
    else
        % If there is no |TrHeatStructDef| file in the first place, it sets
        % ChangesQ to 1.
        ChangesQ = 1;
    end
    
    %%
    % Saving the workspace to the local |mat| file. If requested by the 
    % user, the GUI information is also saved to the file specified by the
    % user.
    save('TrHeatStructDef','L','B','Nx','Ny','EdgesOrder','LoopsOrderC',...
        'LoopsOrderP','NumberGaussPoints','MeshOption','Theta','deltaT',...
        'TotalTime','lambda1','lambda2','lambda3','k','rho','c','Q',...
        'T0','v0','SpecOutput','DirName','FileName','edgesArray',...
        'p','t','ChangesQ');
    
    % If the folder where it looks for pre-load files is different from the
    % current folder, it creates a copy of TrHeatStructDef in the former
    if  ~strcmp(pwd,handles.WorkingFolder)
        save(fullfile(handles.WorkingFolder,'TrHeatStructDef'),...
        'L','B','Nx','Ny','EdgesOrder','LoopsOrderC',...
        'LoopsOrderP','NumberGaussPoints','MeshOption','Theta','deltaT',...
        'TotalTime','lambda1','lambda2','lambda3','k','rho','c','Q',...
        'T0','v0','SpecOutput','DirName','FileName','edgesArray',...
        'p','t','ChangesQ');
    end   
    
    % If requested by the user, saves to file
    if ~isempty(DirName)
        save(fullfile(DirName,FileName),...
        'L','B','Nx','Ny','EdgesOrder','LoopsOrderC',...
        'LoopsOrderP','NumberGaussPoints','MeshOption','Theta','deltaT',...
        'TotalTime','lambda1','lambda2','lambda3','k','rho','c','Q',...
        'T0','v0','SpecOutput','DirName','FileName','edgesArray',...
        'p','t','ChangesQ');
    end    
    
    %%
    % Preparing to exit.   
    close(handles.figure1);
    
    %%
    % If there are relevant changes, it physically deletes the |mat| files
    % associated to the following GUIs.
    if ChangesQ
        delete('TrHeatTriBC1.mat','TrHeatTriBC2.mat');
    end
    TrHeatTriBC1;
    
end



%% RESET BUTTON
% --- Executes on button press in pushbutton_reset.
% * It deletes all entries of the edit fields and resets the pop-up menus;
% * It also deletes the save information so that the save button action is
% reversed.
function pushbutton_reset_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_reset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%%
% Removing all field information
set(handles.editL,'String','');
set(handles.editB,'String','');
set(handles.editNx,'String','');
set(handles.editNy,'String','');
set(handles.editEO,'String','');
set(handles.editLOC,'String','');
set(handles.editLOP,'String','');
set(handles.editNGP,'String','');
set(handles.edittheta,'String','');
set(handles.editdeltaT,'String','');
set(handles.edittotalT,'String','');
set(handles.edit_lmd1,'String','');
set(handles.edit_lmd2,'String','');
set(handles.edit_lmd3,'String','');
set(handles.edit_k,'String','');
set(handles.edit_rho,'String','');
set(handles.edit_c,'String','');
set(handles.edit_Q,'String','');
set(handles.edit_T0,'String','');
set(handles.editv0,'String','');
set(handles.editSpecOut,'String','');

%%
% Reseting |edit_Path|, the pop-up menus, and the properties of the mesh
% definition fields
set(handles.edit_Path,'String',...
    'THE MODEL WILL NOT BE SAVED!  Please press Save if you wish to save the model.',...
    'BackgroundColor',[1.0 0.694 0.392]);
set(handles.popupmenu_mesh,'Value',1);    

set(handles.editL, 'enable', 'on');
set(handles.editB, 'enable', 'on');
set(handles.editNx, 'enable', 'on');
set(handles.editNy, 'enable', 'on');

handles.FileName = '';
handles.DirName = '';

%%
% Updating the |handles| structure
guidata(hObject, handles);


%% SAVE BUTTON
% * Executes on button press in |Save|;
% * Reads the path and filename provided by the user to save the model;
% * Generates the |FileName| and |DirName| members in the |handles|
% structure to store this information and to make it available to other
% functions.
function Save_Callback(hObject, eventdata, handles)
% hObject    handle to Save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% 
% Reads the |handles| structure
handles=guidata(hObject);

%%
% Gets the path and filename to save the model
[FileName,DirName] = uiputfile('*.mat','Save as');

if FileName
    %%
    % Registers the path to the save file in the |edit_Path| box
    set(handles.edit_Path,'string',fullfile(DirName,FileName),...
        'BackgroundColor',[0 1 0]);
    
    %%
    % Generates the |FileName| and |DirName| members in the |handles|
    % structure to store this information and to make it available to other
    % functions.
    handles.FileName = FileName;
    handles.DirName = DirName;
    
    %%
    % If user cancelled the saving process and no save file was given
else
    handles.FileName = '';
    handles.DirName = '';
    
end

%%
% Updating the |handles| structure
guidata(hObject, handles);


%% LOAD BUTTON
% * Executes on button press in |Load|;
% * Loads all relevant variables in the file appointed by the user into the 
% workspace. If the definition of the model was not completed in the
% previous session, some of the variables cannot be loaded (the
% corresponding warning is suppressed);
% * Stores the required variables into the TrHeatStructDef |mat| file. This
% data must always be available;
% * Verifies if the data generated by the second and third GUIs are
% available. If they are, stores them into the local |mat| files, otherwise
% _deletes_ the |mat| files, to force their reinitialization;
% * Deletes the |FileName| and |DirName| variables to avoid overwriting of
% the loaded file;
% * Refreshes the interface, filling in the fields with the newly loaded
% values.
function Load_Callback(hObject, eventdata, handles)
% hObject    handle to Load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% 
% Reads the |handles| structure
handles=guidata(hObject);

%%
% Loads the |mat| file indicated by the user
[FileName,DirName] = uigetfile('*.mat','File to load');

%%
% Deletes the |FileName| and |DirName| variables to avoid overwriting;
handles.FileName = '';
handles.DirName = '';

%%
% Loading all relevant variables into the current workspace
load(fullfile(DirName,FileName),'L','B','Nx','Ny','EdgesOrder','LoopsOrderC',...
        'LoopsOrderP','NumberGaussPoints','MeshOption','Theta','deltaT',...
        'TotalTime','lambda1','lambda2','lambda3','k','rho','c','Q',...
        'T0','v0','SpecOutput','edgesArray','data','dataDir','dataNeu',...
        'dataRobin','edgesDirichlet','edgesNeumann','edgesRobin','p','t','ChangesQ');
    
%%
% Saving the local |StructDef.mat| file 
save('TrHeatStructDef','L','B','Nx','Ny','EdgesOrder','LoopsOrderC',...
        'LoopsOrderP','NumberGaussPoints','MeshOption','Theta','deltaT',...
        'TotalTime','lambda1','lambda2','lambda3','k','rho','c','Q',...
        'T0','v0','SpecOutput','DirName','FileName','edgesArray',...
        'p','t','ChangesQ');
    
% If the folder where it looks for pre-load files is different from the
% current folder, it creates a copy of StructDef in the former
if  ~strcmp(pwd,handles.WorkingFolder)
    save(fullfile(handles.WorkingFolder,'TrHeatStructDef'),'L','B','Nx',...
        'Ny','EdgesOrder','LoopsOrderC',...
        'LoopsOrderP','NumberGaussPoints','MeshOption','Theta','deltaT',...
        'TotalTime','lambda1','lambda2','lambda3','k','rho','c','Q',...
        'T0','v0','SpecOutput','DirName','FileName','edgesArray',...
        'p','t','ChangesQ');
end

%%
% If 'p' exists (is not zero), updates 'p' and 't' in the base workspace
if any(any(p~=0))
    assignin('base','p',p);
    assignin('base','t',t);
end

%%
% Depending on the kind of mesh that is selected ...
if MeshOption == 1 % regular rectangular mesh
    %%
    % ... checks if the variable created by the second GUI exists and if it
    % doesn't, it deletes the local |mat| file to force reinitialization...
    if ~exist('data','var')
        delete('TrHeatRegBC1.mat');
        %%
        % ... or it stores (and overwrites) the |mat| file if the variable
        % exists.
    else
        save('TrHeatRegBC1','data');
    end
    
    %%
    % ... checks if the variables created by the third GUI exist and if not
    % it deletes the local |mat| file to force reinitialization...
    if ~exist('dataDir','var') && ~exist('dataNeu','var') && ...
            ~exist('dataRobin','var')
        
        delete('TrHeatRegBC2.mat');
        %%
        % ... or it stores (and overwrites) the |mat| file if the variable
        % exists.
    else
        delete('TrHeatRegBC2.mat');
        dummy = 0;
        save('TrHeatRegBC2','dummy');
        if exist('dataDir','var')
            save('TrHeatRegBC2','-append','edgesDirichlet','dataDir');
        end
        if exist('dataNeu','var')
            save('TrHeatRegBC2','-append','edgesNeumann','dataNeu');
        end     
        if exist('dataRobin','var')
            save('TrHeatRegBC2','-append','edgesRobin','dataRobin');
        end          
    end
    
%%
% Same operation for the irregular mesh file.
else
    if ~exist('data','var')
        delete('TrHeatTriBC1.mat');
    else
        save('TrHeatTriBC1','data');
    end
    
    if ~exist('dataDir','var') && ~exist('dataNeu','var') && ...
            ~exist('dataRobin','var')
        delete('TrHeatTriBC2.mat');
    else
        delete('TrHeatTriBC2.mat');
        dummy = 0;
        save('TrHeatTriBC2','dummy');
        if exist('dataDir','var')
            save('TrHeatTriBC2','-append','edgesDirichlet','dataDir');
        end
        if exist('dataNeu','var')
            save('TrHeatTriBC2','-append','edgesNeumann','dataNeu');
        end
        if exist('dataRobin','var')
            save('TrHeatTriBC2','-append','edgesRobin','dataRobin');
        end
    end
end
 
set(handles.editL,'String',sprintf('%d',L));
set(handles.editB,'String',sprintf('%d',B));
set(handles.editNx,'String',sprintf('%d',Nx));
set(handles.editNy,'String',sprintf('%d',Ny));
set(handles.editEO,'String',sprintf('%d',EdgesOrder));
set(handles.editLOC,'String',sprintf('%d',LoopsOrderC));
set(handles.editLOP,'String',sprintf('%d',LoopsOrderP));
set(handles.editNGP,'String',sprintf('%d',NumberGaussPoints));
set(handles.edittheta,'String',sprintf('%d',Theta));
set(handles.editdeltaT,'String',sprintf('%d',deltaT));
set(handles.edittotalT,'String',sprintf('%d',TotalTime));
set(handles.edit_lmd1,'String',sprintf('%d',lambda1));
set(handles.edit_lmd2,'String',sprintf('%d',lambda2));
set(handles.edit_lmd3,'String',sprintf('%d',lambda3));
set(handles.edit_k,'String',sprintf('%d',k));
set(handles.edit_rho,'String',sprintf('%d',rho));
set(handles.edit_c,'String',sprintf('%d',c));
set(handles.edit_Q,'String',Q);
set(handles.edit_T0,'String',T0);
set(handles.editv0,'String',v0);
set(handles.editSpecOut,'String',sprintf('%d',SpecOutput));
set(handles.popupmenu_mesh,'Value',MeshOption);
set(handles.edit_Path,'String',...
    'THE MODEL WILL NOT BE SAVED!  Please press Save if you wish to save the model.',...
    'BackgroundColor',[1.0 0.694 0.392]);

%%
% ... and activates or inactivates the regular mesh fields according to the
% type of mesh in the loaded model.
if MeshOption == 1
    set(handles.editL, 'enable', 'on');
    set(handles.editB, 'enable', 'on');
    set(handles.editNx, 'enable', 'on');
    set(handles.editNy, 'enable', 'on');
else
    set(handles.editL, 'enable', 'off');
    set(handles.editB, 'enable', 'off');
    set(handles.editNx, 'enable', 'off');
    set(handles.editNy, 'enable', 'off');
end

%%
% Updates the |handles| structure
guidata(hObject, handles);


%% CLEAR BUTTON
% * Executes on button press in |Clear|;
% * Deletes the path and name of the save file and reinitializes the
% |edit_Path| field.
function Clear_Callback(hObject, eventdata, handles)
% hObject    handle to Clear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% 
% Reads the |handles| structure
handles=guidata(hObject);

%%
% Deletes the |FileName| and |DirName| variables 
handles.FileName = '';
handles.DirName = '';
set(handles.edit_Path,'string',...
    'THE MODEL WILL NOT BE SAVED!  Please press Save if you wish to save the model.',...
    'BackgroundColor',[1.0 0.694 0.392]);

%%
% Updates the |handles| structure
guidata(hObject, handles);














%% INSTRUCTIONS FOR THE CREATION OF THE GUI ENTITIES
% The following functions control the creation, behaviour and occasional
% checks of the entities that populate the GUI. There shouldn't be any need
% to change them.

function editL_Callback(hObject, eventdata, handles)
% hObject    handle to editL (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editL as text
%        str2double(get(hObject,'String')) returns contents of editL as a double
L = str2double(get(hObject,'String'));
if isnan(L) || ~isreal(L) || isequal(L,0) || L < 0
    set(hObject,'String','');
    errordlg('You must enter a real positive value','Invalid input','modal');
    uicontrol(hObject);
    return;
end

% --- Executes during object creation, after setting all properties.
function editL_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editL (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editB_Callback(hObject, eventdata, handles)
% hObject    handle to editB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editB as text
%        str2double(get(hObject,'String')) returns contents of editB as a double
B = str2double(get(hObject,'String'));
if isnan(B) || ~isreal(B) || isequal(B,0) || B < 0
    set(hObject,'String','');
    errordlg('You must enter a real positive value','Invalid input','modal');
    uicontrol(hObject);
    return;
end

% --- Executes during object creation, after setting all properties.
function editB_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editNx_Callback(hObject, eventdata, handles)
% hObject    handle to editNx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editNx as text
%        str2double(get(hObject,'String')) returns contents of editNx as a double
Nx = str2double(get(hObject,'String'));
if isnan(Nx) || ~isreal(Nx) || logical(abs(round(Nx)-Nx)<eps)==0 || isequal(Nx,0) || Nx < 0
    set(hObject,'String','');
    errordlg('You must enter a positive integer value','Invalid input','modal');
    uicontrol(hObject);
    return;
end

% --- Executes during object creation, after setting all properties.
function editNx_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editNx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editNy_Callback(hObject, eventdata, handles)
% hObject    handle to editNy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editNy as text
%        str2double(get(hObject,'String')) returns contents of editNy as a double
Ny = str2double(get(hObject,'String'));
if isnan(Ny) || ~isreal(Ny) || logical(abs(round(Ny)-Ny)<eps)==0 || isequal(Ny,0) || Ny < 0
    set(hObject,'String','');
    errordlg('You must enter a positive integer value','Invalid input','modal');
    uicontrol(hObject);
    return;
end

% --- Executes during object creation, after setting all properties.
function editNy_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editNy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editEO_Callback(hObject, eventdata, handles)
% hObject    handle to editEO (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editEO as text
%        str2double(get(hObject,'String')) returns contents of editEO as a double
EdgesOrder = str2double(get(handles.editEO,'String'));
if isnan(EdgesOrder) || ~isreal(EdgesOrder) ||... 
        logical(abs(round(EdgesOrder)-EdgesOrder)<eps)==0 || EdgesOrder < 0
    set(hObject,'String','');
    errordlg('You must enter either 0 or a natural number','Invalid input','modal');
    uicontrol(hObject);
    return;
end

% --- Executes during object creation, after setting all properties.
function editEO_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editEO (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editLOC_Callback(hObject, eventdata, handles)
% hObject    handle to editLOC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editLOC as text
%        str2double(get(hObject,'String')) returns contents of editLOC as a double
LoopsOrderC = str2double(get(handles.editLOC,'String'));
if isnan(LoopsOrderC) || ~isreal(LoopsOrderC) ||... 
        logical(abs(round(LoopsOrderC)-LoopsOrderC)<eps)==0 || LoopsOrderC < 0
    set(hObject,'String','');
    errordlg('You must enter either 0 or a natural number','Invalid input','modal');
    uicontrol(hObject);
    return;
end

% --- Executes during object creation, after setting all properties.
function editLOC_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editLOC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editLOP_Callback(hObject, eventdata, handles)
% hObject    handle to editLOP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editLOP as text
%        str2double(get(hObject,'String')) returns contents of editLOP as a double
LoopsOrderP = str2double(get(handles.editLOP,'String'));
if isnan(LoopsOrderP) || ~isreal(LoopsOrderP) ||... 
        logical(abs(round(LoopsOrderP)-LoopsOrderP)<eps)==0 || LoopsOrderP < 0
    set(hObject,'String','');
    errordlg('You must enter either 0 or a natural number','Invalid input','modal');
    uicontrol(hObject);
    return;
end

% --- Executes during object creation, after setting all properties.
function editLOP_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editLOP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editNGP_Callback(hObject, eventdata, handles)
% hObject    handle to editNGP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editNGP as text
%        str2double(get(hObject,'String')) returns contents of editNGP as a double
NumberGaussPoints = str2double(get(handles.editNGP,'String'));
if isnan(NumberGaussPoints) || ~isreal(NumberGaussPoints) || ...
        logical(abs(round(NumberGaussPoints)-NumberGaussPoints)<eps)==0 ||... 
        isequal(NumberGaussPoints,0) || NumberGaussPoints < 0
    set(hObject,'String','');
    errordlg('You must enter a positive integer value','Invalid input','modal');
    uicontrol(hObject);
    return;
end

% --- Executes during object creation, after setting all properties.
function editNGP_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editNGP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_k_Callback(hObject, eventdata, handles)
% hObject    handle to edit_k (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_k as text
%        str2double(get(hObject,'String')) returns contents of edit_k as a double
k = str2double(get(handles.edit_k,'String'));
if isnan(k) || ~isreal(k) || k <= 0
    set(hObject,'String','');
    errordlg('You must enter a positive real numerical value','Invalid input','modal');
    uicontrol(hObject);
    return;
end

% --- Executes during object creation, after setting all properties.
function edit_k_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_k (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_rho_Callback(hObject, eventdata, handles)
% hObject    handle to edit_rho (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_rho as text
%        str2double(get(hObject,'String')) returns contents of edit_rho as a double
rho = str2double(get(handles.edit_rho,'String'));
if isnan(rho) || ~isreal(rho) || rho <= 0
    set(hObject,'String','');
    errordlg('You must enter a real numerical value','Invalid input','modal');
    uicontrol(hObject);
    return;
end

% --- Executes during object creation, after setting all properties.
function edit_rho_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_rho (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_c_Callback(hObject, eventdata, handles)
% hObject    handle to edit_c (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_c as text
%        str2double(get(hObject,'String')) returns contents of edit_c as a double
c = str2double(get(handles.edit_c,'String'));
if isnan(c) || ~isreal(c) || c <= 0
    set(hObject,'String','');
    errordlg('You must enter a real numerical value','Invalid input','modal');
    uicontrol(hObject);
    return;
end

% --- Executes during object creation, after setting all properties.
function edit_c_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_c (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edittheta_Callback(hObject, eventdata, handles)
% hObject    handle to edittheta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edittheta as text
%        str2double(get(hObject,'String')) returns contents of edittheta as a double
Theta = str2double(get(handles.edittheta,'String'));
if isnan(Theta) || ~isreal(Theta) || Theta <= 0
    set(hObject,'String','');
    errordlg('You must enter a positive real number','Invalid input','modal');
    uicontrol(hObject);
    return;
end

% --- Executes during object creation, after setting all properties.
function edittheta_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edittheta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editdeltaT_Callback(hObject, eventdata, handles)
% hObject    handle to editdeltaT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editdeltaT as text
%        str2double(get(hObject,'String')) returns contents of editdeltaT as a double
deltaT = str2double(get(handles.edittheta,'String'));
if isnan(deltaT) || ~isreal(deltaT) || deltaT <= 0
    set(hObject,'String','');
    errordlg('You must enter a positive real number','Invalid input','modal')
    uicontrol(hObject)
    return
end

% --- Executes during object creation, after setting all properties.
function editdeltaT_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editdeltaT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edittotalT_Callback(hObject, eventdata, handles)
% hObject    handle to edittotalT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edittotalT as text
%        str2double(get(hObject,'String')) returns contents of edittotalT as a double
TotalTime = str2double(get(handles.edittheta,'String'));
if isnan(TotalTime) || ~isreal(TotalTime) || TotalTime <= 0
    set(hObject,'String','');
    errordlg('You must enter a positive real number','Invalid input','modal')
    uicontrol(hObject)
    return
end

% --- Executes during object creation, after setting all properties.
function edittotalT_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edittotalT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editSpecOut_Callback(hObject, eventdata, handles)
% hObject    handle to editSpecOut (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editSpecOut as text
%        str2double(get(hObject,'String')) returns contents of editSpecOut as a double
SpecOutput = str2double(get(handles.editSpecOut,'String'));
if isnan(SpecOutput) || ~isreal(SpecOutput) || ...
        logical(abs(round(SpecOutput)-SpecOutput)<eps)==0 || SpecOutput < 0
    set(hObject,'String','');
    errordlg('You must enter either 0 or a natural number','Invalid input','modal');
    uicontrol(hObject);
    return;
end

% --- Executes during object creation, after setting all properties.
function editSpecOut_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editSpecOut (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_lmd1_Callback(hObject, eventdata, handles)
% hObject    handle to edit_lmd1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_lmd1 as text
%        str2double(get(hObject,'String')) returns contents of edit_lmd1 as a double
c = str2double(get(handles.edit_lmd1,'String'));
if isnan(c) || ~isreal(c) || isequal(c,0)
    set(hObject,'String','');
    errordlg('You must enter a non-zero numerical value','Invalid input','modal');
    uicontrol(hObject);
    return;
end

% --- Executes during object creation, after setting all properties.
function edit_lmd1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_lmd1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_lmd2_Callback(hObject, eventdata, handles)
% hObject    handle to edit_lmd2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_lmd2 as text
%        str2double(get(hObject,'String')) returns contents of edit_lmd2 as a double
c = str2double(get(handles.edit_lmd2,'String'));
if isnan(c) || ~isreal(c) || isequal(c,0)
    set(hObject,'String','');
    errordlg('You must enter a non-zero numerical value','Invalid input','modal');
    uicontrol(hObject);
    return;
end

% --- Executes during object creation, after setting all properties.
function edit_lmd2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_lmd2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_lmd3_Callback(hObject, eventdata, handles)
% hObject    handle to edit_lmd3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_lmd3 as text
%        str2double(get(hObject,'String')) returns contents of edit_lmd3 as a double
c = str2double(get(handles.edit_lmd3,'String'));
if isnan(c) || ~isreal(c) || isequal(c,0)
    set(hObject,'String','');
    errordlg('You must enter a non-zero numerical value','Invalid input','modal');
    uicontrol(hObject);
    return;
end

% --- Executes during object creation, after setting all properties.
function edit_lmd3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_lmd3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_Path_Callback(hObject, eventdata, handles)
% hObject    handle to edit_Path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_Path as text
%        str2double(get(hObject,'String')) returns contents of edit_Path as a double

% --- Executes during object creation, after setting all properties.
function edit_Path_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_Path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in popupmenu_mesh.
function popupmenu_mesh_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_mesh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_mesh contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_mesh
if get(handles.popupmenu_mesh,'Value') == 1
    set(handles.editL, 'enable', 'on');
    set(handles.editB, 'enable', 'on');
    set(handles.editNx, 'enable', 'on');
    set(handles.editNy, 'enable', 'on');
else
   set(handles.editL, 'enable', 'off');
   set(handles.editB, 'enable', 'off');
   set(handles.editNx, 'enable', 'off');
   set(handles.editNy, 'enable', 'off');
end

% --- Executes during object creation, after setting all properties.
function popupmenu_mesh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_mesh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function text10_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% --- Executes during object creation, after setting all properties.
function text11_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% --- Executes during object creation, after setting all properties.
function text12_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% --- Executes during object creation, after setting all properties.
function text13_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% --- Executes during object creation, after setting all properties.
function text16_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% --- Executes during object creation, after setting all properties.
function text17_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% --- Executes during object creation, after setting all properties.
function text20_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% --- Executes during object creation, after setting all properties.
function text31_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text31 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


function edit_Q_Callback(hObject, eventdata, handles)
% hObject    handle to edit_Q (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_Q as text
%        str2double(get(hObject,'String')) returns contents of edit_Q as a double

% --- Executes during object creation, after setting all properties.
function edit_Q_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_Q (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editv0_Callback(hObject, eventdata, handles)
% hObject    handle to editv0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editv0 as text
%        str2double(get(hObject,'String')) returns contents of editv0 as a double

% --- Executes during object creation, after setting all properties.
function editv0_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editv0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_T0_Callback(hObject, eventdata, handles)
% hObject    handle to edit_T0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_T0 as text
%        str2double(get(hObject,'String')) returns contents of edit_T0 as a double

% --- Executes during object creation, after setting all properties.
function edit_T0_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_T0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes when figure1 is resized.
function figure1_ResizeFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
