function [Loops, TotalTime, deltaT, Ini, InternHeat, Theta] = TimeGen(Loops) 
% TIMEGEN processes the time discretization data.

% TIMEGEN is called by MAIN***. It reads the time discretizaton data 
% inserted in the GUIs, computes the algorithmic wave number ohm, and
% returns handles to the functions that compute the initial fields and heat
% source.
%
% INPUT DATA (read from the TrHeatStructDef.mat file created by the GUI1):
% * Theta: the generalized mid-point algorithmic parameter;
% * deltaT: the time step;
% * TotalTime: the total time of the analysis.
% * Loops structure (from MAIN***).
%
% OUTPUT DATA (returned to MAIN***):
% * the data read from the GUI, described above;
% * Ini: handle to the function GENINI0 that computes the initial
% conditions in a given matrix of points;
% * InternHeat: handle to the function GENINTERNHEAT that computes the 
% heat source in a given matrix of points;
% * Loops structure, updated with the algorithmic wave number.
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
% Details on the time discretization process for parabolic problems can be
% found in Reference [2], Section 5.2.1.

%% Loads time discretization information from TrHeatStructDef.mat
load('TrHeatStructDef','Theta','deltaT','TotalTime');

%% Computes the wave number ohm for all elements
for ii=1:length(Loops.area)
    Loops.ohm(ii) = sqrt(-1/(Loops.material(ii,4)*Theta*deltaT));
end

%% Get handles to the initial conditions and heat source functions
Ini = @GenIni0;
InternHeat = @GenInternHeat;

end

function [T0, V0] = GenIni0(globalx,globaly)
% GENINI0 is a local function that computes the initial conditions in
% points (globalx,globaly). It receives globalx and globaly as input and
% returns two 2D matrices with the T0 and V0 values in the input points.

%% working out the input data

% T0 and v0 are strings loaded from TrHeatStructDef.mat. Each contains an
% analytic expression in 'x' and 'y' that defines the initial conditions
% in the global coordinates.
load('TrHeatStructDef','T0','v0');

% The following lines convert 'x' to 'globalx', 'y' to 'globaly', and the
% generic operations to term-by-term operations.
T0str = strrep(strrep(strrep(strrep(strrep(T0,'x','globalx'),'y','globaly'),'*','.*'),'/','./'),'^','.^');
V0str = strrep(strrep(strrep(strrep(strrep(v0,'x','globalx'),'y','globaly'),'*','.*'),'/','./'),'^','.^');
% The problem is that 'x' is found in the commonly used function 'exp',
% where the substitution of 'x' by 'globalx' must be undone.
T0str = strrep(T0str,'eglobalxp','exp');
V0str = strrep(V0str,'eglobalxp','exp');
% Evaluating the analytic expressions
T0 = eval(T0str);
V0 = eval(V0str);

end

function bb = GenInternHeat(globalx,globaly,crttime)
% GENINITERNHEAT is a local function that computes the heat source in
% points (globalx,globaly), at instant crttime. It receives globalx and
% globaly as input and returns a 2D matrix with the internal heat values in
% the requested points. 

%% working out the input data

% Q is a string loaded from TrHeatStructDef.mat. It contains an analytic 
% expression in 'x' and 'y' that defines the heat source in the global
% coordinates. 
load('TrHeatStructDef','Q');
% The following line converts 'x' to 'globalx', 'y' to 'globaly', and the
% generic operations to term-by-term operations.
bbstr = strrep(strrep(strrep(strrep(strrep(strrep(Q,'x','globalx'),'y','globaly'),'t','crttime'),'*','.*'),'/','./'),'^','.^');
% The problem is that 'x' is found in the commonly used function 'exp',
% where the substitution of 'x' by 'globalx' must be undone.
bbstr = strrep(bbstr,'eglobalxp','exp');
% Evaluating the analytic expression
bb = eval(bbstr);

end
