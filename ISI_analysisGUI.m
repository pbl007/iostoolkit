function varargout = ISI_analysisGUI(varargin)
% ISI_analysisGUI M-file for ISI_analysisGUI.fig
%      ISI_analysisGUI, by itself, creates a new ISI_analysisGUI or raises
%      the existing
%      singleton*.
%
%      H = ISI_analysisGUI returns the handle to a new ISI_analysisGUI or the handle to
%      the existing singleton*.
%
%      ISI_analysisGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ISI_analysisGUI.M with the given input arguments.
%
%      ISI_analysisGUI('Property','Value',...) creates a new ISI_analysisGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ISI_analysisGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property
%      application
%      stop.  All inputs are passed to ISI_analysisGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%

%VERSIONS:
% 2011.10.31    mjp     initial version


% Last Modified by GUIDE v2.5 16-Mar-2012 12:44:46

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @ISI_analysisGUI_OpeningFcn, ...
    'gui_OutputFcn',  @ISI_analysisGUI_OutputFcn, ...
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


% --- Executes just before ISI_analysisGUI is made visible.
function ISI_analysisGUI_OpeningFcn(hObject, eventdata, handles, varargin) %#ok
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ISI_analysisGUI (see VARARGIN)

handles.pathstr = pwd;     % initial file directory
set(handles.path,'string',handles.pathstr);

set(hObject,'Tag', 'ISIanalysisGUI_fig')

% Choose default command line output for ISI_analysisGUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ISI_analysisGUI wait for user response (see UIRESUME)
% uiwait(handles.ISIanalysisGUI_fig);


return

% --- Outputs from this function are returned to the command line.
function varargout = ISI_analysisGUI_OutputFcn(hObject, eventdata, handles) %#ok
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Get default command line output from handles structure
varargout{1} = handles.output;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% file/path info


function path_Callback(hObject, eventdata, handles) %#ok
curpath=get(handles.path,'string');
if ~isempty(curpath) && ~exist(curpath,'dir') %directory doesnt exist
    warndlg(sprintf('The folder: \n%s\ndoes not exist.',curpath));
    set(handles.path,'string',pwd);
    handles.pathstr=pwd;
else
    handles.pathstr=curpath;
end
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function path_CreateFcn(hObject, eventdata, handles)  %#ok
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in getPath.
function getPath_Callback(hObject, eventdata, handles)  %#ok
oldpath=handles.pathstr;
if isempty(oldpath), oldpath=pwd;  end
curpath=uigetdir(oldpath);
if curpath==0, curpath=pwd; end %cancelled out
set(handles.path,'string',curpath);
handles.pathstr=curpath;
guidata(hObject, handles);


function data_filename_Callback(hObject, eventdata, handles)  %#ok
datafile=get(handles.data_filename,'string');
if ~exist(fullfile(handles.pathstr,datafile),'file') %file doesnt exist
    warndlg(sprintf('The file: \n%s\ndoes not exist.',datafile));
    set(handles.data_filename,'string','');
else
    name=regexp(handles.params.data_filename,'(.+)_\w{6}.dat','tokens');
    set(handles.whiskername,'string',name{1});
end
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function data_filename_CreateFcn(hObject, eventdata, handles)  %#ok
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in getdatafile.
function getdatafile_Callback(hObject, eventdata, handles) %#ok

if isempty(eventdata)
    curdir = pwd;
    cd(handles.pathstr); %load to path directory
    oldfile=get(handles.data_filename,'string');
    [datafile,path2file]=uigetfile('*.dat');
    
    if datafile==0, datafile=oldfile; end %cancelled out
    cd(curdir);
    set(handles.data_filename,'string',datafile);
    %update also path to selected file
    set(handles.path,'string',path2file);
    handles.pathstr=path2file;
    
else
    % use passed file name
    datafile = eventdata;
    set(handles.data_filename,'string',datafile);
end

name=regexp(datafile,'(.+)_\w{6}.dat','tokens');

if ~isempty(name)
    set(handles.whiskername,'string',name{1});
else
    %     warning('ISI_analysisGUI:getdatafile','Couldn''t isolate whisker name, please enter it manually.');
    warndlg('Couldn''t isolate whisker name, please enter it manually.','ISI_analysisGUI:getdatafile');
    
    set(handles.whiskername,'string','');
end

guidata(hObject, handles);

hGUI = findobj('Tag', 'ISIanalysisGUI_fig');
set(hGUI, 'UserData', [])

return


function vessel_filename_Callback(hObject, eventdata, handles)  %#ok
vesselfile=get(handles.vessel_filename,'string');
if ~exist(fullfile(handles.pathstr,vesselfile),'file') %file doesnt exist
    warndlg(sprintf('The file: \n%s\ndoes not exist.',vesselfile));
    set(handles.vessel_filename,'string','');
end
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function vessel_filename_CreateFcn(hObject, eventdata, handles)  %#ok
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in getvessel.
function getvessel_Callback(hObject, eventdata, handles) %#ok
curdir=pwd;
cd(handles.pathstr); %load to path directory
oldfile=get(handles.vessel_filename,'string');
vesselfile=uigetfile({'*.png'});
if vesselfile==0, vesselfile=oldfile; end %cancelled out
cd(curdir);
set(handles.vessel_filename,'string',vesselfile);
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function whiskername_CreateFcn(hObject, eventdata, handles) %#ok
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function whiskername_Callback(hObject, eventdata, handles) %#ok
%don't do anything, just scrape raw string, since it's only for display


function nolight_filename_Callback(hObject, eventdata, handles) %#ok
nolightfile=get(handles.nolight_filename,'string');
if ~exist(fullfile(handles.pathstr,nolightfile),'file') %file doesnt exist
    warndlg(sprintf('The file: \n%s\ndoes not exist.',nolightfile));
    set(handles.nolight_filename,'string','');
end
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function nolight_filename_CreateFcn(hObject, eventdata, handles) %#ok
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in getnolight.
function getnolight_Callback(hObject, eventdata, handles) %#ok
nolightfile=uigetfile('*.dat');
set(handles.nolight_filename,'string',nolightfile);
guidata(hObject, handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% analysis parameters
function preStimDurSec_Callback(hObject, eventdata, handles) %#ok

% --- Executes during object creation, after setting all properties.
function preStimDurSec_CreateFcn(hObject, eventdata, handles) %#ok
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function stimDurSec_Callback(hObject, eventdata, handles) %#ok

% --- Executes during object creation, after setting all properties.
function stimDurSec_CreateFcn(hObject, eventdata, handles) %#ok
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function contour_min_Callback(hObject, eventdata, handles)  %#ok

% --- Executes during object creation, after setting all properties.
function contour_min_CreateFcn(hObject, eventdata, handles)  %#ok
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function contour_max_Callback(hObject, eventdata, handles)  %#ok

% --- Executes during object creation, after setting all properties.
function contour_max_CreateFcn(hObject, eventdata, handles)  %#ok
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function contour_step_Callback(hObject, eventdata, handles)  %#ok

% --- Executes during object creation, after setting all properties.
function contour_step_CreateFcn(hObject, eventdata, handles)  %#ok
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function stiminterval_lo_Callback(hObject, eventdata, handles)  %#ok

% --- Executes during object creation, after setting all properties.
function stiminterval_lo_CreateFcn(hObject, eventdata, handles) %#ok
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in chk_savemovie.
function chk_savemovie_Callback(hObject, eventdata, handles) %#ok

% --- Executes on button press in chk_savemat.
function chk_savemat_Callback(hObject, eventdata, handles) %#ok

% --- Executes on button press in chk_vesselmask.
function chk_vesselmask_Callback(hObject, eventdata, handles) %#ok

% --- Executes on button press in chk_nolight.
function chk_nolight_Callback(hObject, eventdata, handles) %#ok

% --- Executes on button press in chk_writesigframe.
function chk_writesigframe_Callback(hObject, eventdata, handles) %#ok

function maskthresh_Callback(hObject, eventdata, handles) %#ok
value=str2double(get(hObject,'string'));
if isnan(value) || ( value<0 || value>1 )
    errordlg('Mask threshold value must be between 0 and 1','Incorrect mask threshold');
    set(hObject,'string','0');
end



% --- Executes during object creation, after setting all properties.
function maskthresh_CreateFcn(hObject, eventdata, handles) %#ok
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on button press in chk_selectROI.
function chk_selectROI_Callback(hObject, eventdata, handles) %#ok



function stiminterval_hi_Callback(hObject, eventdata, handles) %#ok


% --- Executes during object creation, after setting all properties.
function stiminterval_hi_CreateFcn(hObject, eventdata, handles) %#ok
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% buttons

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function btn_run_Callback(hObject, eventdata, handles) %#ok
% THIS is the magic button. It scrapes the form, creates the param struct,
% and passes it to the ISI_analysis function.

%We need to stack the structs, because ISI_analysis expects a struct
%containing an array of file parameter structs. We are only running one
%file at a time, so we just have the filesQueue struct within setParams

% Iterate over files
if get(handles.process_all_files, 'value')
    % Get list of all .dat files in directory
    sFileList = dir(fullfile(handles.pathstr,'*.dat'));
else
    sFileList = get(handles.data_filename, 'string');
end

if isempty(sFileList)
    warndlg('No file has been selected.', 'ISI')
    return;
end

if isstruct(sFileList)
    nLoop = length(sFileList);
else
    nLoop = 1;
end

for nFi = 1:nLoop
    if get(handles.process_all_files, 'value')
        % Set file name in GUI
        getdatafile_Callback(hObject, sFileList(nFi).name, handles)
        % Create waitbar
        if ~exist('hWait')
            hWait = waitbar(nFi/length(sFileList), 'Processing all .dat files. Please wait...');
        end
        if ishandle(hWait)
            waitbar(nFi/length(sFileList), hWait);
        else
            break;
        end
    end
    
    %first set the general parameters.
    setParams.saveMovie=get(handles.chk_savemovie,'value');
    setParams.saveFig=get(handles.chk_writesigframe,'value');
    setParams.saveFigAll = get(handles.chk_writesigframeall,'value');
    setParams.saveToMat=get(handles.chk_savemat,'value');
    setParams.useVesselMask=get(handles.chk_vesselmask,'value');
    
    %create the struct for all file parameters
    setParams.filesQueue.path2dir = handles.pathstr;
    setParams.filesQueue.name = get(handles.data_filename,'string');
    setParams.filesQueue.refImage=get(handles.vessel_filename,'string');
    setParams.filesQueue.precision='int16';
    if ~isempty( [get(handles.trials2use_Start, 'string')] ) && ~isempty( [get(handles.trials2use_End, 'string')] )
        setParams.filesQueue.Trials2Use= str2double(get(handles.trials2use_Start, 'string')):str2double(get(handles.trials2use_End, 'string'));
    else
        setParams.filesQueue.Trials2Use = [];
    end
    setParams.filesQueue.Whisker=get(handles.whiskername,'string');
    setParams.filesQueue.minFaces=50; %ignores contours that would be too small
    setParams.filesQueue.maxFaces=1000; %ignores contours that are too big
    
    contour_min=str2double(get(handles.contour_min,'string'));
    contour_step=str2double(get(handles.contour_step,'string'));
    contour_max=str2double(get(handles.contour_max,'string'));
    setParams.filesQueue.ContourLineVals=contour_min:contour_step:contour_max;
    stiminterval_lo=str2double(get(handles.stiminterval_lo,'string'));
    stiminterval_hi=str2double(get(handles.stiminterval_hi,'string'));
    setParams.filesQueue.stimInterval=[stiminterval_lo stiminterval_hi];
    setParams.filesQueue.smoothSigma = str2double(get(handles.smooth_sigma, 'string'));
    
    climinterval_lo=str2double(get(handles.climinterval_lo,'string'));
    climinterval_hi=str2double(get(handles.climinterval_hi,'string'));
    setParams.filesQueue.climAll = [climinterval_lo climinterval_hi];
    
    imgbin_x = str2double(get(handles.imgbin_x, 'string'));
    imgbin_y = str2double(get(handles.imgbin_y, 'string'));
    setParams.filesQueue.imgBin = [imgbin_x imgbin_y];
    
    setParams.filesQueue.preStimDurSec=str2double(get(handles.preStimDurSec,'string'));
    setParams.filesQueue.stimDurSec=str2double(get(handles.stimDurSec,'string'));
    
    setParams.filesQueue.maskthresh = str2double(get(handles.maskthresh, 'string'));
    setParams.filesQueue.selectROI = get(handles.chk_selectROI, 'value');
    setParams.filesQueue.selectSignalROI = get(handles.chk_selectSignalROI, 'value');
    setParams.filesQueue.selectSignalProfile = get(handles.chk_selectSignalProfile, 'value');
    setParams.filesQueue.runBarrelFinder = get(handles.chk_runBarrelFinder, 'value');
    
    % By default, load data, average trials and run analysis
    setParams.filesQueue.DoLoad = true;
    setParams.filesQueue.DoTrialAverage = true;
    setParams.filesQueue.DoAnalyze = true;
    
    % If this function was called via the Info button, instruct next
    % scripts to only display file info, but not load or analyze data
    if hObject == handles.btn_explore
        setParams.filesQueue.DoLoad = false;
        setParams.filesQueue.DoTrialAverage = false;
        setParams.filesQueue.DoAnalyze = false;
    end
    
    % If this function was called via the Load Data button, instruct next
    % scripts to only load data
    if hObject == handles.btn_loaddata
        setParams.filesQueue.DoLoad = true;
        setParams.filesQueue.DoTrialAverage = false;
        setParams.filesQueue.DoAnalyze = false;
    end
    
    % If this function was called via the Average Trial button, instruct next
    % scripts to load data (if not already done so) and then average trials.
    if hObject == handles.btn_averagetrials
        setParams.filesQueue.DoLoad = true;
        setParams.filesQueue.DoTrialAverage = true;
        setParams.filesQueue.DoAnalyze = false;
    end
    
    % if we are using the nolight normalization, set those parameters
    setParams.useNoLight=get(handles.chk_nolight,'value');
    setParams.Rnolight=setParams.filesQueue;
    setParams.Rnolight.name=get(handles.nolight_filename,'string');
    setParams.Rnolight.Whisker='none';
    
    % Run analysis
    try
        ISI_analysis(setParams); %ISI_analysis returns setParams, but we aren't using it
    catch e
        errordlg({['Error using ==> ', e.stack(1).name, ' at ' num2str(e.stack(1).line)],...
            '',e.message,''},...
            'Oops.');
        %rethrow(e);
        break;
    end
    
end
if exist('hWait')
    if ishandle(hWait)
        close(hWait)
    end
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% --- Executes on button press in btn_explore.
function btn_explore_Callback(hObject, eventdata, handles) %#ok
btn_run_Callback(hObject, eventdata, handles)
return

% --- Executes on button press in btn_viewmask.
function btn_viewmask_Callback(hObject, eventdata, handles) %#ok
if isempty(get(handles.vessel_filename,'string'))
    warndlg('No vessel filename given');
else
    ISI_viewMask(fullfile(handles.pathstr, get(handles.vessel_filename,'string')),...
        str2double(get(handles.maskthresh,'string')));
end

function trials2use_Start_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function trials2use_Start_CreateFcn(hObject, eventdata, handles)


if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function trials2use_End_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function trials2use_End_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function framebinsize_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function framebinsize_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in chk_writesigframeall.
function chk_writesigframeall_Callback(hObject, eventdata, handles)

function edit21_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function edit21_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit22_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function edit22_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit23_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function edit23_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in btn_loaddata.
function btn_loaddata_Callback(hObject, eventdata, handles)
btn_run_Callback(hObject, eventdata, handles)
return


% --- Executes on button press in chk_selectSignalROI.
function chk_selectSignalROI_Callback(hObject, eventdata, handles)

% --- Executes on button press in chk_runBarrelFinder.
function chk_runBarrelFinder_Callback(hObject, eventdata, handles)

function climinterval_lo_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function climinterval_lo_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function climinterval_hi_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function climinterval_hi_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function imgbin_x_Callback(hObject, eventdata, handles)
% hObject    handle to imgbin_x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of imgbin_x as text
%        str2double(get(hObject,'String')) returns contents of imgbin_x as a double


% --- Executes during object creation, after setting all properties.
function imgbin_x_CreateFcn(hObject, eventdata, handles)
% hObject    handle to imgbin_x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function imgbin_y_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function imgbin_y_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function smooth_sigma_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function smooth_sigma_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in btn_apply_clim.
function btn_apply_clim_Callback(hObject, eventdata, handles)
% Apply new colormap limits (clim) to all open figures
hFig = findobj(0, 'type', 'figure');
nCLimLo = str2double(get(handles.climinterval_lo, 'string')); % permil
nCLimHi = str2double(get(handles.climinterval_hi, 'string')); % permil
nCLimLo = nCLimLo / 1000; % percent
nCLimHi = nCLimHi / 1000; % percent
% Iterate over figures
for hf = hFig(:)
    set(findobj(hf, 'type', 'axes'), 'clim', [nCLimLo nCLimHi])
end

return


% --- Executes on button press in checkbox10.
function checkbox10_Callback(hObject, eventdata, handles)
return


% --- Executes on button press in btn_averagetrials.
% Load data, if not already done, and compute trial averages
function btn_averagetrials_Callback(hObject, eventdata, handles)
btn_run_Callback(hObject, eventdata, handles)
return


% --- Executes on button press in chk_selectSignalProfile.
function chk_selectSignalProfile_Callback(hObject, eventdata, handles)

return


% --- Executes on button press in process_all_files.
function process_all_files_Callback(hObject, eventdata, handles)

if get(hObject, 'value')
    set(handles.data_filename, 'enable', 'off')
else
    set(handles.data_filename, 'enable', 'on')
end

return


% --- Executes on selection change in popup_plugins.
function popup_plugins_Callback(hObject, eventdata, handles)
% Run selected plugin
sPlugin = get(hObject, 'string');
sPluginId = get(hObject,'value');
% Get plugin path
sPwd = pwd;
sPath = mfilename('fullpath');
vIndx = findstr(filesep, sPath);
sPath = [sPath(1:vIndx(end)) 'plugins' filesep sPlugin{sPluginId} filesep];

% Change current directory to plugin folder as files therein may not be in the Matlab path
cd(sPath)

% Run plugin
try
    eval(sprintf('%s(hObject, eventdata, handles);',sPlugin{sPluginId}))
catch
    errordlg(sprintf('An error occurred when running the %s plugin:\n\n %s', sPlugin{sPluginId}, lasterr))
end
cd(sPwd)

return

% --- Executes during object creation, after setting all properties.
function popup_plugins_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% Get list of plugins
sPath = mfilename('fullpath');
vIndx = findstr(filesep, sPath);
sPath = [sPath(1:vIndx(end)) 'plugins' filesep];
sPwd = pwd;
cd(sPath)
tDir = dir;
cPlugins = {};
entryIsDir = [tDir.isdir]%keep only directories
tDir = tDir(entryIsDir);
for i = 3:length(tDir)
    if tDir(i).isdir
        cPlugins{i-2} = tDir(i).name;
    end
end
cd(sPwd)

% Update popup menu in GUI
set(hObject, 'string', cPlugins)

return
