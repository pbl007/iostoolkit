function idx = ISI_selectContour(XdataLowest,YdataLowest,prmtsCurrent,dipIm, contour_level)
%label each of the contour lines and let user select by chosing from dialog


imRef = imread(fullfile(prmtsCurrent.path2dir,prmtsCurrent.refImage));


h2RoiSelect = figure('Name','Select Isolated barrel');
%Compute size difference between filtered "dip" image and reference
[deltaRC] = size(imRef) - size(dipIm);
deltaR = fix(deltaRC(1)/2);
deltaC = fix(deltaRC(2)/2);
imshow(imRef)
axis image
hold on

%% create region for each 
nROI = numel(XdataLowest);
str = cell(nROI+1,1);
for iROI = 1 : nROI
    x =  XdataLowest{iROI}+deltaR;
    y = YdataLowest{iROI}+deltaC;
    plot(x,y,'b-','LineWidth',2)
    text(mean(x(~isnan(x))),mean(y(~isnan(y))),num2str(iROI),'color','r','fontSize',12);
    str{iROI} = num2str(iROI);
end
str{end} = 'None';
title(['Contour= ' num2str(contour_level)]);

%% list dialog 

cpos=get(h2RoiSelect,'position'); %contour fig position
figpos=[cpos(1)+cpos(3) cpos(2)+cpos(4) 192 386];

[idx]=showlistbox('ListString',str,'PromptString','Select ROI',...
    'position',figpos,'initialvalue',length(str));


if idx == nROI + 1; idx = -1;end %user selecte None -> set idx to -1
if isempty(idx)
    error('ISI_isolateBarrel: cancelled contour selection');
end

%% clean up
close(h2RoiSelect)





%% showlistbox
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [selection, value]=showlistbox(varargin)
%SHOWLISTBOX list selection dialog box, based on listdlg
%   [SELECTION,OK] = LISTDLG('ListString',S) creates a modal dialog box
%   which allows you to select a string or multiple strings from a list.
%   SELECTION is a vector of indices of the selected strings (length 1 in
%   the single selection mode).  This will be [] when OK is 0.  OK is 1 if
%   you push the OK button, or 0 if you push the Cancel button or close the
%   figure.
%   Double-clicking on an item or pressing <CR> when multiple items are
%   selected has the same effect as clicking the OK button.  Pressing <CR>
%   is the same as clicking the OK button. Pressing <ESC> is the same as
%   clicking the Cancel button.
%
%   Inputs are in parameter,value pairs:
%
%   Parameter       Description
%   'ListString'    cell array of strings for the list box.
%   'ListSize'      [width height] of listbox in pixels; defaults
%                   to [160 300].
%   'InitialValue'  vector of indices of which items of the list box
%                   are initially selected; defaults to the first item.
%   'Name'          String for the figure's title; defaults to ''.
%   'PromptString'  string matrix or cell array of strings which appears 
%                   as text above the list box; defaults to {}.
%   'OKString'      string for the OK button; defaults to 'OK'.
%   'CancelString'  string for the Cancel button; defaults to 'Cancel'.
%   'Position'      [left bottom width height], bases the actual 
%                   width & height on ListSize, default vals are w=192,h=386
%
%mjp 2011.09.30     a lot of the code is borrowed directly from listdlg

figname = '';
listsize = [160 300];
promptstring = {};
liststring = [];
initialvalue = [];
okstring = 'OK';
cancelstring = 'Cancel';
fus = 8; %frame/uicontrol height
ffs = 8; %frame/figure spacing
uh = 22; %uicontrol button height


ex = get(0,'defaultuicontrolfontsize')*1.7;  % height extent per line of uicontrol text (approx)
fp = get(0,'defaultfigureposition');
w = 2*(fus+ffs)+listsize(1);
h = 2*ffs+6*fus+ex*length(promptstring)+listsize(2)+uh;
fp = [fp(1) fp(2)+fp(4)-h w h];  % keep upper left corner fixed

if mod(length(varargin),2) ~= 0
    % input args have not com in pairs, woe is me
    error('MATLAB:listdlg:InvalidArgument', 'Arguments to LISTDLG must come param/value in pairs.')
end
for i=1:2:length(varargin)
    switch lower(varargin{i})
     case 'name'
      figname = varargin{i+1};
     case 'promptstring'
      promptstring = varargin{i+1};
     case 'listsize'
      listsize = varargin{i+1};
     case 'liststring'
      liststring = varargin{i+1};
     case 'initialvalue'
      initialvalue = varargin{i+1};
     case 'position'
      newfp = varargin{i+1}; %taking initial left and bottom
      %default w=192, h=386
      fp = [newfp(1) newfp(2)+fp(4)-h w h];  % keep upper left corner fixed
    end
end
    
%set defaults
if isempty(liststring)
    error('showlistbox:NeedParameter', 'ListString parameter is required.')
end
if isempty(initialvalue)
    initialvalue = 1;
end
if ischar(promptstring)
    promptstring = cellstr(promptstring); 
end

fig_props = { ...
    'name'                   figname ...
    'color'                  get(0,'defaultUicontrolBackgroundColor') ...
    'resize'                 'off' ...
    'numbertitle'            'off' ...
    'menubar'                'none' ...
    'windowstyle'            'modal' ...
    'visible'                'off' ...
    'createfcn'              ''    ...
    'position'               fp   ...
    'closerequestfcn'        'delete(gcbf)' ...
            };

liststring=cellstr(liststring);
fig = figure(fig_props{:});
if ~isempty(promptstring)
    prompt_text = uicontrol('style','text','string',promptstring,...
        'horizontalalignment','left',...
        'position',[ffs+fus fp(4)-(ffs+fus+ex*length(promptstring))+fus ...
        listsize(1) ex*length(promptstring)]); %#ok
end
btn_wid = (fp(3)-2*(ffs+fus)-fus)/2;
listbox = uicontrol('style','listbox',...
                    'position',[ffs+fus ffs+uh+4*fus listsize],...
                    'string',liststring,...
                    'backgroundcolor','w',...
                    'max',2,... %allows multiple selection
                    'tag','listbox',...
                    'value',initialvalue, ...
                    'callback', {@doListboxClick});
ok_btn = uicontrol('style','pushbutton',...
                   'string',okstring,...
                   'position',[ffs+fus ffs+fus btn_wid uh],...
                   'callback',{@doOK,listbox});
cancel_btn = uicontrol('style','pushbutton',...
                       'string',cancelstring,...
                       'position',[ffs+2*fus+btn_wid ffs+fus btn_wid uh],...
                       'callback',{@doCancel,listbox});

set([fig, ok_btn, cancel_btn, listbox], 'keypressfcn', {@doKeypress, listbox});
% make sure we are on screen
movegui(fig)
set(fig, 'visible','on'); drawnow;

try
    % Give default focus to the listbox *after* the figure is made visible
    uicontrol(listbox);
    uiwait(fig);
catch e
    if ishandle(fig)
        delete(fig)
    end
    rethrow(e)
end
if isappdata(0,'ListDialogAppData__')
    ad = getappdata(0,'ListDialogAppData__');
    selection = ad.selection;
    value = ad.value;
    rmappdata(0,'ListDialogAppData__')
else
    % figure was deleted
    selection = [];
    value = 0;
end


%% figure, OK and Cancel KeyPressFcn
function doKeypress(src, evd, listbox) %#ok
switch evd.Key
 case 'escape'
  doCancel([],[],listbox);
end

%% OK callback
function doOK(ok_btn, evd, listbox) %#ok
if (~isappdata(0, 'ListDialogAppData__'))
    ad.value = 1;
    ad.selection = get(listbox,'value');
    setappdata(0,'ListDialogAppData__',ad);
    delete(gcbf);
end

%% Cancel callback
function doCancel(cancel_btn, evd, listbox) %#ok
ad.value = 0;
ad.selection = [];
setappdata(0,'ListDialogAppData__',ad)
delete(gcbf);

%% SelectAll callback
function doSelectAll(selectall_btn, evd, listbox) %#ok
set(selectall_btn,'enable','off')
set(listbox,'value',1:length(get(listbox,'string')));

%% Listbox callback
function doListboxClick(listbox, evd, selectall_btn) %#ok
% if this is a doubleclick, doOK
if strcmp(get(gcbf,'SelectionType'),'open')
    doOK([],[],listbox);
elseif nargin == 3
    if length(get(listbox,'string'))==length(get(listbox,'value'))
        set(selectall_btn,'enable','off')
    else
        set(selectall_btn,'enable','on')
    end
end





