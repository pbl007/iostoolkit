%% plot all IOS contour maps, given file list

% basepath='D:\dklab\Mike IOS data\';
basepath='E:\dklab\Mike IOS data\';
% basepath='C:\Users\Mike\Desktop\dklab\Mike IOS data\';
%basepath='C:\Documents and Settings\Mike\Desktop\dklab\Mike IOS data\';

pathdate='2011.09.29.mjp';

vesselimg='blue focusstack.png';
%have center whisker listed first, then in list in circle around it

plotcentroid=1; %if true, plot the centroid
usepatch=0; %if true, plot a transparent patch, else plot contour line

% %2011.10.11.mjp
% barreldata={...
%     'C1_171221.mat',...
%     'C2_170038.mat',...
%     'C3_172127.mat',...
%     'D1_173137.mat',...
%     'D2_173844.mat',...
%     'D3_180603.mat',...
%     'D4_181843.mat',...
%     'E1_184812.mat',...
%     'E2_190021.mat',...
%     'E3_190941.mat',...
%     };
%


%2011.09.29.mjp
barreldata={ ...
    'B1_211743.mat',...
    'B2_212715.mat',...
    'B3_213852.mat',...
    'C1_200416.mat',...
    'C2_202042.mat',...
    'C3_204626.mat',...
    'C4_210150.mat',...
    'D1_195451.mat',...
    'D2_190713.mat',...
    'D3_191847.mat',...
    'D4_193007.mat',...
    'E1_184553.mat',...
    'E2_180324.mat',...
    'E3_181416.mat',...
    };

%barreldata={ ...
%    'C2_153342.mat',...       %C2
%    'D3_151305.mat',...       %C1
%};

% barreldata={ ...
%     'C2_131731.mat',...       %C2
%     'C1_132222.mat',...       %C1
%     'B1_132856.mat',...       %B1
%     '',...       %B2
%     'B3_134344.mat',...       %B3
%     'C3_135022.mat',...       %C3
%     '',...       %D3
%     'D2_131205.mat',...       %D2
%     '',...       %D1
%     };
%

% barreldata={ ...
%     'C2_122303.mat',...       %C2
%     'C3_130654.mat',...       %C3
%     'D3_131124.mat',...       %D3
%     'D2_124428.mat',...       %D2
%     'D1_123605.mat',...       %D1
% 	'gamma_125245.mat',...
%     'C1_123039.mat',...       %C1
%     };

% %%2011.07.11
% barreldata={ ...
%     'C2_4s_162231.mat',...       %C2
%     'C1_4s_162941.mat',...       %C1
%     'B1_4s_164922.mat',...       %B1
%     'B2_4s_164244.mat',...       %B2
%     '',...       %B3
%     'C3_4s_165625.mat',...       %C3
%     'D3_4s_161713.mat',...       %D3
%     'D2_4s_155956.mat',...       %D2
%     'D1_4s_161029.mat',...       %D1
%     };

%2011.06.21
% barreldata={ ...
%     'C2_4stim_160358.mat',...       %C2
% 	'C1_4stim_151508.mat',...
%     'B1_4stim_161533.mat',...       %B1
%     'B2_4stim_163149.mat',...       %B2
%     'B3_4stim_164336.mat',...       %B3
%     'C3_4stim_155736.mat',...       %C3
% 	'D4_4stim_172002.mat',...
%     'D3_4stim_171429.mat',...       %D3
%     'D2_4stim_170707.mat',...       %D2
%     'D1_4stim_170128.mat',...       %D1
% 	'gamma_4stim_165531.mat',...
%     };


% barreldata={'placeholder','rat6_c1.mat'};

nfile=length(barreldata);
colorM=hsv(nfile);

%choose 9 colors
% white (center), red, orange, yellow, green, cyan, teal, blue, purple
% colorM=[...
%     %1 1 1; %white
%     1 0 0; %red
%     1 0.5 0; %orange
%     1 1 0; %yellow
%     0 1 0; %green
%     0 1 1; %cyan
%     0 0.6 0.6; %teal
%     0 0 0.8; %blue
%     0.75 0 0.75; %purple
%     ];


%% load image

%load image
imfile=fullfile(basepath,pathdate,vesselimg);
if ~exist(imfile,'file')
    error('Incorrect blood vessel image filename.');
end
Ivessel=imread(imfile);
figure;
imshow(Ivessel);
drawnow;

%% run through contours and plot

for fi=1:nfile
    filename=fullfile(basepath,pathdate,barreldata{fi});
    if ~exist(filename,'file') || isempty(barreldata{fi})
        disp(['File skipped: ' barreldata{fi} ' doesn''t exist']);
        continue;
    end
    
    load(filename);
    if ~exist('ISIdata','var')
        error(['Invalid data, check contents of file ' barreldata{fi}]);
    end
    
    %contour line plot; contour coordinates are based on reference image
    x=ISIdata.contourM(:,1);
    y=ISIdata.contourM(:,2);
    name=strtok(barreldata{fi},'_');
    
    if usepatch
        alpha=0.5;
        patch(x(~isnan(x)), y(~isnan(y)), colorM(mod(fi-1,size(colorM,1))+1,:),'edgecolor','none', 'edgealpha',alpha,'facealpha',alpha);
        text(mean(x(~isnan(x)))+5,mean(y(~isnan(y))),name(1:2),'color','k','fontSize',10)
    else
        line(x, y, 'color', colorM(mod(fi-1,size(colorM,1))+1,:), 'linewidth',2);
        text(mean(x(~isnan(x)))+5,mean(y(~isnan(y))),name(1:2),'color',colorM(mod(fi-1,size(colorM,1))+1,:),'fontSize',10)
    end
    
    
    if plotcentroid
        %make mask from contour and find centroid
        maskx=x(~isnan(x));
        masky=y(~isnan(y));
        BW=poly2mask(maskx,masky,size(Ivessel,1),size(Ivessel,2));
        stats=regionprops(BW,'Centroid');
        if length(stats)>1 %we have multiple contours, so cheat
            clear stats;
            stats.Centroid(1)=mean(x(~isnan(x)));
            stats.Centroid(2)=mean(y(~isnan(y)));
        end
        h=line(stats.Centroid(1),stats.Centroid(2),'color','k','marker','o',...
            'markersize',3,'MarkerFaceColor','k');
        %uistack(h,'bottom');
    end
    
    %patch(ISIdata.contourM(:,1),ISIdata.contourM(:,2),...
    %   colorM(mod(fi-1,size(colorM,1))+1,:), 'edgealpha',alpha,'facealpha',alpha);
    
    %transparent patch
    %     if isfield(ISIdata,'contourMROI')
    %         plotcontourmatrix(ISIdata.contourMROI, ISIdata.contouroffset+ISIdata.ROIcoords(1:2),...
    %             colorM(mod(fi-1,size(colorM,1))+1,:), 1, 0.6 );
    %     else
    %         plotcontourmatrix(ISIdata.contourM, ISIdata.contouroffset,...
    %             colorM(mod(fi-1,size(colorM,1))+1,:), 1, 0.6 );
    %     end
    pause(0.01); drawnow;
    clear ISIdata;
    
end


title(pathdate,'interpreter','none');


