% function MAP_fig_1_multi(mname,config)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ################################################################# %% DESCRIPTION
% Script to analyse data for:
% Grieves et al. (2020) 
% Quantification of firing rate map procedures
%
% This script 
% NOTE: This script is dependent on graphtest to load and prepare data and indices
%
% Makes:
% 
%

% HISTORY:
% version 1.0.0, Release 12/07/20 Initial release/comments added (script generated before this date)
%
% Author: Roddy Grieves
% UCL, 26 Bedford Way
% eMail: r.grieves@ucl.ac.uk
% Copyright 2020 Roddy Grieves
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ################################################################# %% FUNCTION BODY
%% mapping settings
    mapset                                              = struct;
    mapset.padding                                      = 0; % in mm
    mapset.mindwell                                     = 0; % (default 0.1) total number of seconds that rat has to be in a bin for it to be filled, for adaptive method this controls how big the bins will be expanded
    mapset.mindist                                      = 40; % (mm, default 1) used by adaptive and gaussian method, only bins within this distance of some position data will be filled, bigger means more filled (blue) bins
    mapset.maxdist                                      = 640; % (mm, default 50) used by kadaptive, adaptive and gaussian method, only bins within this distance of some position data will be filled, bigger means more filled (blue) bins
    mapset.srate                                        = 50; % (default 50) sampling rate of data in Hz, used to calculate time
    mapset.steps                                        = 32; % the number of convolution size steps to use for kadaptive    
    mapset.kern                                         = 'biweight'; % kernel
    mapset.smethod                                      = 1; % smoothing method, 1 = before division, 2 = after, 3 = no smoothing
    mapset.twindow                                      = 0.25; % time window (s) over which to estimate instantaneous firing for temporal methods

%% ##################################### Heading 2
%% #################### Heading 3
%% Heading 4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ################################################################# %% FUNCTION BODY
%% #################### Heatmap of error by binsize and smoothing
    % Figure settings
    fig1 = figure('visible','on','Position',[50,60,1100,920]); 
    set(gcf,'InvertHardCopy','off'); % this stops white lines being plotted black      
    set(gcf,'color','w'); % makes the background colour white
    colormap(jet(256)); % to make sure the colormap is not the horrible default one

    % method settings
    ynow = 720;
    xnow = 50;
    psiz = 90;
    xbuff = psiz+20+psiz;
    xvec = xnow:xbuff:10000;
    ybuff = psiz+20;
    yvec = ynow:-ybuff:-1000;
    pvec = [3 4 7 8 9];
    mbuff = 540;

    d = 25; % duration of sessions, mins
    r = 300; % average radius of fields, mm

    ename = 'arena120cm'; 
    mnames = {'histogram','ash','kadaptive','ksde','fyhn','kyadaptive'};    
    dnames = {'Histogram','ASH','Adaptive smoothing','KSDE','tKSDE','Adaptive binning'}; 

    mem = 1;
    load('example_pcells.mat','dat');    
    maps_for_plotting = cell(size(dat,1),2,length(mnames));
    fname = ['fig_ex_pcells.mat'];
    override_maps = 0;
    if exist(fname,'file') && ~override_maps
        mem = 1;
        disp(sprintf('\t\t...loading old maps'))
        load(fname,'maps_for_plotting')
    else
        disp(sprintf('\t\t...generating maps'))            
        mem = 0;
    end

    %% plot spikes and positions
    for pp = 1:length(pvec) % for every example cell 
        pos = dat.pos{pvec(pp)}(:,1:2);
        spk = dat.spk{pvec(pp)}(:,1:2);
        ax = axes('Units','pixels','Position',[xvec(pp)+35 yvec(1)+25 120 120]);  
            plot(pos(:,1),pos(:,2),'k'); hold on;
            plot(spk(:,1),spk(:,2),'r.','MarkerSize',6);
            axis xy off
            daspect([1 1 1])
            ax.XLim = [min(pos(:,1)) max(pos(:,1))];
            ax.YLim = [min(pos(:,2)) max(pos(:,2))];
            text(0,1.1,sprintf('Cell %d',pp),'Units','normalized','FontSize',12,'HorizontalAlignment','left')
            
            % if pp==3
            %     text(0.5,-0.25,sprintf('Balanced solution parameters'),'Units','normalized','FontSize',12,'HorizontalAlignment','center')
            % end
            if pp==1
                % text(0.5,-0.1,sprintf('120 mm'),'Units','normalized','FontSize',10,'HorizontalAlignment','center')
                text(-0.12,0.5,sprintf('1200 mm'),'Units','normalized','FontSize',10,'HorizontalAlignment','center','Rotation',90)
            end

        % ax = axes('Units','pixels','Position',[xvec(pp)+mbuff yvec(1)+30 psiz psiz]);  
        %     plot(pos(:,1),pos(:,2),'k'); hold on;
        %     plot(spk(:,1),spk(:,2),'r.','MarkerSize',6);
        %     axis xy off
        %     daspect([1 1 1])
        %     ax.XLim = [min(pos(:,1)) max(pos(:,1))];
        %     ax.YLim = [min(pos(:,2)) max(pos(:,2))];
        %     text(0,1.1,sprintf('Cell %d',pp),'Units','normalized','FontSize',12,'HorizontalAlignment','left')            
            % if pp==3
            %     text(0.5,-0.25,sprintf('Minimum error solution parameters'),'Units','normalized','FontSize',12,'HorizontalAlignment','center')
            % end            
    end

    %% generate and plot firing rate maps
    for mm = 1:length(mnames) % for every method
        disp(sprintf('\t%s...',mnames{mm}))
        mname = mnames{mm};

        % collect regression results
        fs = [4000 8000 16000];
        rfs = [127 179 253];
        switch mname
            case {'ksde'}
                ds = [4 16 24];
            otherwise
                ds = [4 16 64];
        end          
        datn = [];
        for f = 1:length(fs)
            for d = 1:length(ds)
                nme = [ename '_' mname '_sigma' num2str(fs(f)) '_duration' num2str(ds(d)) '_gamultiobj.mat'];
                if exist(nme,'file')
                    load(nme);
                    datn = [datn; rfs(f) ds(d) balance_solution1(1:2) minerr_solution(1:2)];
                end
            end
        end

        % balanced settings
        [b,s,r] = mvregress([ones(size(datn,1),1),datn(:,1:2)],datn(:,3:4)); 
        f = @(xy) [b(1,1) + xy(:,1).*b(2,1) + xy(:,2).*b(3,1), b(1,2) + xy(:,1).*b(2,2) + xy(:,2).*b(3,2)]; % [field size, duration] > [bin size, smoothing]
        balanced_params = f([300,25]);

        % minimum settings
        if strcmp(mname,'histogram')
            [b,s,r] = mvregress([ones(size(datn,1),1),datn(:,1:2)],datn(:,6)); 
            f = @(xy) [b(1,1) + xy(:,1).*b(2,1) + xy(:,2).*b(3,1)]; % [field size, duration] > [smoothing]
            minimum_params = [1 f([300,25])];
        elseif strcmp(mname,'ksde')
            [b,s,r] = mvregress([ones(size(datn,1),1),datn(:,1:2)],datn(:,6)); 
            f = @(xy) [b(1,1) + xy(:,1).*b(2,1) + xy(:,2).*b(3,1)]; % [field size, duration] > [smoothing]
            minimum_params = [2.5 f([300,25])];
        else
            [b,s,r] = mvregress([ones(size(datn,1),1),datn(:,1:2)],datn(:,5:6)); 
            f = @(xy) [b(1,1) + xy(:,1).*b(2,1) + xy(:,2).*b(3,1), b(1,2) + xy(:,1).*b(2,2) + xy(:,2).*b(3,2)]; % [field size, duration] > [bin size, smoothing]
            minimum_params = f([300,25]);
        end

        for pp = 1:length(pvec) % for every example cell
            disp(sprintf('\t\tcell %d...',pp))
            
            if mem
                ratemap = maps_for_plotting{pp,1,mm};
            else    
                pos = dat.pos{pvec(pp)}(:,1:2);
                spk = dat.spk{pvec(pp)}(:,1:2);
                rmset = mapset; % rate mapper settings structure - passing mapset directly to ratemapper causes a memory leak
                rmset.method = mname;
                rmset.binsize = balanced_params(1);
                rmset.ssigma = balanced_params(2);
                rmset.pot = dat.pos{pvec(pp)}(:,3);
                rmset.spt = dat.spk{pvec(pp)}(:,3);
                
                if strcmp(mname,'fyhn')
                    rmset.maxdist = 320;
                    rmset.mindist = 50;                                
                end
                [ratemap,dwellmap,spikemap,rmset,speedlift] = rate_mapper(pos,spk,rmset);
                maps_for_plotting{pp,1,mm} = ratemap;
            end            

            ax2 = axes('Units','pixels','Position',[xvec(pp) yvec(mm+1) psiz psiz]);  
                im = imagesc(ratemap,'alphadata',~isnan(ratemap)); hold on;
                daspect([1 1 1])
                axis xy off tight
                colormap(gca,turbo) 
                if pp==1
                    if mm==1
                        text(-0.1,1.25,sprintf('%s',dnames{mm}),'Units','normalized','FontSize',12)
                    else
                        text(-0.1,1.1,sprintf('%s',dnames{mm}),'Units','normalized','FontSize',12)
                    end
                end
                if mm==1
                    text(0,1.1,sprintf('Balanced'),'Units','normalized','FontSize',10)                        
                end

                % colorbar
                if mm==length(mnames) && pp==1
                    axc = axes('Units','pixels','Position',[ax2.Position(1) ax2.Position(2)-30 ax2.Position(3)+40 12]); 
                        mat = linspace(ax2.CLim(1),ax2.CLim(2),100);
                        imagesc([ax2.CLim(1),ax2.CLim(2)],[1 1],mat);
                        colormap(axc,ax2.Colormap);
                        axis xy
    
                        axc.XTick = [ax2.CLim(1),ax2.CLim(2)];
                        axc.XTickLabel = {'0','Max'};
                        axc.YTick = [];
                        axc.XAxisLocation = 'bottom';
                        text(0.5,1.5,sprintf('Firing rate (Hz)'),'FontSize',8,'HorizontalAl','center','Units','normalized')
                end 

            if mem
                ratemap = maps_for_plotting{pp,2,mm};
            else    
                pos = dat.pos{pvec(pp)}(:,1:2);
                spk = dat.spk{pvec(pp)}(:,1:2);
                rmset = mapset; % rate mapper settings structure - passing mapset directly to ratemapper causes a memory leak
                rmset.method = mname;
                rmset.binsize = minimum_params(1);
                rmset.ssigma = minimum_params(2);
                rmset.pot = dat.pos{pvec(pp)}(:,3);
                rmset.spt = dat.spk{pvec(pp)}(:,3);
                [ratemap,dwellmap,spikemap,rmset,speedlift] = rate_mapper(pos,spk,rmset);
                maps_for_plotting{pp,2,mm} = ratemap;
            end            

            ax2 = axes('Units','pixels','Position',[xvec(pp)+psiz+5 yvec(mm+1) psiz psiz]);  
                im = imagesc(ratemap,'alphadata',~isnan(ratemap)); hold on;
                daspect([1 1 1])
                axis xy off tight
                colormap(gca,turbo)  

                if mm==1
                    text(0,1.1,sprintf('Minimum'),'Units','normalized','FontSize',10)                        
                end


        end
% return
    end
    if ~mem
        save(fname,'maps_for_plotting')
    end

        % return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ################################################################# %% Save the figure
    % Save the figure    
    disp(sprintf('\tSaving figure...'))    
    figname = [fig_dir2 '\S9_Fig.png'];
    [~,~,~] = mkdir(fig_dir2);    
    exportgraphics(fig1,figname,'BackgroundColor','w','ContentType','image','Resolution',350);  
    close(fig1)    
    
    


  
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    