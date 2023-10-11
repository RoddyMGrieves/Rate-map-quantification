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
    % get position data (random walks)
    nwalks = config.nwalks;
    all_pos = cell(nwalks,1);
    for ww = 1:nwalks
        wlength = 64;
        wdir = [scratch_space '\' config.environs{1} '_' num2str(wlength) '_walk' num2str(ww) '.mat'];
        load(wdir,'pox','poy','pot','epoly','emat','wlength'); % load walk data       
        
        pox = fillmissing(pox(:),'linear');
        poy = fillmissing(poy(:),'linear');
        
        duration = numel(pox)*(1/50);
        all_pos(ww,1) = {[pox(:) poy(:) pot(:)]};
    end

    if ~exist('sdata','var')
        sdir = [scratch_space '\' config.environs{1} '_sigma16000_duration64_sdata.mat'];            
        disp(sprintf('\t\t...loading %s',sdir));            
        load(sdir,'sdata'); 
    end

    ename = 'arena120cm'; 
    override_pareto = 0;
    override_maps = 0;
    rmap_colormap = 'turbo';

%% ##################################### Heading 2
%% #################### Heading 3
%% Heading 4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ################################################################# %% FUNCTION BODY
%% #################### Heatmap of error by binsize and smoothing
    % Figure settings
    fig1 = figure('visible','on','Position',[50,60,1400,1050]); 
    set(gcf,'InvertHardCopy','off'); % this stops white lines being plotted black      
    set(gcf,'color','w'); % makes the background colour white
    colormap(jet(256)); % to make sure the colormap is not the horrible default one

%% #################### Heatmap of error by binsize and smoothing
    %%%%%%%%%%%%%%%%%%%%%
    % excell = 69;
    excell = 83;
    % for ee = 93:128
    % % excell = 92;
    % % Figure settings
    % fig1 = figure('visible','on','Position',[50,60,1400,1050]); 
    % set(gcf,'InvertHardCopy','off'); % this stops white lines being plotted black      
    % set(gcf,'color','w'); % makes the background colour white
    % colormap(jet(256)); % to make sure the colormap is not the horrible default one
    % excell = ee;

    xnow = 70;
    ynow = 820;
    siz = 120; 
    spac = 100;
    xvec = xnow : siz+spac : 1000;
    yvec = ones(size(xvec))*ynow;

    ax1 = axes('Units','pixels','Position',[xvec(1) yvec(1) siz siz],'Color','none','Clipping','off');         
        pmap = sdata.pmap_1mm{excell,1};
        im = imagesc(pmap); hold on;
        daspect([1 1 1])
        axis xy off tight
        colormap(gca,rmap_colormap)

        ah = add_panel_title('a',sprintf('Spike probability map'),'yoffset',-5,'xoffset',0,'width',300);    

        % colorbar
        axc = axes('Units','pixels','Position',[ax1.Position(1)+150 ax1.Position(2)+20 12 ax1.Position(4)-20]);
            mat = (linspace(0,100,100))';
            imagesc(mat,ones(size(mat)),mat);
            colormap(axc,rmap_colormap);
            axis xy

            axc.YTick = [0.01 100];
            axc.YTickLabel = {'0','Max'};
            axc.XTick = [];
            axc.YAxisLocation = 'right';
            axc2 = axes('Units','pixels','Position',[axc.Position(1) axc.Position(2)-20 axc.Position(3) axc.Position(3)]); 
            axis on; box on; axc2.XTick = []; axc2.YTick = 0.5; axc2.YTickLabel = {'Unvisited'}; axc2.YAxisLocation = 'right';

    % cut the position and spike data to the trial length
    % we want
    pindx = sdata.pos_index(excell);
    pos_now = all_pos{pindx,1};
    pox = pos_now(:,1);
    poy = pos_now(:,2);
    pot = pos_now(:,3); 

    spk = sdata.spk{excell};
    spx = spk(:,1);
    spy = spk(:,2);
    spt = spk(:,3);

    tcut = 4*60;
    ppox = pox(pot<tcut);
    ppoy = poy(pot<tcut);
    ppot = pot(pot<tcut);            
    pspx = spx(spt<tcut);
    pspy = spy(spt<tcut);
    pspt = spt(spt<tcut);            
    pos = [ppox ppoy]; % positions in mm
    spk = [pspx pspy]; % spikes in mm

    ax1 = axes('Units','pixels','Position',[xvec(2) yvec(1) siz siz],'Color','none','Clipping','off');         
        plot(ppox,ppoy,'Color',[.5 .5 .5 .5]); hold on;
        scatter(pspx,pspy,20,'r','filled','MarkerEdgeColor','none','MarkerFaceAlpha',0.5);
        daspect([1 1 1])
        axis xy off  

        ah = add_panel_title('b',sprintf('Positions and spikes'),'yoffset',-5,'xoffset',10,'width',300);    

%% #################### 
    %%%%%%%%%%%%%%%%%%%%%
    ax1 = axes('Units','pixels','Position',[xnow 400 siz siz],'Color','none','Clipping','off');         
        rmset = mapset; % rate mapper settings structure - passing mapset directly to ratemapper causes a memory leak
        rmset.method = 'histogram';
        rmset.binsize = 40;
        rmset.ssigma = 30;
        rmset.bmethod = 22;
        rmset.maplims = [-max(abs(epoly(:,1))) -max(abs(epoly(:,2))) max(abs(epoly(:,1))) max(abs(epoly(:,2)))]; 
        [ratemap,dwellmap,spikemap,rmset,speedlift] = rate_mapper(pos,spk,rmset);        
% rmset.maplims
        im = imagesc(ratemap,'alphadata',~isnan(ratemap)); hold on;
        daspect([1 1 1])
        axis xy off tight
        colormap(gca,'turbo')  

    ax1 = axes('Units','pixels','Position',[xnow+siz+50 400 siz siz],'Color','none','Clipping','off');   
        bs = range(epoly(:,1))/rmset.bmethod;
        rmset.maplims = [-max(abs(epoly(:,1))) -max(abs(epoly(:,2))) max(abs(epoly(:,1))) max(abs(epoly(:,2)))]-[bs 0 bs 0]./2; 
        [ratemap,dwellmap,spikemap,rmset,speedlift] = rate_mapper(pos,spk,rmset);        
% rmset.maplims
        im = imagesc(ratemap,'alphadata',~isnan(ratemap)); hold on;
        daspect([1 1 1])
        axis xy off tight
        colormap(gca,'turbo')  

    ax1 = axes('Units','pixels','Position',[xnow+siz*2+100 400 siz siz],'Color','none','Clipping','off');         
        rmset.maplims = [-max(abs(epoly(:,1))) -max(abs(epoly(:,2))) max(abs(epoly(:,1))) max(abs(epoly(:,2)))]-[0 bs 0 bs]./2; 
        [ratemap,dwellmap,spikemap,rmset,speedlift] = rate_mapper(pos,spk,rmset);        
% rmset.maplims
        im = imagesc(ratemap,'alphadata',~isnan(ratemap)); hold on;
        daspect([1 1 1])
        axis xy off tight
        colormap(gca,'turbo')  
    % end











