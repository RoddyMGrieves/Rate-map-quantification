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
    font_siz = 12;

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
    tfonts = [20 12];
    fsiz = 12;

%% #################### Heatmap of error by binsize and smoothing
    %%%%%%%%%%%%%%%%%%%%%
    excell = 71;

    xnow = 70;
    ynow = 820;
    siz = 120; 
    spac = 100;
    xvec = xnow : siz+spac : 1000;
    yvec = ones(size(xvec))*ynow;

    ax1 = axes('Units','pixels','Position',[xvec(1) yvec(1) siz siz],'Color','none','Clipping','off','FontSize',fsiz);         
        pmap = sdata.pmap_1mm{excell,1};
        im = imagesc(pmap); hold on;
        daspect([1 1 1])
        axis xy off tight
        colormap(gca,rmap_colormap)

        ah = add_panel_title('A',sprintf('Spike probability map'),'yoffset',-5,'xoffset',0,'width',300,'fontsize',tfonts);    

        % colorbar
        axc = axes('Units','pixels','Position',[ax1.Position(1)+150 ax1.Position(2)+20 12 ax1.Position(4)-20],'FontSize',fsiz);
            mat = (linspace(0,100,100))';
            imagesc(mat,ones(size(mat)),mat);
            colormap(axc,rmap_colormap);
            axis xy

            axc.YTick = [0.01 100];
            axc.YTickLabel = {'0','Max'};
            axc.XTick = [];
            axc.YAxisLocation = 'right';
            axc.FontSize = 8;
            axc2 = axes('Units','pixels','Position',[axc.Position(1) axc.Position(2)-20 axc.Position(3) axc.Position(3)]); 
            axis on; box on; axc2.XTick = []; axc2.YTick = 0.5; axc2.YTickLabel = {'Unvisited'}; axc2.YAxisLocation = 'right';
            axc2.FontSize = 8;

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

    tcut = 16*60;
    ppox = pox(pot<tcut);
    ppoy = poy(pot<tcut);
    ppot = pot(pot<tcut);            
    pspx = spx(spt<tcut);
    pspy = spy(spt<tcut);
    pspt = spt(spt<tcut);            
    pos = [ppox ppoy]; % positions in mm
    spk = [pspx pspy]; % spikes in mm

    ax1 = axes('Units','pixels','Position',[xvec(2) yvec(1) siz siz],'Color','none','Clipping','off');         
        p1 = plot(ppox,ppoy,'Color',[.5 .5 .5 .5]); hold on;
        s1 = scatter(pspx,pspy,20,'r','filled','MarkerEdgeColor','none','MarkerFaceAlpha',0.5);
        daspect([1 1 1])
        axis xy off  

        ah = add_panel_title('B',sprintf('Positions and spikes'),'yoffset',-5,'xoffset',10,'width',300,'fontsize',tfonts);    

        [~,leg] = legendflex([p1 s1]...
            ,{'trajectory','spikes'}...
            ,'anchor',{'e','e'},'ncol',1,'box','off','buffer',[120,0],'xscale',1,'fontsize',font_siz); 

%% #################### Heatmap of error by binsize and smoothing
    % method settings
    xnow = 70;
    ynow = ynow-100;
    mnames = {'histogram','ash','kadaptive','ksde','fyhn','kyadaptive'};
    dnames = {'Histogram','Averaged Shifted Histogram (ASH)','Adaptive Smoothing','Kernel Smoothed Density Estimate (KSDE)','Temporal KSDE','Adaptive Binning'}; 
    ls = {'C','D','E','F','G','H','I'};
    mapidx = 4;
    xsiz = 200;
    xbuff = 200;
    xvec = xnow : (xsiz+xbuff) : 1000;
    ysiz = 200;
    ybuff = 195;
    yvec = ynow : -(ysiz+ybuff) : 0;
    [xx,yy] = meshgrid(xvec,yvec);
    xx = xx';
    yy = yy';

    % plot settings
    n_color_levels = 64;
    error_colormap = flipud(inferno(n_color_levels));    
    clims = [1*10^-12 1*10^-11];   
    interp_method = 'nearest';  
  
    maps_for_plotting = cell(length(mnames),7,7);
    fname = [fig_dir 'histogram_fig_maps_maps.mat'];
    override_maps = 0;
    if exist(fname,'file') && ~override_maps
        mem = 1;
        disp(sprintf('\t\t...loading old maps'))
        load(fname,'maps_for_plotting')
    else
        disp(sprintf('\t\t...generating maps'))            
        mem = 0;
    end

    for mm = 1:length(mnames) % for every method
        mname = mnames{mm};
        sigma_now = '4000';
        time_now = '16'; 

        psiz = 42;
        buff = 47;
        bvec = xx(mm):buff:10000;
        svec = yy(mm):-buff:-1000;
        smoo_spac = -0.85;        
        switch mname
            case {'histogram'}
                bvnow = [1 10 20 40 80 160 320];
                svnow = [1 10 20 40 80 160 320];
                smoo_text = sprintf('Smoothing %c (mm)',963);
                grid_text = sprintf('Bin size (mm)');
            case {'ash'}
                bvnow = [20 40 60 80 160 320 640];
                svnow = [2 4 8 16 32 48 64]; 
                smoo_text = 'Smoothing {\itm}';
                grid_text = sprintf('Bin size (mm)');                
            case {'ksde'}
                bvnow = [2.5 5 10 20 40 160 320];
                svnow = [10 20 30 40 80 160 320];  
                smoo_text = sprintf('Smoothing bandwidth (mm)');
                grid_text = sprintf('Bin size (mm)');                
            case {'kyadaptive'}
                bvnow = [2.5 5 10 20 80 160 320];
                svnow = [0.5 1 1.5 2 4 6 10]; 
                smoo_text = 'Smoothing {\itt} (s)';
                grid_text = sprintf('Bin size (mm)');         
            case {'kadaptive'}
                bvnow = [1 3 5 10 40 160 320];
                svnow = [100 1000 6000 10000 16000 24000 32000];  
                smoo_text = sprintf('Smoothing %c',945);
                smoo_spac = -1.40;                
                grid_text = sprintf('Bin size (mm)');  
            case {'fyhn'}
                bvnow = [5 10 20 40 80 160 320];
                svnow = [10 20 30 40 80 160 320];
                smoo_text = sprintf('Smoothing %c (mm)',963);
                grid_text = sprintf('Bin size (mm)');                   
            otherwise
                keyboard
        end

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

        tcut = 16*60;
        ppox = pox(pot<tcut);
        ppoy = poy(pot<tcut);
        ppot = pot(pot<tcut);            
        pspx = spx(spt<tcut);
        pspy = spy(spt<tcut);
        pspt = spt(spt<tcut);            
        pos = [ppox ppoy]; % positions in mm
        spk = [pspx pspy]; % spikes in mm

        loopout = looper(length(bvnow)*length(svnow));        
        for bb = 1:length(bvnow)
            for ss = 1:length(svnow)
                ax = axes('Units','pixels','Position',[bvec(bb) svec(ss) psiz psiz],'FontSize',fsiz);    
                    if bb==1 && ss==1
                        ah = add_panel_title(ls{mm},sprintf(dnames{mm}),'yoffset',-10,'xoffset',0,'width',400,'fontsize',tfonts);
                    end

                    if mem
                        ratemap = maps_for_plotting{mm,bb,ss};
                    else                    
                        rmset = mapset; % rate mapper settings structure - passing mapset directly to ratemapper causes a memory leak
                        rmset.method = mname;
                        rmset.binsize = bvnow(bb);
                        rmset.ssigma = svnow(ss);
                        rmset.ash = svnow(ss);
                        rmset.maplims = [-max(abs(epoly(:,1))) -max(abs(epoly(:,2))) max(abs(epoly(:,1))) max(abs(epoly(:,2)))]; 
                        rmset.pot = ppot;
                        rmset.spt = pspt;
                        if strcmp(mname,'fyhn')
                            rmset.maxdist = 320;
                            rmset.mindist = 50;                                
                        end
                        [ratemap,dwellmap,spikemap,rmset,speedlift] = rate_mapper(pos,spk,rmset);
                        maps_for_plotting{mm,bb,ss} = ratemap;
                    end

                    im = imagesc(ratemap,'alphadata',~isnan(ratemap)); hold on;
                    daspect([1 1 1])
                    axis xy off tight
                    colormap(gca,rmap_colormap)     
                    ax.CLim = [0 prctile(ratemap(:),99)];

                    if bb==1 && ss==floor(length(svnow)/2)
                        % if mm==1 || mm==4
                            text(smoo_spac,-0.5,smoo_text,'HorizontalAl','center','Units','normalized','FontSize',fsiz,'rotation',90)
                        % end
                    elseif bb==floor(length(bvnow)/2) && ss==length(svnow)
                        % if mm>3
                            text(1.5,-0.52,grid_text,'HorizontalAl','center','Units','normalized','FontSize',fsiz)    
                        % end
                    end
                    if ss==length(svnow)
                        text(0.5,-0.16,sprintf('%s',num2str(bvnow(bb))),'HorizontalAl','center','Units','normalized','FontSize',fsiz,'rotation',0)   
                    end
                    if bb==1
                        text(-0.06,0.5,sprintf('%s',num2str(svnow(ss))),'HorizontalAl','right','Units','normalized','FontSize',fsiz,'rotation',0)                             
                    end                    
                    loopout = looper(loopout);        
            end
        end
% return
    end
    if ~mem
        save(fname,'maps_for_plotting');
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ################################################################# %% Save the figure
    % Save the figure    
    disp(sprintf('\tSaving figure...'))    
    figname = [fig_dir2 '\Fig 6.png'];
    [~,~,~] = mkdir(fig_dir2);    
    exportgraphics(fig1,figname,'BackgroundColor','w','ContentType','image','Resolution',350);  
    close(fig1)    
    
    


































