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
    
%     if ~exist('sdata','var')
%         sdir = [scratch_space '\' config.environs{1} '_sigma16000_duration64_sdata.mat'];            
%         disp(sprintf('\t\t...loading %s',sdir));            
%         load(sdir,'sdata'); 
%     end

    ename = 'arena120cm'; 
    override_maps = 0;
    rmap_colormap = 'turbo';
    
%% ##################################### Heading 2
%% #################### Heading 3
%% Heading 4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ################################################################# %% FUNCTION BODY
    % Figure settings
    fig1 = figure('visible','on','Position',[50,60,1400,900]); 
    set(gcf,'InvertHardCopy','off'); % this stops white lines being plotted black      
    set(gcf,'color','w'); % makes the background colour white
    colormap(jet(256)); % to make sure the colormap is not the horrible default one

    excell = [70 3 33 86 28 22];            
    switch mname
        case {'histogram'}
            smoo_text = sprintf('Smoothing %c (mm)',963);
            grid_text = sprintf('Bin size (mm)');
            excell = [45 13 86 114 124 128];            
        case {'ash'}
            smoo_text = 'Smoothing {\itm}';
            grid_text = sprintf('Bin size (mm)');      
            excell = [67 27 59 92 113 132];                        
            
        case {'ksde'}
            smoo_text = sprintf('Smoothing bandwidth (mm)');
            grid_text = sprintf('Bin size (mm)');             
            excell = [70 3 33 86 28 22];  
            
        case {'kyadaptive'}
            smoo_text = 'Smoothing {\itt}';
            grid_text = sprintf('Bin size (mm)');              
            excell = [90 85 41 86 128 38];            

        case {'kadaptive'} 
            smoo_text = sprintf('Smoothing %c',945);
            grid_text = sprintf('Bin size (mm)');  
            excell = [70 3 33 86 28 22];  
            
        case {'fyhn'} 
            smoo_text = sprintf('Smoothing %c (mm)',963);
            grid_text = sprintf('Bin size (mm)');   
            excell = [70 3 33 86 28 22];              
    end
    
%% #################### Example cells
    xnow = 70;
    ynow = 550;
        
    if 1
        ncells = length(excell);
        
        maps_for_plotting = cell(ncells,2);
        fname = [fig_dir mname '_ex_maps.mat'];
        if exist(fname,'file') && ~override_maps
            mem = 1;
            disp(sprintf('\t\t...loading old maps'))
            load(fname,'maps_for_plotting')
        else
            disp(sprintf('\t\t...generating maps'))            
            mem = 0;
        end
        
        sigma_now = '4000';
        time_now = '16';                                
        load([scratch_space '\' ename '_sigma' num2str(sigma_now) '_duration64_sdata.mat'])            

        wdir = [scratch_space '\' ename '_64_walk1.mat'];
        load(wdir,'epoly');  
        tcut = 16*60;

        load([scratch_space '\' ename '_' mname '_sigma' num2str(sigma_now) '_duration' num2str(time_now) '_gamultiobj.mat'],'balance_solution1','minerr_solution')        
        
        msiz = 100;
        ybuff = msiz+40;      
        xbuff = msiz+5;

        for cc = 1:ncells
            % true PDF for the cell at 1mm resolution
            xtemp = xnow+(xbuff*(cc-1));
            axpos = [xtemp ynow msiz msiz];
            ax = axes('Units','pixels','Position',axpos);        
                pmap = sdata.pmap_1mm{excell(cc),1};
                im = imagesc(pmap); hold on;
                daspect([1 1 1])
                axis xy off tight
                colormap(gca,rmap_colormap)  

                text(0.5,1.1,sprintf('cell %d',cc),'HorizontalAl','center','Units','normalized','FontSize',8);  
                if cc==1
                    text(-0.1,1.30,'Spike probability map','HorizontalAl','left','Units','normalized','FontSize',10,'rotation',0);
                end                
                
                if cc==ncells
                    axc = axes('Units','pixels','Position',[ax.Position(1)+ax.Position(3)+20 ax.Position(2) 12 ax.Position(4)-10]);
                        mat = (linspace(0,100,100))';
                        imagesc(mat,ones(size(mat)),mat);
                        colormap(axc,rmap_colormap);
                        axis xy

                        axc.YTick = [];
                        axc.XTick = [];
                        text(0.5,1.30,sprintf('Spike\nprobability'),'FontSize',8,'HorizontalAl','center','Units','normalized')
                        text(0.5,1.1,sprintf('Max'),'FontSize',8,'HorizontalAl','center','Units','normalized')
                        text(0.5,-0.1,sprintf('0'),'FontSize',8,'HorizontalAl','center','Units','normalized')  
                end                
                
            % spikes and position data
            axpos = [xtemp ynow-ybuff msiz msiz];
            ax = axes('Units','pixels','Position',axpos);              
                pindx = sdata.pos_index(excell(cc));
                pos_now = all_pos{pindx,1};
                pox = pos_now(:,1);
                poy = pos_now(:,2);
                pot = pos_now(:,3);                 
                spk = sdata.spk{excell(cc)};
                spx = spk(:,1);
                spy = spk(:,2);
                spt = spk(:,3);
                
                plot(pox,poy,'Color',[.5 .5 .5 .5]); hold on;
                scatter(spx,spy,5,'r','filled','MarkerEdgeColor','none','MarkerFaceAlpha',0.5);
                daspect([1 1 1])
                axis xy off  
                
                if cc==1
                    text(-0.1,1.15,'Positions and spikes','HorizontalAl','left','Units','normalized','FontSize',10,'rotation',0);
                end
                
            % pareto optimal solution map
            axpos = [xtemp ynow-ybuff*2 msiz msiz];                
            ax = axes('Units','pixels','Position',axpos); 
                pindx = sdata.pos_index(excell(cc));
                pos_now = all_pos{pindx,1};
                pox = pos_now(:,1);
                poy = pos_now(:,2);
                pot = pos_now(:,3);
                ppox = pox(pot<tcut);
                ppoy = poy(pot<tcut);
                ppot = pot(pot<tcut);
                
                spk = sdata.spk{excell(cc),1};
                spx = spk(:,1);
                spy = spk(:,2);
                spt = spk(:,3);                
                pspx = spx(spt<tcut);
                pspy = spy(spt<tcut);
                pspt = spt(spt<tcut);                
                pos = [ppox ppoy]; % positions in mm
                spk = [pspx pspy]; % spikes in mm
        
                if mem
                    ratemap = maps_for_plotting{cc,1};
                else                    
                    rmset = mapset; % rate mapper settings structure - passing mapset directly to ratemapper causes a memory leak
                    rmset.method = mname;
                    rmset.binsize = balance_solution1(1);
                    rmset.ssigma = balance_solution1(2);
                    rmset.ash = balance_solution1(2);
                    rmset.maplims = [-max(abs(epoly(:,1))) -max(abs(epoly(:,2))) max(abs(epoly(:,1))) max(abs(epoly(:,2)))];   
                    rmset.pot = ppot;
                    rmset.spt = pspt;
                    if strcmp(mname,'fyhn')
                        rmset.maxdist = 320;
                        rmset.mindist = 50;                                
                    end                    
                    [ratemap,~,~,~,~] = rate_mapper(pos,spk,rmset);
                    maps_for_plotting{cc,1} = ratemap;
                end                

                im = imagesc(ratemap,'alphadata',~isnan(ratemap)); hold on;
                daspect([1 1 1])
                axis xy off tight
                colormap(gca,rmap_colormap) 
                if cc==1
                    text(-0.1,1.15,sprintf('Balanced solution %s = %.1f, %s = %.1f',grid_text,balance_solution1(1),smoo_text,balance_solution1(2)),'HorizontalAl','left','Units','normalized','FontSize',10,'rotation',0);
                end

                if cc==ncells
                    axc = axes('Units','pixels','Position',[ax.Position(1)+ax.Position(3)+20 ax.Position(2) 12 ax.Position(4)-10]);
                        mat = (linspace(0,100,100))';
                        imagesc(mat,ones(size(mat)),mat);
                        colormap(axc,rmap_colormap);
                        axis xy

                        axc.YTick = [];
                        axc.XTick = [];
                        text(0.5,1.30,sprintf('Firing\nrate (Hz)'),'FontSize',8,'HorizontalAl','center','Units','normalized')
                        text(0.5,1.1,sprintf('Max'),'FontSize',8,'HorizontalAl','center','Units','normalized')
                        text(0.5,-0.1,sprintf('0'),'FontSize',8,'HorizontalAl','center','Units','normalized')  
                end
                
            % minimum error solution map
            axpos = [xtemp ynow-ybuff*3 msiz msiz];                
            ax = axes('Units','pixels','Position',axpos);        
                if mem
                    ratemap = maps_for_plotting{cc,2};
                else                    
                    rmset = mapset; % rate mapper settings structure - passing mapset directly to ratemapper causes a memory leak
                    rmset.method = mname;
                    rmset.binsize = minerr_solution(1);           
                    rmset.ssigma = minerr_solution(2);
                    rmset.ash = minerr_solution(2);                    
                    rmset.maplims = [-max(abs(epoly(:,1))) -max(abs(epoly(:,2))) max(abs(epoly(:,1))) max(abs(epoly(:,2)))];  
                    rmset.pot = ppot;
                    rmset.spt = pspt;
                    if strcmp(mname,'fyhn')
                        rmset.maxdist = 320;
                        rmset.mindist = 50;                                
                    end                       
                    [ratemap,~,~,~,~] = rate_mapper(pos,spk,rmset);
                    maps_for_plotting{cc,2} = ratemap;
                end                  
                
                im = imagesc(ratemap,'alphadata',~isnan(ratemap)); hold on;
                daspect([1 1 1])
                axis xy off tight
                colormap(gca,rmap_colormap)  
                if cc==1
                    text(-0.1,1.15,sprintf('Minimum error solution %s = %.1f, %s = %.1f',grid_text,minerr_solution(1),smoo_text,minerr_solution(2)),'HorizontalAl','left','Units','normalized','FontSize',10,'rotation',0);
                end
 
        end
        if ~mem
            save(fname,'maps_for_plotting');
        end        
    end
    
%     set(gcf,'visible','on')
%     return
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ################################################################# %% Save the figure
    % Save the figure    
    disp(sprintf('\tSaving figure...'))    
    figname = [fig_dir2 '\' mname '_fig_examples.png'];
    [~,~,~] = mkdir(fig_dir2);    
    exportgraphics(fig1,figname,'BackgroundColor','w','ContentType','image','Resolution',250);  
    close(fig1)    
    
    
    
    
  
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    