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
        sdir = [scratch_space '\' config.environs{1} '_sigma8000_duration64_sdata.mat'];            
        disp(sprintf('\t\t...loading %s',sdir));            
        load(sdir,'sdata'); 
    end

    ename = 'arena120cm'; 
    mname = 'ksde';
    fmap_colormap = 'turbo';
    other_colormap = 'hot';
    
%% ##################################### Heading 2
%% #################### Heading 3
%% Heading 4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ################################################################# %% FUNCTION BODY
    % Figure settings
    fig1 = figure('visible',fig_vis,'Position',[50,60,1400,800]); 
    set(gcf,'InvertHardCopy','off'); % this stops white lines being plotted black      
    set(gcf,'color','w'); % makes the background colour white
    colormap(jet(256)); % to make sure the colormap is not the horrible default one

    binsize = 50;
    
%% #################### Spikes and trajectory
    xnow = 50;
    ynow = 500;
    excell = 12;

    ax1 = axes('Units','pixels','Position',[150 ynow 200 200]);    
        %ah = add_panel_title('a',sprintf('Random walk & simulated spikes'),'yoffset',5,'xoffset',0,'width',400);  

        pindx = sdata.pos_index(excell);
        pos_now = all_pos{pindx,1};
        pox = pos_now(:,1);
        poy = pos_now(:,2);
        pot = pos_now(:,3);
        spk = sdata.spk{excell};
        spx = spk(:,1);
        spy = spk(:,2);
        spt = spk(:,3);
        pox = pox(pot<16*60);
        poy = poy(pot<16*60);
        spx = spx(spt<16*60);
        spy = spy(spt<16*60);

        p1 = plot(pox,poy,'Color',[.5 .5 .5 .5]); hold on;
        s1 = scatter(spx,spy,12,'r','filled','MarkerEdgeColor','none','MarkerFaceAlpha',0.75);
        plot(epoly(:,1),epoly(:,2),'k','LineWidth',1.5);
        
        daspect([1 1 1])
        axis xy off  
        text_offset = 1.07;
        text(0,text_offset,'16 minute session','Units','normalized','FontSize',8,'HorizontalAl','left');

        [~,leg] = legendflex([p1 s1]...
            ,{'positions','spikes'}...
            ,'anchor',{'s','s'},'ncol',1,'box','off','buffer',[-30,-40],'xscale',1,'fontsize',9); 

    % create ratemap
    rmset = mapset; % rate mapper settings structure - passing mapset directly to ratemapper causes a memory leak
    rmset.method = mname;
    rmset.binsize = 40;
    rmset.ssigma = 100;
    rmset.maplims = [-max(abs(epoly(:,1))) -max(abs(epoly(:,2))) max(abs(epoly(:,1))) max(abs(epoly(:,2)))];                     
    [ratemap,dwellmap,spikemap,rmset,~] = rate_mapper([pox poy],[spx spy],rmset);    
    
    % trajectory + spikes + bins
    vbuff = 20;    
    ax2 = axes('Units','pixels','Position',[80 ax1.Position(2)-240 140 140]);
        map_pos = rmset.points_to_map([pox poy]);
        map_epoly = rmset.points_to_map(epoly); 
        map_grid = rmset.points_to_map([rmset.xgrid(:) rmset.ygrid(:)]);
        map_spk = rmset.points_to_map([spx spy]);    
          
        p1 = plot(map_pos(:,1),map_pos(:,2),'Color',[.5 .5 .5 .5]); hold on;
        plot(map_epoly(:,1),map_epoly(:,2),'k','LineWidth',1.5); hold on;
        s1 = scatter(map_spk(:,1),map_spk(:,2),12,'r','filled','MarkerEdgeColor','none','MarkerFaceAlpha',0.25); hold on;        
        daspect([1 1 1])
        axis xy off  

        [Xq,Yq] = meshgrid(movmean(map_grid(:,1),2,'EndPoints','discard'),movmean(map_grid(:,2),2,'EndPoints','discard'));
        scatter(Xq,Yq,5,'b','filled','MarkerEdgeColor','none','MarkerFaceAlpha',0.5); hold on;

        text(0,text_offset,sprintf('All query points'),'Units','normalized','FontSize',8,'HorizontalAl','left')
        ax2.XLim = [min(map_grid(:,1)) max(map_grid(:,1))];
        ax2.YLim = [min(map_grid(:,2)) max(map_grid(:,2))];  
        
    % trajectory + spikes + bins + radius
    ax6 = axes('Units','pixels','Position',[ax2.Position(1)+ax2.Position(3)+vbuff ax2.Position(2) ax2.Position(3) ax2.Position(4)]);  
        p1 = plot(map_pos(:,1),map_pos(:,2),'Color',[.5 .5 .5 .25]); hold on;
        plot(map_epoly(:,1),map_epoly(:,2),'Color','k','LineWidth',1.5); hold on;
        s1 = scatter(map_spk(:,1),map_spk(:,2),12,'r','filled','MarkerEdgeColor','none','MarkerFaceAlpha',0.25); hold on;
        scatter(Xq,Yq,5,'b','filled','MarkerEdgeColor','none','MarkerFaceAlpha',0.25); hold on;        
        daspect([1 1 1])
        axis xy off  

        text(0,text_offset,sprintf('Example query point'),'Units','normalized','FontSize',8,'HorizontalAl','left')
        ax6.XLim = ax2.XLim;
        ax6.YLim = ax2.YLim;

        idx = find(Xq(:)==13 & Yq(:)==11,1);
        exbin = idx;
        scatter(Xq(exbin),Yq(exbin),8,'b','filled','MarkerEdgeColor','none','MarkerFaceAlpha',1); hold on;
        plot(Xq(exbin),Yq(exbin),'ko','MarkerSize',10)

    % weighted distances to all positions
    ax7 = axes('Units','pixels','Position',[ax6.Position(1)+ax6.Position(3)+vbuff+20 ax6.Position(2) ax6.Position(3) ax6.Position(4)]);         
        % prepare smoothing coefficients
        coeff1 = (sqrt(2*pi) .* rmset.ssigma);
        coeff2 = (1 ./ (sqrt(2*pi) .* rmset.ssigma));
        k = @(x) ( (exp(-0.5 * (x./rmset.ssigma).^2) ./ coeff1) ./ coeff2 ); % gaussian function        
        
        pos_now = map_pos .* rmset.binsize;
        bin_now = [Xq(exbin),Yq(exbin)] .* rmset.binsize;
        
        rindx = true(size(pos_now,1),1);
        dp = sqrt(sum((pos_now(rindx,:)-bin_now).^2,2)); % calculate the distance to every position data point
        dp_norm = k(dp(:)); % gaussian weight the distance values
        
        scatter(map_pos(:,1),map_pos(:,2),7,dp_norm,'filled','MarkerEdgeColor','none','MarkerFaceAlpha',1); hold on;
        plot(map_epoly(:,1),map_epoly(:,2),'Color','k','LineWidth',1.5); hold on;
        daspect([1 1 1])
        axis xy off          
        ax7.XLim = ax2.XLim;
        ax7.YLim = ax2.YLim;    
        scatter(Xq(exbin),Yq(exbin),8,'b','filled','MarkerEdgeColor','none','MarkerFaceAlpha',1); hold on;
        plot(Xq(exbin),Yq(exbin),'ko','MarkerSize',10)        
        text(0,1.12,sprintf('Kernel weighted distance\nto all positions'),'Units','normalized','FontSize',8,'HorizontalAl','left')
        
        text(1.25,0.5,sprintf('For every\nquery point\n%c(distances)',931),'HorizontalAl','center','FontSize',8,'Units','normalized')        
        
    % weighted distances to all spikes
    ax8 = axes('Units','pixels','Position',[ax6.Position(1)+ax6.Position(3)+vbuff+20 ax6.Position(2)-200 ax6.Position(3) ax6.Position(4)]);         
        pos_now = map_spk .* rmset.binsize;
        bin_now = [Xq(exbin),Yq(exbin)] .* rmset.binsize;
        
        rindx = true(size(pos_now,1),1);
        dp = sqrt(sum((pos_now(rindx,:)-bin_now).^2,2)); % calculate the distance to every position data point
        dp_norm = k(dp(:)); % gaussian weight the distance values
        
        scatter(map_spk(:,1),map_spk(:,2),10,dp_norm,'filled','MarkerEdgeColor','none','MarkerFaceAlpha',1); hold on;
        plot(map_epoly(:,1),map_epoly(:,2),'Color','k','LineWidth',1.5); hold on;
        daspect([1 1 1])
        axis xy off          
        ax8.XLim = ax2.XLim;
        ax8.YLim = ax2.YLim;
        scatter(Xq(exbin),Yq(exbin),8,'b','filled','MarkerEdgeColor','none','MarkerFaceAlpha',1); hold on;
        plot(Xq(exbin),Yq(exbin),'ko','MarkerSize',10)        
        text(0,1.12,sprintf('Kernel weighted distance\nto all spikes'),'Units','normalized','FontSize',8,'HorizontalAl','left')
              
        text(1.25,0.5,sprintf('For every\nquery point\n%c(distances)',931),'HorizontalAl','center','FontSize',8,'Units','normalized')        
        
   % dwellmap
    ax10 = axes('Units','pixels','Position',[ax7.Position(1)+ax7.Position(3)+vbuff+55 ax7.Position(2) ax7.Position(3) ax7.Position(4)]);  
        imagesc(dwellmap,'alphadata',~isnan(dwellmap)); hold on;
        %plot(map_epoly(:,1),map_epoly(:,2),'k','LineWidth',1.5); hold on;
        daspect([1 1 1])
        axis xy off  
        view(0,90)
        ax10.XLim = ax2.XLim;
        ax10.YLim = ax2.YLim;
        text(0,text_offset,sprintf('KSDE of positions'),'Units','normalized','FontSize',8,'HorizontalAl','left')
        colormap(ax10,other_colormap)
        
        axc = axes('Units','pixels','Position',[ax10.Position(1)+ax10.Position(3)+10 ax10.Position(2)+0 12 ax10.Position(4)*0.9]); 
            mat = linspace(ax10.CLim(1),ax10.CLim(2),100)';
            imagesc(ones(size(mat)),mat,mat);
            colormap(axc,ax10.Colormap);
            axis xy

            axc.YTick = ax10.CLim;
            axc.YTickLabel = {'0','Max'};
            axc.XTick = [];
            axc.YAxisLocation = 'right';
        
    % spikemap
    ax11 = axes('Units','pixels','Position',[ax8.Position(1)+ax8.Position(3)+vbuff+55 ax8.Position(2) ax8.Position(3) ax8.Position(4)]);    
        imagesc(spikemap,'alphadata',~isnan(spikemap)); hold on;
        %plot(map_epoly(:,1),map_epoly(:,2),'k','LineWidth',1.5); hold on;
        daspect([1 1 1])
        axis xy off  
        view(0,90)
        ax11.XLim = ax2.XLim;
        ax11.YLim = ax2.YLim;    
        text(0,text_offset,sprintf('KSDE of spikes'),'Units','normalized','FontSize',8,'HorizontalAl','left')
        colormap(ax11,other_colormap)

    % ratemap
    ax12 = axes('Units','pixels','Position',[ax1.Position(1)+ax1.Position(3)+90 ax1.Position(2) ax1.Position(3) ax1.Position(4)]);    
        imagesc(ratemap,'alphadata',~isnan(ratemap)); hold on;
        daspect([1 1 1])
        axis xy off  
        view(0,90)
        ax12.XLim = ax2.XLim;
        ax12.YLim = ax2.YLim;     
        ax12.CLim(1) = 0;        
        text(0,text_offset,sprintf('Firing rate map'),'Units','normalized','FontSize',8,'HorizontalAl','left')
        colormap(ax12,fmap_colormap)
    
        axc = axes('Units','pixels','Position',[ax12.Position(1)+ax12.Position(3)+10 ax12.Position(2)+0 12 ax12.Position(4)*0.9]); 
            mat = linspace(ax12.CLim(1),ax12.CLim(2),100)';
            imagesc(ones(size(mat)),mat,mat);
            colormap(axc,fmap_colormap);
            axis xy

            axc.YTick = ax12.CLim;
            axc.YTickLabel = {'0 Hz','Max'};
            axc.XTick = [];
            axc.YAxisLocation = 'right';
            text(0,1.12,sprintf('Firing\nrate (Hz)'),'FontSize',8,'HorizontalAl','left','Units','normalized')    

        lwid = 1.3;
        annotation('line',[0.042 0.10],[0.73 0.73],'LineWidth',lwid) % to position/spike plot
        
        pos_y = 0.510;
        annotation('line',[0.042 0.042],[0.73 pos_y],'LineWidth',lwid) % left side
        annotation('arrow',[0.042 0.054],[pos_y pos_y],'LineWidth',lwid); % to pos query points 
        annotation('arrow',[0.115 0.166],[pos_y pos_y],'LineWidth',lwid) % pos query points to example radii
        annotation('arrow',[0.248 0.2935],[pos_y pos_y],'LineWidth',lwid); % example query point to pos dists
        annotation('line',[0.279 0.279],[pos_y 0.18],'LineWidth',lwid); % line down to spike dists
        annotation('arrow',[0.279 0.2935],[0.18 0.18],'LineWidth',lwid); % example query point to pos dists
        annotation('arrow',[0.36 0.445],[pos_y pos_y],'LineWidth',lwid) % pos dists to dwellmap
        annotation('arrow',[0.36 0.445],[0.26 0.26],'LineWidth',lwid) % spike dists to spikemap
        annotation('line',[0.518 0.59],[pos_y pos_y],'LineWidth',lwid) % dwellmap to right side
        annotation('line',[0.518 0.59],[0.26 0.26],'LineWidth',lwid) % spikemap to right side
        annotation('line',[0.59 0.59],[0.26 0.73],'LineWidth',lwid) % right side                
        annotation('arrow',[0.59 0.48],[0.73 0.73],'LineWidth',lwid) % right side to ratemap

%         delete(an2)
        an2 = annotation('textbox','String',sprintf('For each query point:\nrate = spike KSDE / pos KSDE\n\\times sample interval'),'Position',[0.54 0.52 0.1 0.1],'EdgeColor','none','BackgroundColor','w','HorizontalAl','center');         
                
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ################################################################# %% Save the figure
    % Save the figure    
    disp(sprintf('\tSaving figure...'))    
    figname = [fig_dir2 '\' mname '_fig_meth.png'];
    [~,~,~] = mkdir(fig_dir2);    
    exportgraphics(gcf,figname,'BackgroundColor','w','ContentType','image','Resolution',250);  
    close(gcf)    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    