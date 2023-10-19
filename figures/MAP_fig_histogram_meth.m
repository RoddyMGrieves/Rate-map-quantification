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
    
    % if ~exist('sdata','var')
        sdir = [scratch_space '\' config.environs{1} '_sigma8000_duration64_sdata.mat'];            
        disp(sprintf('\t\t...loading %s',sdir));            
        load(sdir,'sdata'); 
    % end

    ename = 'arena120cm'; 
    mname = 'histogram';
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

    % create unsmoothed ratemap
    rmset = mapset; % rate mapper settings structure - passing mapset directly to ratemapper causes a memory leak
    rmset.method = mname;
    rmset.binsize = 50;
    rmset.ssigma = 50;
    rmset.smethod = 3; % no smoothing
    rmset.maplims = [-max(abs(epoly(:,1))) -max(abs(epoly(:,2))) max(abs(epoly(:,1))) max(abs(epoly(:,2)))]; 
    rmset.padding = 0;
    [ratemap,dwellmap,spikemap,rmset,~] = rate_mapper([pox poy],[spx spy],rmset);    
    
    % trajectory
    ax2 = axes('Units','pixels','Position',[90 ax1.Position(2)-240 140 140]);
        map_pos = rmset.points_to_map([pox poy]);
        map_epoly = rmset.points_to_map(epoly); 
        map_grid = rmset.points_to_map([rmset.xgrid(:) rmset.ygrid(:)]);            
        p1 = plot(map_pos(:,1),map_pos(:,2),'Color',[.5 .5 .5 .5]); hold on;
        plot(map_epoly(:,1),map_epoly(:,2),'k','LineWidth',1.5); hold on;
        daspect([1 1 1])
        axis xy off  
        ax2.XLim = [min(map_grid(:,1)) max(map_grid(:,1))];
        ax2.YLim = [min(map_grid(:,2)) max(map_grid(:,2))];        
        text(0,text_offset,sprintf('Position data'),'Units','normalized','FontSize',8,'HorizontalAl','left')
    
    % spikes
    ax3 = axes('Units','pixels','Position',[ax2.Position(1) ax2.Position(2)-200 ax2.Position(3) ax2.Position(4)]);   
        map_spk = rmset.points_to_map([spx spy]);    
        s1 = scatter(map_spk(:,1),map_spk(:,2),12,'r','filled','MarkerEdgeColor','none','MarkerFaceAlpha',0.75); hold on;
        plot(map_epoly(:,1),map_epoly(:,2),'k','LineWidth',1.5); hold on;
        daspect([1 1 1])
        axis xy off 
        ax3.XLim = [min(map_grid(:,1)) max(map_grid(:,1))];
        ax3.YLim = [min(map_grid(:,2)) max(map_grid(:,2))];         
        text(0,text_offset,sprintf('Spikes'),'Units','normalized','FontSize',8,'HorizontalAl','left')
    
    % trajectory + bins
    vbuff = 10;
    ax4 = axes('Units','pixels','Position',[ax2.Position(1)+ax2.Position(3)+vbuff ax2.Position(2)-0 ax2.Position(3) ax2.Position(4)]);  
        p1 = plot(map_pos(:,1),map_pos(:,2),'Color',[.5 .5 .5 .5]); hold on;
        plot(map_epoly(:,1),map_epoly(:,2),'k','LineWidth',1.5); hold on;
        daspect([1 1 1])
        axis xy off  
    
        [Xq,Yq] = meshgrid(movmean(map_grid(:,1),2,'EndPoints','discard'),movmean(map_grid(:,2),2,'EndPoints','discard'));
        Vq = ones(size(Xq));
        surf(Xq,Yq,Vq,'EdgeColor','k','FaceAlpha',0); hold on;
        view(0,90)
        text(0,text_offset,sprintf('Position bins'),'Units','normalized','FontSize',8,'HorizontalAl','left')
        ax4.XLim = ax2.XLim;
        ax4.YLim = ax2.YLim;
        
    % spikes + bins
    ax5 = axes('Units','pixels','Position',[ax3.Position(1)+ax3.Position(3)+vbuff ax3.Position(2)-0 ax3.Position(3) ax3.Position(4)]);    
        s1 = scatter(map_spk(:,1),map_spk(:,2),12,'r','filled','MarkerEdgeColor','none','MarkerFaceAlpha',0.75); hold on;
        plot(map_epoly(:,1),map_epoly(:,2),'k','LineWidth',1.5); hold on;
        daspect([1 1 1])
        axis xy off         
        surf(Xq,Yq,Vq,'EdgeColor','k','FaceAlpha',0); hold on;
        view(0,90)
        text(0,text_offset,sprintf('Spike bins'),'Units','normalized','FontSize',8,'HorizontalAl','left')
        ax5.XLim = ax3.XLim;
        ax5.YLim = ax3.YLim;
        
    % dwellmap
    ax6 = axes('Units','pixels','Position',[ax4.Position(1)+ax4.Position(3)+vbuff ax4.Position(2) ax4.Position(3) ax4.Position(4)]);  
        imagesc(dwellmap,'alphadata',~isnan(dwellmap)); hold on;
        %plot(map_epoly(:,1),map_epoly(:,2),'k','LineWidth',1.5); hold on;
        daspect([1 1 1])
        axis xy off  
        view(0,90)
        ax6.XLim = ax4.XLim;
        ax6.YLim = ax4.YLim;
        text(0,text_offset,sprintf('Dwell map'),'Units','normalized','FontSize',8,'HorizontalAl','left')
        colormap(ax6,other_colormap)
        
    % spikemap
    ax7 = axes('Units','pixels','Position',[ax5.Position(1)+ax5.Position(3)+vbuff ax5.Position(2) ax5.Position(3) ax5.Position(4)]);    
        imagesc(spikemap,'alphadata',~isnan(spikemap)); hold on;
        %plot(map_epoly(:,1),map_epoly(:,2),'k','LineWidth',1.5); hold on;
        daspect([1 1 1])
        axis xy off  
        view(0,90)
        ax7.XLim = ax5.XLim;
        ax7.YLim = ax5.YLim;
        text(0,text_offset,sprintf('Spike map'),'Units','normalized','FontSize',8,'HorizontalAl','left')
        colormap(ax7,other_colormap)
        
    % dwellmap convolution
    ax8 = axes('Units','pixels','Position',[ax6.Position(1)+ax6.Position(3)+15 ax6.Position(2)-(ax6.Position(4)*0.4) ax6.Position(3) ax6.Position(4)]);
        H = fspecial('gaussian',rmset.ssize,rmset.ssigma);
        imagesc(H); hold on;
        daspect([1 1 1])
        axis xy off
        view(90,90)
        ax8.XLim = ax6.XLim;
        ax8.YLim = ax6.YLim;
        text(-0.11,0.92,sprintf('{\\times}'),'Units','normalized','FontSize',15,'HorizontalAl','left')
        text(0,1.1,sprintf('Gaussian\nkernel'),'Units','normalized','FontSize',8,'HorizontalAl','left')
        colormap(ax8,other_colormap)
    
    % spikemap convolution
    ax9 = axes('Units','pixels','Position',[ax7.Position(1)+ax7.Position(3)+15 ax7.Position(2)-(ax7.Position(4)*0.4) ax7.Position(3) ax7.Position(4)]);           
        imagesc(H); hold on;
        daspect([1 1 1])
        axis xy off
        view(90,90)
        ax9.XLim = ax7.XLim;
        ax9.YLim = ax7.YLim;    
        text(-0.11,0.92,sprintf('{\\times}'),'Units','normalized','FontSize',15,'HorizontalAl','left')
        text(0,1.1,sprintf('Gaussian\nkernel'),'Units','normalized','FontSize',8,'HorizontalAl','left')
        colormap(ax9,other_colormap)
    
    % create smoothed ratemap
    rmset = mapset; % rate mapper settings structure - passing mapset directly to ratemapper causes a memory leak
    rmset.method = mname;
    rmset.binsize = 50;
    rmset.ssigma = 50;
    rmset.smethod = 1; % apply smoothing
    rmset.maplims = [-max(abs(epoly(:,1))) -max(abs(epoly(:,2))) max(abs(epoly(:,1))) max(abs(epoly(:,2)))];                     
    [ratemap,dwellmap,spikemap,rmset,~] = rate_mapper([pox poy],[spx spy],rmset);     
    
    % dwellmap
    ax10 = axes('Units','pixels','Position',[ax6.Position(1)+ax6.Position(3)+80 ax6.Position(2) ax6.Position(3) ax6.Position(4)]);  
        imagesc(dwellmap,'alphadata',~isnan(dwellmap)); hold on;
        %plot(map_epoly(:,1),map_epoly(:,2),'k','LineWidth',1.5); hold on;
        daspect([1 1 1])
        axis xy off  
        view(0,90)
        ax10.XLim = ax4.XLim;
        ax10.YLim = ax4.YLim;
        text(0,text_offset,sprintf('Smoothed dwell map'),'Units','normalized','FontSize',8,'HorizontalAl','left')
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
    ax11 = axes('Units','pixels','Position',[ax7.Position(1)+ax7.Position(3)+80 ax7.Position(2) ax7.Position(3) ax7.Position(4)]);    
        imagesc(spikemap,'alphadata',~isnan(spikemap)); hold on;
        %plot(map_epoly(:,1),map_epoly(:,2),'k','LineWidth',1.5); hold on;
        daspect([1 1 1])
        axis xy off  
        view(0,90)
        ax11.XLim = ax5.XLim;
        ax11.YLim = ax5.YLim;    
        text(0,text_offset,sprintf('Smoothed spike map'),'Units','normalized','FontSize',8,'HorizontalAl','left')
        colormap(ax11,other_colormap)
    
    % ratemap
    ax12 = axes('Units','pixels','Position',[ax1.Position(1)+ax1.Position(3)+90 ax1.Position(2) ax1.Position(3) ax1.Position(4)]);    
        imagesc(ratemap,'alphadata',~isnan(ratemap)); hold on;
        %plot(map_epoly(:,1),map_epoly(:,2),'k','LineWidth',1.5); hold on;
        daspect([1 1 1])
        axis xy off  
        view(0,90)
        ax12.XLim = ax5.XLim;
        ax12.YLim = ax5.YLim;     
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
        annotation('line',[0.03 0.10],[0.73 0.73],'LineWidth',lwid) % to position/spike plot
        
        spike_y = 0.260;        
        annotation('line',[0.03 0.03],[0.73 spike_y],'LineWidth',lwid) % down side
        annotation('arrow',[0.03 0.06],[spike_y spike_y],'LineWidth',lwid) % to spike plot 
        annotation('arrow',[0.0925 0.164],[spike_y spike_y],'LineWidth',lwid) % spikes to spike bins
        annotation('arrow',[0.215 0.275],[spike_y spike_y],'LineWidth',lwid) % spike bins to spike map
        annotation('arrow',[0.318 0.432],[spike_y spike_y],'LineWidth',lwid) % spike map to smoothed spike map
        annotation('line',[0.513 0.58],[spike_y spike_y],'LineWidth',lwid) % smoothed spike map to side
        annotation('line',[0.58 0.58],[spike_y 0.73],'LineWidth',lwid) % right side
        annotation('arrow',[0.58 0.48],[0.73 0.73],'LineWidth',lwid) % right side to ratemap
    
        pos_y = 0.510;
        annotation('arrow',[0.03 0.06],[pos_y pos_y],'LineWidth',lwid) % to pos plot 
        annotation('arrow',[0.115 0.164],[pos_y pos_y],'LineWidth',lwid) % pos to pos bins
        annotation('arrow',[0.217 0.275],[pos_y pos_y],'LineWidth',lwid) % pos bins to pos map
        annotation('arrow',[0.318 0.432],[pos_y pos_y],'LineWidth',lwid) % pos map to smoothed pos map
        annotation('line',[0.513 0.58],[pos_y pos_y],'LineWidth',lwid) % smoothed pos map to side    

        % an2 = annotation('textbox','String',sprintf('rate =\nspike map / dwell map\n\\times sample interval'),'Position',[0.52 0.60 0.12 0.078],'EdgeColor','none','BackgroundColor','w','HorizontalAl','center'); 
        an2 = annotation('textbox','String',sprintf('rate =\nspike map / dwell map'),'Position',[0.52 0.60 0.12 0.06],'EdgeColor','none','BackgroundColor','w','HorizontalAl','center'); 
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ################################################################# %% Save the figure
    % Save the figure    
    disp(sprintf('\tSaving figure...'))    
    figname = [fig_dir2 '\' mname '_fig_meth.png'];
    [~,~,~] = mkdir(fig_dir2);    
    exportgraphics(gcf,figname,'BackgroundColor','w','ContentType','image','Resolution',250);  
    close(gcf)        
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    