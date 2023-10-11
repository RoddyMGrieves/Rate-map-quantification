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
    mname = 'adaptive';
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
    excell = 4;
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
    rmset.ssigma = 1600;
    rmset.maplims = [-max(abs(epoly(:,1))) -max(abs(epoly(:,2))) max(abs(epoly(:,1))) max(abs(epoly(:,2)))];                     
    [ratemap,dwellmap,spikemap,rmset,~] = rate_mapper([pox poy],[spx spy],rmset);    
    
    % trajectory + spikes + bins
    vbuff = 20;    
    ax2 = axes('Units','pixels','Position',[160 ax1.Position(2)-240 140 140]);
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

        text(0,text_offset,sprintf('Position query points'),'Units','normalized','FontSize',8,'HorizontalAl','left')
        ax2.XLim = [min(map_grid(:,1)) max(map_grid(:,1))];
        ax2.YLim = [min(map_grid(:,2)) max(map_grid(:,2))];  
        
    % trajectory + spikes + bins + radius
    ax6 = axes('Units','pixels','Position',[ax2.Position(1)+ax2.Position(3)+vbuff ax2.Position(2) ax2.Position(3) ax2.Position(4)]);  
        p1 = plot(map_pos(:,1),map_pos(:,2),'Color',[.5 .5 .5 .25]); hold on;
        plot(map_epoly(:,1),map_epoly(:,2),'Color',[1 1 1 .25],'LineWidth',1.5); hold on;
        s1 = scatter(map_spk(:,1),map_spk(:,2),12,'r','filled','MarkerEdgeColor','none','MarkerFaceAlpha',0.25); hold on;
        scatter(Xq,Yq,5,'b','filled','MarkerEdgeColor','none','MarkerFaceAlpha',0.25); hold on;        
        daspect([1 1 1])
        axis xy off  

        text(0,text_offset,sprintf('Example radii'),'Units','normalized','FontSize',8,'HorizontalAl','left')
        ax6.XLim = ax2.XLim;
        ax6.YLim = ax2.YLim;

        idx = find(Xq(:)==8 & Yq(:)==12,1);
        exbin = idx;
        r = rmset.radmap(exbin)/rmset.binsize; % radius in bins
        scatter(Xq(exbin),Yq(exbin),8,'b','filled','MarkerEdgeColor','none','MarkerFaceAlpha',1); hold on;
        rectangle('Position',[Xq(exbin)-r Yq(exbin)-r r*2 r*2],'Curvature',[1 1])
        plot([Xq(exbin) Xq(exbin)+r],[Yq(exbin) Yq(exbin)],'Color','k')
        text(Xq(exbin)+r+1,Yq(exbin),sprintf('%.f mm',rmset.radmap(exbin)),'HorizontalAl','left','FontSize',8)
       
        idx2 = find(Xq(:)==17 & Yq(:)==28,1);        
        exbin2 = idx2;
        r2 = rmset.radmap(exbin2)/rmset.binsize; % radius in bins
        scatter(Xq(exbin2),Yq(exbin2),8,'b','filled','MarkerEdgeColor','none','MarkerFaceAlpha',1); hold on;
        rectangle('Position',[Xq(exbin2)-r2 Yq(exbin2)-r2 r2*2 r2*2],'Curvature',[1 1],'Clipping','on')
        plot([Xq(exbin2) Xq(exbin2)+r2],[Yq(exbin2) Yq(exbin2)],'Color','k')        
        text(Xq(exbin2)+r2+1,Yq(exbin2),sprintf('%.f mm',rmset.radmap(exbin2)),'HorizontalAl','left','FontSize',8)

        text(0.01,-0.15,sprintf('Expand until:\nradius {\\geq} %c / (n {\\times} {\\surd}s)',945),'HorizontalAl','left','FontSize',10,'Units','normalized')

    % distance map
    ax10 = axes('Units','pixels','Position',[ax6.Position(1)+ax6.Position(3)+vbuff+10 ax6.Position(2) ax6.Position(3) ax6.Position(4)]); 
        radmap = rmset.radmap;
        imagesc(radmap,'alphadata',~isnan(radmap)); hold on;
        daspect([1 1 1])
        axis xy off  
        view(0,90)
        ax10.XLim = ax2.XLim;
        ax10.YLim = ax2.YLim;
        ax10.CLim(1) = 0;
        text(0,text_offset,sprintf('All radii (mm)'),'Units','normalized','FontSize',8,'HorizontalAl','left')
        colormap(ax10,other_colormap)
        
        axc = axes('Units','pixels','Position',[ax10.Position(1)+ax10.Position(3)+10 ax10.Position(2)+0 12 ax10.Position(4)*0.9]); 
            mat = linspace(ax10.CLim(1),ax10.CLim(2),100)';
            imagesc(ones(size(mat)),mat,mat);
            colormap(axc,ax10.Colormap);
            axis xy

            axc.YTick = floor(ax10.CLim);
            axc.XTick = [];
            axc.YAxisLocation = 'right';        
            text(0,1.12,sprintf('Radius (mm)'),'FontSize',8,'HorizontalAl','left','Units','normalized')    
        
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
        annotation('line',[0.09 0.10],[0.73 0.73],'LineWidth',lwid) % to position/spike plot
        
        pos_y = 0.510;
        annotation('line',[0.09 0.09],[0.73 pos_y],'LineWidth',lwid) % left side
        annotation('arrow',[0.09 0.11],[pos_y pos_y],'LineWidth',lwid) % to pos query points 
        annotation('arrow',[0.188 0.225],[pos_y pos_y],'LineWidth',lwid) % pos query points to example radii
        annotation('arrow',[0.277 0.345],[pos_y pos_y],'LineWidth',lwid) % example radii to all radii
        annotation('line',[0.40 0.51],[pos_y pos_y],'LineWidth',lwid) % all radii to side  
        annotation('line',[0.51 0.51],[pos_y 0.73],'LineWidth',lwid) % right side                
        annotation('arrow',[0.51 0.48],[0.73 0.73],'LineWidth',lwid) % right side to ratemap

%         delete(an2)
        an2 = annotation('textbox','String',sprintf('For each query point:\nrate = s / n \\times sample interval'),'Position',[0.41 0.54 0.20 0.05],'EdgeColor','none','BackgroundColor','w','HorizontalAl','center');         

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ################################################################# %% Save the figure
    % Save the figure    
    disp(sprintf('\tSaving figure...'))    
    figname = [fig_dir2 '\' mname '_fig_meth.png'];
    [~,~,~] = mkdir(fig_dir2);    
    exportgraphics(gcf,figname,'BackgroundColor','w','ContentType','image','Resolution',250);  
    close(gcf)        
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    