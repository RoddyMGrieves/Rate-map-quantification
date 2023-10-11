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
        sdir = [scratch_space '\' config.environs{1} '_sigma16000_duration64_sdata.mat'];            
        disp(sprintf('\t\t...loading %s',sdir));            
        load(sdir,'sdata'); 
    % end

    ename = 'arena120cm'; 
    mname = 'kadaptive';
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
    rmset.method = 'histogram';
    rmset.binsize = 40;
    rmset.ssigma = 0;
    rmset.maplims = [-max(abs(epoly(:,1))) -max(abs(epoly(:,2))) max(abs(epoly(:,1))) max(abs(epoly(:,2)))];                     
    [ratemap,dwellmap,spikemap,rmset,~] = rate_mapper([pox poy],[spx spy],rmset);    
    
    % trajectory
    ax2 = axes('Units','pixels','Position',[70 ax1.Position(2)-240 140 140]);
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
    
    vbuff = 10;
    % raw dwellmap
    ax6 = axes('Units','pixels','Position',[ax2.Position(1)+ax2.Position(3)+vbuff ax2.Position(2) ax2.Position(3) ax2.Position(4)]);  
        dwellmap = histcounts2(poy,pox,rmset.ygrid,rmset.xgrid);  
        kvect = [1 2 3 4 5];
        dmaps = NaN([size(dwellmap) length(kvect)]);
        for kk = 1:length(kvect)
            SE = strel('disk',kvect(kk));
            kern = double(SE.Neighborhood);
            dmaps(:,:,kk) = imfilter(dwellmap,kern,0,'same');
        end     

        imagesc(dwellmap,'alphadata',~isnan(dwellmap)); hold on;
        daspect([1 1 1])
        axis xy off  
        view(0,90)
        ax6.XLim = ax2.XLim;
        ax6.YLim = ax2.YLim;
        ax6.CLim(1) = 0;                        
        text(0,text_offset,sprintf('Dwell map'),'Units','normalized','FontSize',8,'HorizontalAl','left')
        colormap(ax6,other_colormap)
        
    % raw spikemap
    ax7 = axes('Units','pixels','Position',[ax3.Position(1)+ax3.Position(3)+vbuff ax3.Position(2) ax3.Position(3) ax3.Position(4)]);   
        spikemap = histcounts2(spy,spx,rmset.ygrid,rmset.xgrid);  
        smaps = NaN([size(spikemap) length(kvect)]);
        for kk = 1:length(kvect)
            SE = strel('disk',kvect(kk));
            kern = double(SE.Neighborhood);
            smaps(:,:,kk) = imfilter(spikemap,kern,0,'same');
        end     
    
        imagesc(spikemap,'alphadata',~isnan(spikemap)); hold on;
        daspect([1 1 1])
        axis xy off  
        view(0,90)
        ax7.XLim = ax2.XLim;
        ax7.YLim = ax2.YLim;
        ax7.CLim(1) = 0;                        
        text(0,text_offset,sprintf('Spike map'),'Units','normalized','FontSize',8,'HorizontalAl','left')
        colormap(ax7,other_colormap)
        
    % dwellmap convolution
    ax8 = axes('Units','pixels','Position',[ax6.Position(1)+ax6.Position(3)+15 ax6.Position(2)-(ax6.Position(4)*0.4)+40 200 200],'Color','none'); 
        [Xq,Yq] = meshgrid(rmset.ygrid,rmset.xgrid);
        for zz = 1:length(kvect)
            Xqt = Xq-zz*350;
            Yqt = Yq-zz*400;
            Zqt = -ones(size(Xq))*zz*100;
            Cqt = dmaps(:,:,zz)./max(max(squeeze(dmaps(:,:,zz))));
            surf(Xqt,Yqt,Zqt,Cqt,'EdgeColor','none'); hold on;
            plot3([min(Xqt(:)) min(Xqt(:)) max(Xqt(:)) max(Xqt(:)) min(Xqt(:))],[min(Yqt(:)) max(Yqt(:)) max(Yqt(:)) min(Yqt(:)) min(Yqt(:))],-[zz zz zz zz zz]*100+10,'Color','w','LineWidth',1);            
        end
        daspect([1 1 1])        
        axis xy off
        view(0,90)
        colormap(ax8,other_colormap)

        nbuff = 300;
        arrow([-950 160+nbuff -100],[-2390 -1360+nbuff -490])
        text(-1690-nbuff*1,-560+nbuff*1,-290,sprintf('kernel radius_{i}'),'FontSize',8,'Rotation',45,'HorizontalAl','center')
        
    % spikemap convolution
    ax9 = axes('Units','pixels','Position',[ax7.Position(1)+ax7.Position(3)+15 ax7.Position(2)-(ax7.Position(4)*0.4)+40 200 200],'Color','none');           
        for zz = 1:length(kvect)
            Xqt = Xq-zz*350;
            Yqt = Yq-zz*400;
            Zqt = -ones(size(Xq))*zz*100;
            Cqt = smaps(:,:,zz)./max(max(squeeze(smaps(:,:,zz))));
            surf(Xqt,Yqt,Zqt,Cqt,'EdgeColor','none'); hold on;
            plot3([min(Xqt(:)) min(Xqt(:)) max(Xqt(:)) max(Xqt(:)) min(Xqt(:))],[min(Yqt(:)) max(Yqt(:)) max(Yqt(:)) min(Yqt(:)) min(Yqt(:))],-[zz zz zz zz zz]*100+10,'Color','w','LineWidth',1);            
        end
        daspect([1 1 1])        
        axis xy off
        view(0,90)
        colormap(ax8,other_colormap)
    
        nbuff = 300;
        arrow([-950 160+nbuff -100],[-2390 -1360+nbuff -490])
        text(-1690-nbuff*1,-560+nbuff*1,-290,sprintf('kernel radius_{i}'),'FontSize',8,'Rotation',45,'HorizontalAl','center')
        
        text(0.2,-0.05,sprintf('Find radius_{i} {\\geq} %c / (n_{i} {\\times} {\\surd}s_{i})',945),'HorizontalAl','left','FontSize',10,'Units','normalized')
        
    % distance map
    ax10 = axes('Units','pixels','Position',[ax9.Position(1)+ax9.Position(3)+vbuff+70 ax7.Position(2)+(ax7.Position(3)/2)+20 ax2.Position(3) ax2.Position(4)]); 
        % create ratemap
        rmset = mapset; % rate mapper settings structure - passing mapset directly to ratemapper causes a memory leak
        rmset.method = 'kadaptive';
        rmset.binsize = 40;
        rmset.ssigma = 1600;
        rmset.maplims = [-max(abs(epoly(:,1))) -max(abs(epoly(:,2))) max(abs(epoly(:,1))) max(abs(epoly(:,2)))];                     
        [ratemap,dwellmap,spikemap,rmset,~] = rate_mapper([pox poy],[spx spy],rmset);  
    
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
        annotation('line',[0.03 0.10],[0.73 0.73],'LineWidth',lwid) % to position/spike plot
        
        spike_y = 0.260;        
        an0 = annotation('line',[0.03 0.03],[0.73 spike_y],'LineWidth',lwid); % down side
        an1 = annotation('arrow',[0.03 0.047],[spike_y spike_y],'LineWidth',lwid); % to spike data 
        an2 = annotation('arrow',[0.075 0.154],[spike_y spike_y],'LineWidth',lwid); % spike data to spike map
        an3 = annotation('line',[0.196 0.27],[spike_y spike_y],'LineWidth',lwid); % spike map to spike convolutions

        pos_y = 0.510;
        an8 = annotation('arrow',[0.03 0.047],[pos_y pos_y],'LineWidth',lwid); % to pos data 
        an9 = annotation('arrow',[0.097 0.154],[pos_y pos_y],'LineWidth',lwid); % pos data to dwell map
        an10 = annotation('line',[0.196 0.27],[pos_y pos_y],'LineWidth',lwid); % dwell map to convolutions
        
        fun_y = 0.045;
        an11 = annotation('line',[0.27 0.27],[pos_y fun_y],'LineWidth',lwid); % dwell map down
        an12 = annotation('arrow',[0.27 0.292],[fun_y fun_y],'LineWidth',lwid); % to adaptive equation
        an13 = annotation('line',[0.415 0.436],[fun_y fun_y],'LineWidth',lwid); % away from adaptive equation
        an14 = annotation('line',[0.436 0.436],[fun_y 0.285],'LineWidth',lwid); % convolutions up
        an15 = annotation('arrow',[0.436 0.457],[0.285 0.285],'LineWidth',lwid); % away from adaptive equation
        
        an16 = annotation('line',[0.52 0.52],[0.375 0.73],'LineWidth',lwid); % distance map up
        an17 = annotation('arrow',[0.52 0.48],[0.73 0.73],'LineWidth',lwid); % left to ratemap
        annotation('textbox','String',sprintf('For every bin:\nrate = s / n'),'Position',[0.465 0.455 0.1 0.05],'EdgeColor','none','BackgroundColor','w','HorizontalAl','center')         
        
%         delete(an2)
        an2 = annotation('textbox','String',sprintf('For every bin:\nrate = s / n \\times sample interval'),'Position',[0.420 0.455 0.20 0.05],'EdgeColor','none','BackgroundColor','w','HorizontalAl','center');       
  % keyboard
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ################################################################# %% Save the figure
    % Save the figure    
    disp(sprintf('\tSaving figure...'))    
    figname = [fig_dir2 '\' mname '_fig_meth.png'];
    [~,~,~] = mkdir(fig_dir2);    
    exportgraphics(gcf,figname,'BackgroundColor','w','ContentType','image','Resolution',250);  
    close(gcf)      
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    