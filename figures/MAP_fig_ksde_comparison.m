%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> DESCRIPTION
% DIS_rev2_fig_1  
% Script to analyse data and produce Fig. 1 for:
% Grieves, Shinder, Rosow, Kenna and Taube (2022) 
% The neural correlates of spatial disorientation in head direction cells. eNeuro
%
% USAGE:
%       MAP_generate_maps process with default settings
%
% See also: graphtest MAP_get_random_walks rate_mapper

% HISTORY:
% version 1.0.0, Release 22/11/21 Initial release
%
% Author: Roddy Grieves
% Dartmouth College, Moore Hall
% eMail: roddy.m.grieves@dartmouth.edu
% Copyright 2021 Roddy Grieves

%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Heading 3
%% >>>>>>>>>>>>>>>>>>>> Heading 2
%% >>>>>>>>>> Heading 1
%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> INPUT ARGUMENTS CHECK
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
        wlength = 64;
        sdir = [scratch_space '\' config.environs{1} '_sigma8000_duration64_sdata.mat'];            
        disp(sprintf('\t\t...loading %s',sdir));            
        load(sdir,'sdata'); 
    end

    ename = 'arena120cm';
    mname = 'ksde';
    rmap_colormap = 'turbo';
    
%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> FUNCTION BODY
    % create a figure
    figure('visible',fig_vis,'Position',[100 50 1400 800]);
    set(gcf,'InvertHardCopy','off'); % this stops white lines being plotted black      
    set(gcf,'color','w'); % makes the background colour white

%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Protocol Exp 1
    xnow = 50;
    ynow = 600;
    if 1
        xvec = xnow+300:110:10000;
        yvec = ynow:-110:-1000;
        m = {'leutgeb','leutgeb_pixelwise','ksde'};
        m2 = {'KSDE\n(Leutgeb kernel)','Leutgeb method\n(pixelwise)','KSDE\n(Matlab kernel)'};
        
        bnow = 20;
        snow = 40;  
        excell = [2 8 6 5];        
        
        psiz = 100;
        nmeths = 3;
        nplots = 4;
        maps_for_plotting = cell(nplots,nmeths);
        override_maps = 0;
        fname = [fig_dir 'ksde_comparison_maps.mat'];
        if exist(fname,'file') && ~override_maps
            mem = 1;
            disp(sprintf('\t...loading old maps'))
            load(fname,'maps_for_plotting')
        else
            mem = 0;
        end
        loopout = looper(nmeths*nplots);           
        for xx = 1:nmeths
            ts = NaN(1,nplots);
            for yy = 1:nplots
                % cut the position and spike data to the trial length
                pindx = sdata.pos_index(excell(yy));
                pos_now = all_pos{pindx,1};
                pox = pos_now(:,1);
                poy = pos_now(:,2);
                pot = pos_now(:,3); 

                spk = sdata.spk{excell(yy)};
                spx = spk(:,1);
                spy = spk(:,2);
                spt = spk(:,3);

                tcut = 16*60;
                ppox = pox(pot<tcut);
                ppoy = poy(pot<tcut);
                pspx = spx(spt<tcut);
                pspy = spy(spt<tcut);                
                pos = [ppox ppoy]; % positions in mm
                spk = [pspx pspy]; % spikes in mm
        
                if xx==1
                    ax = axes('Units','pixels','Position',[xvec(xx)-250 yvec(yy) psiz psiz]);    
                        im = imagesc(sdata.pmap_1mm{excell(yy)}); hold on;
                        daspect([1 1 1])
                        axis xy off tight
                        colormap(gca,rmap_colormap)   
                        if yy==1
                            text(0.5,1.3,sprintf('Actual spike\nprobability'),'HorizontalAl','center','VerticalAl','top','FontSize',8,'Units','normalized')
                        end
                        
                    ax = axes('Units','pixels','Position',[xvec(xx)-125 yvec(yy) psiz psiz]);    
                        pindx = sdata.pos_index(excell(yy));
                        pos_now = all_pos{pindx,1};
                        pox = pos_now(:,1);
                        poy = pos_now(:,2);
                        sp = sdata.spk{excell(yy)};
                        spx = sp(:,1);
                        spy = sp(:,2);

                        plot(pox,poy,'Color',[.5 .5 .5 .5]); hold on;
                        scatter(spx,spy,5,'r','filled','MarkerEdgeColor','none','MarkerFaceAlpha',0.5);
                        daspect([1 1 1])
                        axis xy off                          
                    
                        daspect([1 1 1])
                        axis xy off tight
                        colormap(gca,flipud(viridis))  
                        if yy==1
                            text(0.5,1.3,sprintf('Spikes &\nPositions'),'HorizontalAl','center','VerticalAl','top','FontSize',8,'Units','normalized')  
                        end
                end

                ax = axes('Units','pixels','Position',[xvec(xx) yvec(yy) psiz psiz]);    
                    if mem
                        ratemap = maps_for_plotting{xx,yy,1};
                        ts(yy) = maps_for_plotting{xx,yy,2};                       
                    else                     
                        rmset = mapset; % rate mapper settings structure - passing mapset directly to ratemapper causes a memory leak
                        rmset.method = m{xx};
                        rmset.binsize = bnow;
                        rmset.ssigma = snow;
                        rmset.maplims = [-max(abs(epoly(:,1))) -max(abs(epoly(:,2))) max(abs(epoly(:,1))) max(abs(epoly(:,2)))];   
                        tic;
                        [ratemap,dwellmap,spikemap,rmset,speedlift] = rate_mapper(pos,spk,rmset);
                        maps_for_plotting{xx,yy,1} = ratemap;
                        maps_for_plotting{xx,yy,2} = toc; 
                        ts(yy) = maps_for_plotting{xx,yy,2};                                               
                    end
                    
                    im = imagesc(ratemap); hold on;
                    daspect([1 1 1])
                    axis xy off tight
                    colormap(gca,rmap_colormap)  
                    if yy==1
                        text(0.5,1.3,sprintf(m2{xx}),'HorizontalAl','center','VerticalAl','top','FontSize',8,'Units','normalized')
                    end                    
                    if yy==nplots
                        %text(0,-0.1,sprintf('Average computation\ntime: %.1fms',mean(ts(:),'omitnan')*1000),'HorizontalAl','left','VerticalAl','top','FontSize',8,'Units','normalized')                        
                    end
                    
                    loopout = looper(loopout);                            
            end
        end
        if ~mem
            save(fname,'maps_for_plotting');
        end

        axc = axes('Units','pixels','Position',[ax.Position(1)+120 ax.Position(2)+200 12 ax.Position(4)*2]); 
            mat = linspace(0,100,100)';
            imagesc(ones(size(mat)),mat,mat);
            colormap(axc,ax.Colormap);
            axis xy

            axc.YTick = [0 100];
            axc.YTickLabel = {'0 Hz','Max'};
            axc.XTick = [];
            axc.YAxisLocation = 'right';
            text(0,1.15,sprintf('Firing\nrate (Hz)'),'FontSize',8,'HorizontalAl','left','Units','normalized')
            % add univisted box
            axc2 = axes('Units','pixels','Position',[axc.Position(1) axc.Position(2)-20 axc.Position(3) 10]); 
            axis on; box on; axc2.XTick = []; axc2.YTick = 0.5; axc2.YTickLabel = {'Unvisited'}; axc2.YAxisLocation = 'right';
            
        ax = axes('Units','pixels','Position',[ax.Position(1)-220 ax.Position(2)-110 320 80]);
            v1 = mean([maps_for_plotting{1,:,2}],'omitnan');
            v2 = mean([maps_for_plotting{2,:,2}],'omitnan');
            v3 = mean([maps_for_plotting{3,:,2}],'omitnan');
            bar(1:3,[v1 v2 v3],0.5,'k');
            
            ax.XTick = [];
            ax.XLim = [0.5 3.5];
            ax.YLim = [0 5];
            ax.XTick = 1:3;
            ax.XTickLabel = [];
            for ii = 1:3
                text(ii,-0.5,sprintf(m2{ii}),'VerticalAlignment','top','HorizontalAlignment','center','FontSize',8)
            end
            ylabel(sprintf('Computation time (s)'))
            box off  
    end

% return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ################################################################# %% Save the figure
    % Save the figure    
    disp(sprintf('\tSaving figure...'))    
    figname = [fig_dir2 '\S3_Fig.png'];
    [~,~,~] = mkdir(fig_dir2);    
    exportgraphics(gcf,figname,'BackgroundColor','w','ContentType','image','Resolution',350);  
    close(gcf)





    
    