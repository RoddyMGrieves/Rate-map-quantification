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

    ename = 'arena120cm'; 
    mname = 'fyhn';
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

%% #################### Example cells
    ynow = 700;
    xnow = 90;
    if 1
        disp(sprintf('\tExample rate maps...'))        
        sdir = [scratch_space '\' ename '_sigma8000_duration64_sdata.mat'];            
        disp(sprintf('\t\t...loading %s',sdir));            
        load(sdir,'sdata'); 
        maps_for_plotting = cell(4,4);
        tvals = [0.02 0.125 0.25 0.5 1 1.5 2 4];
        
        for aa = 1:length(tvals)
            psiz = 120;
            buff = 95;
            xnowp = xnow+((psiz*(aa-1))+35);

            hbuff = 15;
            xvec = [xnow xnow+psiz+hbuff xnow+(psiz*2)+(hbuff*2) xnow+(psiz*3)+(hbuff*3) xnow+(psiz*4)+(hbuff*4) xnow+(psiz*5)+(hbuff*5) xnow+(psiz*6)+(hbuff*6) xnow+(psiz*7)+(hbuff*7)];
            
            excell = [1 20 38];
            
            override_maps = 1;
            fname = [fig_dir mname '_time_maps.mat'];
            if exist(fname,'file') && ~override_maps
                mem = 1;
                disp(sprintf('\t\t...loading old maps'))
                load(fname,'maps_for_plotting')
            else
                disp(sprintf('\t\t...generating maps'))            
                mem = 0;
            end
            loopout = looper(length(excell));        
            for xx = 1:3
                % cut the position and spike data to the trial length
                % we want
                pindx = sdata.pos_index(excell(xx));
                pos_now = all_pos{pindx,1};
                pox = pos_now(:,1);
                poy = pos_now(:,2);
                pot = pos_now(:,3); 

                spk = sdata.spk{excell(xx)};
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
                  
                vbuff = 20;
                yvec = [ynow ynow-psiz-vbuff ynow-(psiz*2)-(vbuff*2) ynow-(psiz*3)-(vbuff*3)];
                ax = axes('Units','pixels','Position',[xvec(aa) yvec(xx) psiz psiz]);    
                    if aa==1 && xx==1
                        ah = add_panel_title('a',sprintf('Example cells mapped using various temporal smoothing'),'yoffset',15,'xoffset',0,'width',300);    
                    end

                    if mem
                        ratemap = maps_for_plotting{aa,xx};
                    else                    
                        rmset = mapset; % rate mapper settings structure - passing mapset directly to ratemapper causes a memory leak
                        rmset.method = mname;
                        rmset.binsize = 20;
                        rmset.ssigma = 30;
                        rmset.pot = ppot;
                        rmset.spt = pspt;
                        rmset.maxdist = 320;
                        rmset.mindist = 50;
                        rmset.twindow = tvals(aa);
                        rmset.maplims = [-max(abs(epoly(:,1))) -max(abs(epoly(:,2))) max(abs(epoly(:,1))) max(abs(epoly(:,2)))];                     
                        [ratemap,dwellmap,spikemap,rmset,speedlift] = rate_mapper(pos,spk,rmset);
                        maps_for_plotting{aa,xx} = ratemap;
                    end

                    im = imagesc(ratemap,'alphadata',~isnan(ratemap)); hold on;
                    daspect([1 1 1])
                    axis xy off tight
                    colormap(gca,rmap_colormap)    
                    if aa==1
                        text(-0.1,0.5,sprintf('cell %d',xx),'HorizontalAl','center','Units','normalized','FontSize',12,'rotation',90)   
                    end
                    if xx==1 % if first cell
                        if aa==1 || aa==3
                            text(1,1.07,sprintf('%.2fs window',tvals(aa)),'HorizontalAl','right','Units','normalized','FontSize',12,'rotation',0)    
                        elseif aa==2
                            text(1,1.07,sprintf('%.3fs window',tvals(aa)),'HorizontalAl','right','Units','normalized','FontSize',12,'rotation',0)                            
                        else
                            text(1,1.07,sprintf('%.1fs window',tvals(aa)),'HorizontalAl','right','Units','normalized','FontSize',12,'rotation',0)                            
                        end
                    end
                    loopout = looper(loopout);        
            end

        end
        
        axc = axes('Units','pixels','Position',[ax.Position(1)+ax.Position(3)+20 ax.Position(2)+45 12 ax.Position(4)*2]); 
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
                   
        if ~mem
            save(fname,'maps_for_plotting');
        end        
    end

%% #################### Heatmap of error by binsize and smoothing
    mapidx = 4;
%     mapidx = 2

    xnow = 70;
    ynow = ynow-560;
    
    n_color_levels = 64;
    error_colormap = flipud(inferno(n_color_levels));    
    clims = [1*10^-12 1*10^-11];
%     clims = 'auto'
    interp_method = 'nearest';   
    
    if 1
        disp(sprintf('\tMISE map...'))
        
        hbuff = 225;
        xvec = [xnow xnow+hbuff xnow+(hbuff*2) xnow+(hbuff*3) xnow+(hbuff*4) xnow+(hbuff*5)];
        tvec = {'_0-0','_0-125','_0-25','_1-0','_2-0'};
        ttext = {'0.02','0.125','0.25','1','2'};
        for ii = 1:length(tvec)
            % main panel showing 16 minute error results
            ax1 = axes('Units','pixels','Position',[xvec(ii) ynow 190 190]);   
                if ii==1
                    ah = add_panel_title('b',sprintf('Bin size \\times smoothing MISE maps for different temporal smoothing values'),'yoffset',10,'xoffset',0,'width',400);
                end
                sigma_now = '4000';
                time_now = '16';              
                load([scratch_space '\' ename '_' mname '_sigma' num2str(sigma_now) '_duration' time_now tvec{ii} '_datamats.mat'],'data_matrix_error','data_matrix_params')  
                bmat = data_matrix_params(:,:,1); % binsizes used when generating maps
                smat = data_matrix_params(:,:,2); % smoothing used when generating maps
                new_bvec = min(bmat(:)):1:max(bmat(:));
                new_svec = min(smat(:)):1:max(smat(:));             
                [Xq,Yq] = meshgrid(new_bvec,new_svec);     

                dme = data_matrix_error(:,:,:,mapidx);
                mat1 = mean(dme,3,'double','omitnan');    

                Vq1 = interp2(bmat,smat,mat1,Xq,Yq,interp_method);
                surf(Xq,Yq,Vq1,'EdgeColor','none'); hold on;
                colormap(gca,error_colormap) 
                set(gca,'ColorScale','log')
                caxis(clims) 

                axis ij
                ax1.XScale = 'log';
                ax1.YScale = 'log';
                ax1.XLim = [min(new_bvec) max(new_bvec)];
                ax1.YLim = [min(new_svec) max(new_svec)];     
                ax1.XTick = [5 10 20 40 80 160 320 640];
                ax1.YTick = [10 20 40 80 160 320 640];
                view(0,90)
                if ii==1
                    ylabel(sprintf('Spatial smoothing %c (mm)',963));
                end
                xlabel(sprintf('Bin size (mm)'));
                text(0,1.05,sprintf('%ss window',ttext{ii}),'FontSize',12,'HorizontalAl','left','Units','normalized')
                text(1,.05,sprintf('%e',min(mat1(:))),'FontSize',12,'HorizontalAl','right','Units','normalized','Color','w')
                
                if ii==length(tvec)
                    axc = axes('Units','pixels','Position',[ax1.Position(1)+ax1.Position(3)+10 ax1.Position(2)+0 12 ax1.Position(4)*0.9]); 
                        mat = linspace(ax1.CLim(1),ax1.CLim(2),100)';
                        imagesc(ones(size(mat)),mat,mat);
                        colormap(axc,error_colormap);
                        axis xy

                        axc.YTick = linspace(ax1.CLim(1),ax1.CLim(2),4);
                        axc.XTick = [];
                        axc.YAxisLocation = 'right';
                        text(0.5,1.1,sprintf('MISE'),'FontSize',8,'HorizontalAl','center','Units','normalized')
                end  
        end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ################################################################# %% Save the figure
    % Save the figure    
    disp(sprintf('\tSaving figure...'))    
    figname = [fig_dir2 '\S4_Fig.png'];
    [~,~,~] = mkdir(fig_dir2);    
    exportgraphics(gcf,figname,'BackgroundColor','w','ContentType','image','Resolution',350);  
    close(gcf) 
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    






