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
    mname = 'histogram';
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
    ynow = 650;
    xnow = 90;
    if 1
        disp(sprintf('\tExample rate maps...'))        
        sdir = [scratch_space '\' ename '_sigma8000_duration64_sdata.mat'];            
        disp(sprintf('\t\t...loading %s',sdir));            
        load(sdir,'sdata'); 
        maps_for_plotting = cell(2,9);
        
        for aa = 1:2
            psiz = 90;
            buff = 95;
            if aa==2
                xnow = xnow+350;
            end
            xvec = repmat([xnow xnow+buff xnow+buff*2],1,3);
            yvec = [ynow ynow ynow ynow-buff ynow-buff ynow-buff ynow-buff*2 ynow-buff*2 ynow-buff*2];
%             bvnow = [5 5 5 5 5 5 5 5 5];
%             svnow = [20 20 20 20 20 20 20 20 20];       
            excell = [1 2 3 10 5 6 7 8 9];
            bvnow = ones(1,length(excell)).*20;
            svnow = ones(1,length(excell)).*30;
            
            override_maps = 1;
            fname = [fig_dir mname '_smoo_maps.mat'];
            if exist(fname,'file') && ~override_maps
                mem = 1;
                disp(sprintf('\t\t...loading old maps'))
                load(fname,'maps_for_plotting')
            else
                disp(sprintf('\t\t...generating maps'))            
                mem = 0;
            end
            loopout = looper(length(excell));        
            for xx = 1:length(excell)
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
                pspx = spx(spt<tcut);
                pspy = spy(spt<tcut);
                pos = [ppox ppoy]; % positions in mm
                spk = [pspx pspy]; % spikes in mm
                
                ax = axes('Units','pixels','Position',[xvec(xx) yvec(xx) psiz psiz]);    
                    if aa==1 && xx==1
                        ah = add_panel_title('a',sprintf('Smoothing before division'),'yoffset',15,'xoffset',0,'width',300);  
                    elseif aa==2 && xx==1 
                        ah = add_panel_title('b',sprintf('Smoothing after division'),'yoffset',15,'xoffset',0,'width',300);  
                    end

                    if mem
                        ratemap = maps_for_plotting{aa,xx};
                    else                    
                        rmset = mapset; % rate mapper settings structure - passing mapset directly to ratemapper causes a memory leak
                        rmset.method = mname;
                        rmset.binsize = bvnow(xx);
                        rmset.ssigma = svnow(xx);
                        if aa==1
                            rmset.smethod = 1;
                        else
                            rmset.smethod = 2;                                
                        end
                        rmset.maplims = [-max(abs(epoly(:,1))) -max(abs(epoly(:,2))) max(abs(epoly(:,1))) max(abs(epoly(:,2)))];                     
                        [ratemap,dwellmap,spikemap,rmset,speedlift] = rate_mapper(pos,spk,rmset);
                        maps_for_plotting{aa,xx} = ratemap;
                    end

                    im = imagesc(ratemap,'alphadata',~isnan(ratemap)); hold on;
                    daspect([1 1 1])
                    axis xy off tight
                    colormap(gca,rmap_colormap)     
                    text(0,1.05,sprintf('cell %d',xx),'HorizontalAl','left','Units','normalized','FontSize',8,'rotation',0)                                                
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

    xnow = 70;
    ynow = ynow-460;
    
    n_color_levels = 64;
    error_colormap = flipud(inferno(n_color_levels));
    clims = [1*10^-12 1*10^-11];
    interp_method = 'linear';   
    
    if 1
        disp(sprintf('\tMISE map...'))
        
        % main panel showing 16 minute error results
        ax1 = axes('Units','pixels','Position',[xnow ynow 190 190]);    
            ah = add_panel_title('c',sprintf('Smoothing before or after MISE'),'yoffset',10,'xoffset',0,'width',400);
            sigma_now = '4000';
            time_now = '16';              
            load([scratch_space '\' ename '_' mname '_sigma' num2str(sigma_now) '_duration' num2str(time_now) '_datamats.mat'],'data_matrix_error','data_matrix_params')  
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
            ax1.YLim = [1 max(new_svec)];     
            ax1.XTick = [1 2 5 10 20 40 80 160 320 640];
            ax1.YTick = [1 2 5 10 20 40 80 160 320 640];
            view(0,90)
            ylabel(sprintf('Smoothing %c (mm)',963));
            xlabel(sprintf('Bin size (mm)'));
            text(0,1.05,sprintf('Smoothing before division MISE'),'FontSize',8,'HorizontalAl','left','Units','normalized')

       % main panel showing 16 minute error results
        ax1 = axes('Units','pixels','Position',[xnow+225 ynow 190 190]);                
            load([scratch_space '\' ename '_' mname '_sigma' num2str(sigma_now) '_duration' num2str(time_now) '_sm2_datamats.mat'],'data_matrix_error','data_matrix_params')  
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
            ax1.YLim = [1 max(new_svec)];     
            ax1.XTick = [1 2 5 10 20 40 80 160 320 640];
            ax1.YTick = [1 2 5 10 20 40 80 160 320 640];
            view(0,90)
            %ylabel(sprintf('Smoothing %c (mm)',963));
            xlabel(sprintf('Bin size (mm)'));
            text(0,1.05,sprintf('Smoothing after division MISE'),'FontSize',8,'HorizontalAl','left','Units','normalized')

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
        
%% #################### Smoothing before vs after
    if 1
        xnow = xnow+520;   
        ax1 = axes('Units','pixels','Position',[xnow ynow 190 190]);   
            ah = add_panel_title('d',sprintf('Smoothing before vs after'),'yoffset',10,'xoffset',0,'width',400);

            sigma_now = config.analysis_pairs{1,1};
            time_now = config.analysis_pairs{1,2}(2); 

            load([scratch_space '\' ename '_histogram_sigma' num2str(sigma_now) '_duration' num2str(time_now) '_datamats.mat'],'data_matrix_params','data_matrix_error')    
            dme = data_matrix_error(:,:,:,mapidx);
            mat1 = mean(dme,3,'double','omitnan');  
            bmat = data_matrix_params(:,:,1); % binsizes used when generating maps
            smat = data_matrix_params(:,:,2); % smoothing used when generating maps          
            new_bvec = min(bmat(:)):1:max(bmat(:));
            new_svec = min(smat(:)):1:max(smat(:));             
            [Xq,Yq] = meshgrid(new_bvec,new_svec); 
            
            Vq1 = interp2(bmat,smat,mat1,Xq,Yq,interp_method);

            load([scratch_space '\' ename '_histogram_sigma' num2str(sigma_now) '_duration' num2str(time_now) '_sm2_datamats.mat'],'data_matrix_params','data_matrix_error')    
            dme = data_matrix_error(:,:,:,mapidx);
            mat1 = mean(dme,3,'double','omitnan');  
            bmat = data_matrix_params(:,:,1); % binsizes used when generating maps
            smat = data_matrix_params(:,:,2); % smoothing used when generating maps          
            new_bvec = min(bmat(:)):1:max(bmat(:));
            new_svec = min(smat(:)):1:max(smat(:));             
            [Xq,Yq] = meshgrid(new_bvec,new_svec); 
            
            Vq2 = interp2(bmat,smat,mat1,Xq,Yq,interp_method);        

            ediff = Vq1 - Vq2;            
            im = surf(Xq,Yq,ediff,'EdgeColor','none'); hold on;
            colormap(ax1,cmocean('balance'))
            cmap = cmocean('balance');

            axis ij
            ax1.XScale = 'log';
            ax1.YScale = 'log';
            ax1.XLim = [min(new_bvec) max(new_bvec)];
            ax1.YLim = [1 max(new_svec)];     
            ax1.XTick = [1 2 5 10 20 40 80 160 320 640];
            ax1.YTick = [1 2 5 10 20 40 80 160 320 640];
            view(0,90)
            ylabel(sprintf('Smoothing %c (mm)',963));
            xlabel(sprintf('Bin size (mm)'));
            caxis([-1*10^-10 1*10^-10])    

        axc = axes('Units','pixels','Position',[ax1.Position(1)+ax1.Position(3)+20 ax1.Position(2)+10 12 ax1.Position(4)*0.9]); 
            mat = linspace(ax1.CLim(1),ax1.CLim(2),100)';
            imagesc(ones(size(mat)),mat,mat);
            colormap(axc,ax1.Colormap);
            axis xy

            axc.YTick = linspace(ax1.CLim(1),ax1.CLim(2),3);
            axc.XTick = [];
            axc.YAxisLocation = 'right';
            text(0.5,1.15,sprintf('before - after'),'FontSize',8,'HorizontalAl','center','Units','normalized')
                
    end

% return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ################################################################# %% Save the figure
    % Save the figure    
    disp(sprintf('\tSaving figure...'))    
    figname = [fig_dir2 '\' mname '_fig_smoo.png'];
    [~,~,~] = mkdir(fig_dir2);    
    exportgraphics(gcf,figname,'BackgroundColor','w','ContentType','image','Resolution',250);  
    close(gcf) 
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    






