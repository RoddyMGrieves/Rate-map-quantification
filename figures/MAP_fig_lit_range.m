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
    mname = 'histogram';
    rmap_colormap = 'turbo';
    
    fname = 'C:\Users\F004KS7\OneDrive - Dartmouth College\Projects in prep\2019 Mapping project\associated data\papers.csv';
    dat = readtable(fname);
    dat = table2array(dat);

    fname = 'C:\Users\F004KS7\OneDrive - Dartmouth College\Projects in prep\2019 Mapping project\associated data\boxcar.csv';
    dat_box = readtable(fname);
    dat_box = table2array(dat_box);

    fname = 'C:\Users\F004KS7\OneDrive - Dartmouth College\Projects in prep\2019 Mapping project\associated data\adaptive.csv';
    dat_adapt = readtable(fname);
    dat_adapt = table2array(dat_adapt);    
    
    fname = 'C:\Users\F004KS7\OneDrive - Dartmouth College\Projects in prep\2019 Mapping project\associated data\ksde.csv';
    dat_ksde = readtable(fname);
    dat_ksde = table2array(dat_ksde);      
    
    switch mname
        case {'histogram'}
            bvnow = [1 10 20 40 80 160 320];
            svnow = [1 10 20 40 80 160 320];
            smoo_text = sprintf('Smoothing %c (mm)',963);
            grid_text = sprintf('Bin size (mm)');
        case {'ash'}
            bvnow = [20 40 60 80 160 320 640];
            svnow = [2 4 8 16 32 48 64]; 
            smoo_text = sprintf('Smoothing %c',948);
            grid_text = sprintf('Bin size (mm)');                
        case {'ksde'}
            bvnow = [2.5 5 10 20 40 160 320];
            svnow = [10 20 30 40 80 160 320];  
            smoo_text = sprintf('Smoothing bandwidth (mm)');
            grid_text = sprintf('Bin size (mm)');                
        case {'kyadaptive'}
            bvnow = [2.5 5 10 20 80 160 320];
            svnow = [0.5 1 1.5 2 4 6 10]; 
            smoo_text = 'Smoothing {\itt}';
            grid_text = sprintf('Bin size (mm)');         
        case {'kadaptive'}
            bvnow = [1 3 5 10 40 160 320];
            svnow = [100 1000 6000 10000 16000 24000 32000];  
            smoo_text = sprintf('Smoothing %c',945);
            grid_text = sprintf('Bin size (mm)');  
        case {'fyhn'}
            bvnow = [5 10 20 40 80 160 320];
            svnow = [10 20 30 40 80 160 320];
            smoo_text = sprintf('Smoothing %c (mm)',963);
            grid_text = sprintf('Bin size (mm)');                   
        otherwise
            keyboard
    end

%% ##################################### Heading 2
%% #################### Heading 3
%% Heading 4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ################################################################# %% FUNCTION BODY
        % Figure settings
        fig1 = figure('visible','on','Position',[50,50,1300,1050]); 
        set(gcf,'InvertHardCopy','off'); % this stops white lines being plotted black      
        set(gcf,'color','w'); % makes the background colour white
        colormap(jet(256)); % to make sure the colormap is not the horrible default one


%% #################### Literature values on MISE heatmap


    ynow = 150;
    ynow = 500;
    xsiz = 200;
    ysiz = 200;

    disp(sprintf('\tMISE map...'))
    n_color_levels = 64;
    error_colormap = flipud(inferno(n_color_levels));    
    clims = [1*10^-12 1*10^-11];
    interp_method = 'nearest';  
    mapidx = 4;
    transp = 0.5;
    xbuff = xsiz+70;
    xvec = [xnow xnow+xbuff xnow+xbuff*2]+30;

    ms = {'histogram','kadaptive','ksde'};
    ds = {'Histogram','Adaptive smoothing','KSDE'};
    for ii = 1:3
        % main panel showing 16 minute error results
        mname = ms{ii};
        ax1 = axes('Units','pixels','Position',[xvec(ii) ynow-5 xsiz ysiz]);    
            if ii==1
                ah = add_panel_title('b',sprintf('Reported bin size and smoothing combinations'),'yoffset',35,'xoffset',-30,'width',400);
            end
            sigma_now = '4000';  
            time_now = '16';                        
            load([scratch_space '\' ename '_' mname '_sigma' num2str(sigma_now) '_duration' num2str(time_now) '_datamats.mat'],'data_matrix_error','data_matrix_fields','data_matrix_cputime','data_matrix_params')    
            dme = data_matrix_error(:,:,:,mapidx);
            bmat = double(data_matrix_params(:,:,1)); % binsizes used when generating maps
            smat = double(data_matrix_params(:,:,2)); % smoothing used when generating maps             
            switch mname
                case {'histogram','fyhn'}
                    smoo_text = sprintf('Smoothing %c (mm)',963);
                    grid_text = sprintf('Bin size (mm)');                    
                    new_bvec = min(bmat(:)):1:max(bmat(:));
                    new_svec = min(smat(:)):1:max(smat(:));  
                case {'ash'}
                    smoo_text = 'Smoothing {\itm}';
                    grid_text = sprintf('Bin size (mm)');                      
                    new_bvec = min(bmat(:)):1:max(bmat(:));
                    new_svec = min(smat(:)):0.1:max(smat(:));  
                case {'kadaptive'}
                    smoo_text = sprintf('Smoothing %c',945);
                    grid_text = sprintf('Bin size (mm)');                     
                    new_bvec = min(bmat(:)):1:max(bmat(:));
                    new_svec = min(smat(:)):100:max(smat(:)); 
                case {'kyadaptive'}
                    smoo_text = 'Smoothing {\itt}';
                    grid_text = sprintf('Bin size (mm)');                    
                    new_bvec = min(bmat(:)):1:max(bmat(:));
                    new_svec = min(smat(:)):0.1:max(smat(:));                        
                case {'ksde'}
                    smoo_text = sprintf('Smoothing bandwidth (mm)');
                    grid_text = sprintf('Bin size (mm)');                      
                    new_bvec = min(bmat(:)):1:max(bmat(:));
                    new_svec = min(smat(:)):0.5:max(smat(:));                    
                otherwise
                    keyboard 
            end             
            [Xq,Yq] = meshgrid(new_bvec,new_svec); 
    
            mat1 = mean(dme,3,'double','omitnan');    
            Vq1 = interp2(bmat,smat,mat1,Xq,Yq,interp_method);
    
            surf(Xq,Yq,Vq1,'EdgeColor','none','FaceAlpha',transp); hold on;
            colormap(gca,error_colormap) 
            set(gca,'ColorScale','log')
            grid off
            caxis(clims) 
    
            % find the 'best' solution (closest to minimising everything)
            load([scratch_space '\' ename '_' mname '_sigma' num2str(sigma_now) '_duration' num2str(time_now) '_gamultiobj.mat'],'balance_solution2','balance_solution1','minerr_solution','x')    
    
            % plot literature values
            zidx = dat(:,3)>0;
    
            jitter = 2;
            switch mname
                case {'histogram'}
                    x = dat(zidx,2) + normrnd(0,jitter,size(dat(zidx,2)));   
                    y = dat(zidx,3) + normrnd(0,jitter,size(dat(zidx,2)));  
                    pc2 = plot3(x,y,repmat(1000,size(dat(zidx,2))),'k*','MarkerSize',10,'LineStyle','none');     
                    text(0,1.14,sprintf(ds{ii}),'Units','normalized','HorizontalAlignment','left','VerticalAlignment','bottom','FontSize',9)                    
                case {'kadaptive'}
                    x = dat_adapt(:,2) + normrnd(0,jitter,size(dat_adapt(:,2)));   
                    y = dat_adapt(:,3) + normrnd(0,jitter,size(dat_adapt(:,2)));  
                    pc2 = plot3(x,y,repmat(1000,size(dat_adapt(:,2))),'k*','MarkerSize',10,'LineStyle','none');   
                    text(0,1.02,sprintf(ds{ii}),'Units','normalized','HorizontalAlignment','left','VerticalAlignment','bottom','FontSize',9)                    
                case {'ksde'}
                    x = dat_ksde(:,2) + normrnd(0,jitter,size(dat_ksde(:,2)));   
                    y = dat_ksde(:,3) + normrnd(0,jitter,size(dat_ksde(:,2)));  
                    pc2 = plot3(x,y,repmat(1000,size(dat_ksde(:,2))),'k*','MarkerSize',10,'LineStyle','none');    
                    text(0,1.02,sprintf(ds{ii}),'Units','normalized','HorizontalAlignment','left','VerticalAlignment','bottom','FontSize',9)                    
            end
            pb = plot3(balance_solution1(1),balance_solution1(2),1000,'ko','MarkerSize',8,'LineStyle','none','MarkerFaceColor','k');            
            pc = plot3(minerr_solution(1),minerr_solution(2),1000,'kv','MarkerSize',8,'LineStyle','none','MarkerFaceColor','k');  

            if ii==2
                [~,leg] = legendflex([pb pc pc2]...
                    ,{'Balanced solution','Minimum error solution','Literature reported values'}...
                    ,'anchor',{'n','n'},'ncol',1,'box','off','buffer',[30,80],'xscale',.5,'fontsize',9);  
            end
    
            axis ij
            ax1.XScale = 'log';
            ax1.YScale = 'log';  
            ax1.XLim = [min(new_bvec) max(new_bvec)];
            ax1.YLim = [min(new_svec) max(new_svec)];             
            switch mname
                case {'histogram'}
                    ax1.YLim = [1 max(new_svec)];    
                    ax1.XTick = [1 2 5 10 20 40 80 160 320 640];
                    ax1.YTick = [1 2 5 10 20 40 80 160 320 640];
                case {'ash'}                     
                    ax1.XTick = [1 2 5 10 20 40 80 160 320 640];
                    ax1.YTick = [0.5 1 2 4 8 16 32 64];     
                case {'kadaptive'}
                    ax1.XTick = [1 2 5 10 20 40 80 160 320 640];
                    ax1.YTick = [128 256 512 1024 2048 4096 8192 16000 32000];   
                case {'kyadaptive'}
                    ax1.XTick = [1 2 5 10 20 40 80 160 320 640];
                    ax1.YTick = [0.5 1 2 4 8];                        
                case {'ksde'}     
                    ax1.XTick = [2.5 5 10 20 40 80 160 320 640];
                    ax1.YTick = [10 20 40 80 160 320];    
                case {'fyhn'}
                    ax1.XTick = [5 10 20 40 80 160 320 640];
                    ax1.YTick = [10 20 40 80 160 320 640];                        
                otherwise
                    keyboard
            end               
            view(0,90)         
            ylabel(smoo_text);
            xlabel(grid_text);
    
            if strcmp(mname,'histogram')
                ax1b = axes('Units','pixels','Position',[ax1.Position(1) ax1.Position(2)+ax1.Position(4)+6 ax1.Position(3) 20]); 
                    surf(Xq,Yq,Vq1,'EdgeColor','none','FaceAlpha',transp); hold on;
                    pc = plot(dat(~zidx,2),dat(~zidx,3)+0.5,'k*','MarkerSize',10,'LineStyle','none');          
    
                    colormap(gca,error_colormap) 
                    set(gca,'ColorScale','log')
    
                    axis ij
                    ax1b.XLim = [min(new_bvec) max(new_bvec)];
                    ax1b.YLim = [0 1];   
                    ax1b.XScale = 'log';                  
                    caxis(ax1.CLim);
                    ax1b.YTick = 0.5;
                    ax1b.YTickLabel = {'none'};
                    ax1b.XTick = [];      
                    view(0,90)
            end

    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ################################################################# %% Save the figure
    % Save the figure    
    disp(sprintf('\tSaving figure...'))    
    figname = [fig_dir2 '\v2_fig_lit_review.png'];
    [~,~,~] = mkdir(fig_dir2);    
    exportgraphics(fig1,figname,'BackgroundColor','w','ContentType','image','Resolution',250);  
    close(fig1)    
    
    

% 
%     return
% 
% 
% 
% 
% 
% %% #################### Literature values on MISE heatmap
%         xnow = 70;
%         ynow = 300;
% 
%         disp(sprintf('\tMISE map...'))
%         n_color_levels = 64;
%         error_colormap = flipud(inferno(n_color_levels));    
%         clims = [1*10^-12 1*10^-11];
%         interp_method = 'nearest';  
%         mapidx = 4;
%         transp = 0.5;
% 
%         % main panel showing 16 minute error results
%         ax1 = axes('Units','pixels','Position',[xnow ynow-5 300 300]);    
%             ah = add_panel_title('a',sprintf('Reported bin size and smoothing combinations'),'yoffset',25,'xoffset',-20,'width',400);
%             sigma_now = '4000';  
%             time_now = '16';                        
%             load([scratch_space '\' ename '_' mname '_sigma' num2str(sigma_now) '_duration' num2str(time_now) '_datamats.mat'],'data_matrix_error','data_matrix_fields','data_matrix_cputime','data_matrix_params')    
%             dme = data_matrix_error(:,:,:,mapidx);
%             bmat = double(data_matrix_params(:,:,1)); % binsizes used when generating maps
%             smat = double(data_matrix_params(:,:,2)); % smoothing used when generating maps         
%             switch mname
%                 case {'histogram','fyhn'}
%                     new_bvec = min(bmat(:)):1:max(bmat(:));
%                     new_svec = min(smat(:)):1:max(smat(:));  
%                 case {'ash'}
%                     new_bvec = min(bmat(:)):1:max(bmat(:));
%                     new_svec = min(smat(:)):0.1:max(smat(:));  
%                 case {'kadaptive'}
%                     new_bvec = min(bmat(:)):1:max(bmat(:));
%                     new_svec = min(smat(:)):100:max(smat(:)); 
%                 case {'kyadaptive'}
%                     new_bvec = min(bmat(:)):1:max(bmat(:));
%                     new_svec = min(smat(:)):0.1:max(smat(:));                        
%                 case {'ksde'}
%                     new_bvec = min(bmat(:)):1:max(bmat(:));
%                     new_svec = min(smat(:)):0.5:max(smat(:));                    
%                 otherwise
%                     keyboard 
%             end                          
%             [Xq,Yq] = meshgrid(new_bvec,new_svec); 
% 
%             mat1 = mean(dme,3,'double','omitnan');    
%             Vq1 = interp2(bmat,smat,mat1,Xq,Yq,interp_method);
% 
%             surf(Xq,Yq,Vq1,'EdgeColor','none','FaceAlpha',transp); hold on;
%             colormap(gca,error_colormap) 
%             set(gca,'ColorScale','log')
%             grid off
%             caxis(clims) 
% 
%             % find the 'best' solution (closest to minimising everything)
%             load([scratch_space '\' ename '_' mname '_sigma' num2str(sigma_now) '_duration' num2str(time_now) '_gamultiobj.mat'],'balance_solution2','balance_solution1','minerr_solution','x')    
% 
%             % plot literature values
%             zidx = dat(:,3)>0;
% 
%             jitter = 2;
%             x = dat(zidx,2) + normrnd(0,jitter,size(dat(zidx,2)));   
%             y = dat(zidx,3) + normrnd(0,jitter,size(dat(zidx,2)));  
%             pc2 = plot3(x,y,repmat(1000,size(dat(zidx,2))),'k*','MarkerSize',10,'LineStyle','none');          
% 
%             pb1 = plot3(balance_solution1(1),balance_solution1(2),1000,'ko','MarkerSize',8,'LineStyle','none','MarkerFaceColor','k');            
%             pc1 = plot3(minerr_solution(1),minerr_solution(2),1000,'kv','MarkerSize',8,'LineStyle','none','MarkerFaceColor','k');             
% 
%             [~,leg] = legendflex([pb1 pc1 pc2]...
%                 ,{'Pareto-optimal','Minimum error','Reported values'}...
%                 ,'anchor',{'se','se'},'ncol',1,'box','off','buffer',[-20,-120],'xscale',.5,'fontsize',9);  
% 
%             axis ij
%             ax1.XScale = 'log';
%             ax1.YScale = 'log';  
%             ax1.XLim = [min(new_bvec) max(new_bvec)];
%             ax1.YLim = [min(new_svec) max(new_svec)];             
%             switch mname
%                 case {'histogram'}
%                     ax1.YLim = [1 max(new_svec)];    
%                     ax1.XTick = [1 2 5 10 20 40 80 160 320 640];
%                     ax1.YTick = [1 2 5 10 20 40 80 160 320 640];
%                 case {'ash'}                     
%                     ax1.XTick = [1 2 5 10 20 40 80 160 320 640];
%                     ax1.YTick = [0.5 1 2 4 8 16 32 64];     
%                 case {'kadaptive'}
%                     ax1.XTick = [1 2 5 10 20 40 80 160 320 640];
%                     ax1.YTick = [128 256 512 1024 2048 4096 8192 16000 32000];   
%                 case {'kyadaptive'}
%                     ax1.XTick = [1 2 5 10 20 40 80 160 320 640];
%                     ax1.YTick = [0.5 1 2 4 8];                        
%                 case {'ksde'}     
%                     ax1.XTick = [2.5 5 10 20 40 80 160 320 640];
%                     ax1.YTick = [10 20 40 80 160 320];    
%                 case {'fyhn'}
%                     ax1.XTick = [5 10 20 40 80 160 320 640];
%                     ax1.YTick = [10 20 40 80 160 320 640];                        
%                 otherwise
%                     keyboard
%             end               
%             view(0,90)         
%             ylabel(smoo_text);
%             xlabel(grid_text);
% 
%             if strcmp(mname,'histogram')
%                 ax1b = axes('Units','pixels','Position',[ax1.Position(1) ax1.Position(2)+ax1.Position(4)+6 ax1.Position(3) 20]); 
%                     surf(Xq,Yq,Vq1,'EdgeColor','none','FaceAlpha',transp); hold on;
%                     pc = plot(dat(~zidx,2),dat(~zidx,3)+0.5,'k*','MarkerSize',10,'LineStyle','none');          
% 
%                     colormap(gca,error_colormap) 
%                     set(gca,'ColorScale','log')
% 
%                     axis ij
%                     ax1b.XLim = [min(new_bvec) max(new_bvec)];
%                     ax1b.YLim = [0 1];   
%                     ax1b.XScale = 'log';                  
%                     caxis(ax1.CLim);
%                     ax1b.YTick = 0.5;
%                     ax1b.YTickLabel = {'none'};
%                     ax1b.XTick = [];      
%                     view(0,90)
%             end
% 
% 
% %% #################### Second figure part
%     % Figure settings
%     fig1 = figure('visible','on','Position',[50,50,1200,700]); 
%     set(gcf,'InvertHardCopy','off'); % this stops white lines being plotted black      
%     set(gcf,'color','w'); % makes the background colour white
%     colormap(jet(256)); % to make sure the colormap is not the horrible default one
% 
% %% #################### Literature values on MISE heatmap (adaptive)
%     xnow = 70;
%     ynow = 300;
% 
%     disp(sprintf('\tMISE map...'))
%     n_color_levels = 64;
%     error_colormap = flipud(inferno(n_color_levels));    
%     clims = [1*10^-12 1*10^-11];
%     interp_method = 'nearest';  
%     mapidx = 4;
%     transp = 0.5;
% 
%     % main panel showing 16 minute error results
%     ax1 = axes('Units','pixels','Position',[xnow ynow-5 300 300]);    
%         ah = add_panel_title('c',sprintf('Adaptive smoothing reported bin size and smoothing'),'yoffset',5,'xoffset',-20,'width',400);
%         sigma_now = '4000';  
%         time_now = '16';             
%         mname = 'kadaptive';
%         load([scratch_space '\' ename '_' mname '_sigma' num2str(sigma_now) '_duration' num2str(time_now) '_datamats.mat'],'data_matrix_error','data_matrix_fields','data_matrix_cputime','data_matrix_params')    
%         dme = data_matrix_error(:,:,:,mapidx);
%         bmat = double(data_matrix_params(:,:,1)); % binsizes used when generating maps
%         smat = double(data_matrix_params(:,:,2)); % smoothing used when generating maps         
%         switch mname
%             case {'histogram','fyhn'}
%                 new_bvec = min(bmat(:)):1:max(bmat(:));
%                 new_svec = min(smat(:)):1:max(smat(:));  
%             case {'ash'}
%                 new_bvec = min(bmat(:)):1:max(bmat(:));
%                 new_svec = min(smat(:)):0.1:max(smat(:));  
%             case {'kadaptive'}
%                 new_bvec = min(bmat(:)):1:max(bmat(:));
%                 new_svec = min(smat(:)):100:max(smat(:)); 
%                 smoo_text = sprintf('Smoothing %c',945);
%                 grid_text = sprintf('Bin size (mm)');                  
%             case {'kyadaptive'}
%                 new_bvec = min(bmat(:)):1:max(bmat(:));
%                 new_svec = min(smat(:)):0.1:max(smat(:));                        
%             case {'ksde'}
%                 new_bvec = min(bmat(:)):1:max(bmat(:));
%                 new_svec = min(smat(:)):0.5:max(smat(:));  
%                 smoo_text = sprintf('Smoothing bandwidth (mm)');
%                 grid_text = sprintf('Bin size (mm)');                  
%             otherwise
%                 keyboard 
%         end                          
%         [Xq,Yq] = meshgrid(new_bvec,new_svec); 
% 
%         mat1 = mean(dme,3,'double','omitnan');    
%         Vq1 = interp2(bmat,smat,mat1,Xq,Yq,interp_method);
% 
%         surf(Xq,Yq,Vq1,'EdgeColor','none','FaceAlpha',transp); hold on;
%         colormap(gca,error_colormap) 
%         set(gca,'ColorScale','log')
%         grid off
%         caxis(clims) 
% 
%         % find the 'best' solution (closest to minimising everything)
%         load([scratch_space '\' ename '_' mname '_sigma' num2str(sigma_now) '_duration' num2str(time_now) '_gamultiobj.mat'],'balance_solution2','balance_solution1','minerr_solution','x')    
% 
%         % plot literature values        
%         jitter = 2;
%         x = dat_adapt(:,2) + normrnd(0,jitter,size(dat_adapt(:,2)));   
%         y = dat_adapt(:,3) + normrnd(0,jitter,size(dat_adapt(:,2)));  
%         pc2 = plot3(x,y,repmat(1000,size(dat_adapt(:,2))),'k*','MarkerSize',10,'LineStyle','none');          
% 
%         pb1 = plot3(balance_solution1(1),balance_solution1(2),1000,'ko','MarkerSize',8,'LineStyle','none','MarkerFaceColor','k');            
%         pc1 = plot3(minerr_solution(1),minerr_solution(2),1000,'kv','MarkerSize',8,'LineStyle','none','MarkerFaceColor','k');             
% 
% %         [~,leg] = legendflex([pb1 pc1 pc2]...
% %             ,{'Pareto-optimal','Minimum error','Reported values'}...
% %             ,'anchor',{'se','se'},'ncol',1,'box','off','buffer',[-20,-100],'xscale',.5,'fontsize',9);  
% 
%         axis ij
%         ax1.XScale = 'log';
%         ax1.YScale = 'log';  
%         ax1.XLim = [min(new_bvec) max(new_bvec)];
%         ax1.YLim = [min(new_svec) max(new_svec)];             
%         switch mname
%             case {'histogram'}
%                 ax1.YLim = [1 max(new_svec)];    
%                 ax1.XTick = [1 2 5 10 20 40 80 160 320 640];
%                 ax1.YTick = [1 2 5 10 20 40 80 160 320 640];
%             case {'ash'}                     
%                 ax1.XTick = [1 2 5 10 20 40 80 160 320 640];
%                 ax1.YTick = [0.5 1 2 4 8 16 32 64];     
%             case {'kadaptive'}
%                 ax1.XTick = [1 2 5 10 20 40 80 160 320 640];
%                 ax1.YTick = [128 256 512 1024 2048 4096 8192 16000 32000];   
%                 smoo_text = sprintf('Smoothing %c',945);
%                 grid_text = sprintf('Bin size (mm)');                  
%             case {'kyadaptive'}
%                 ax1.XTick = [1 2 5 10 20 40 80 160 320 640];
%                 ax1.YTick = [0.5 1 2 4 8];                        
%             case {'ksde'}     
%                 ax1.XTick = [2.5 5 10 20 40 80 160 320 640];
%                 ax1.YTick = [10 20 40 80 160 320]; 
%                 smoo_text = sprintf('Smoothing bandwidth (mm)');
%                 grid_text = sprintf('Bin size (mm)');                  
%             case {'fyhn'}
%                 ax1.XTick = [5 10 20 40 80 160 320 640];
%                 ax1.YTick = [10 20 40 80 160 320 640];                        
%             otherwise
%                 keyboard
%         end               
%         view(0,90)         
%         ylabel(smoo_text);
%         xlabel(grid_text);
% 
% %% #################### Literature values on MISE heatmap (KSDE)
%     xnow = xnow+450;
% 
%     disp(sprintf('\tMISE map...'))
%     n_color_levels = 64;
%     error_colormap = flipud(inferno(n_color_levels));    
%     clims = [1*10^-12 1*10^-11];
%     interp_method = 'nearest';  
%     mapidx = 4;
%     transp = 0.5;
% 
%     % main panel showing 16 minute error results
%     ax1 = axes('Units','pixels','Position',[xnow ynow-5 300 300]);    
%         ah = add_panel_title('d',sprintf('KSDE reported bin size and smoothing bandwidth'),'yoffset',5,'xoffset',-20,'width',400);
%         sigma_now = '4000';  
%         time_now = '16';             
%         mname = 'ksde';
%         load([scratch_space '\' ename '_' mname '_sigma' num2str(sigma_now) '_duration' num2str(time_now) '_datamats.mat'],'data_matrix_error','data_matrix_fields','data_matrix_cputime','data_matrix_params')    
%         dme = data_matrix_error(:,:,:,mapidx);
%         bmat = double(data_matrix_params(:,:,1)); % binsizes used when generating maps
%         smat = double(data_matrix_params(:,:,2)); % smoothing used when generating maps         
%         switch mname
%             case {'histogram','fyhn'}
%                 new_bvec = min(bmat(:)):1:max(bmat(:));
%                 new_svec = min(smat(:)):1:max(smat(:));  
%             case {'ash'}
%                 new_bvec = min(bmat(:)):1:max(bmat(:));
%                 new_svec = min(smat(:)):0.1:max(smat(:));  
%             case {'kadaptive'}
%                 new_bvec = min(bmat(:)):1:max(bmat(:));
%                 new_svec = min(smat(:)):100:max(smat(:)); 
%                 smoo_text = sprintf('Smoothing %c',945);
%                 grid_text = sprintf('Bin size (mm)');                 
%             case {'kyadaptive'}
%                 new_bvec = min(bmat(:)):1:max(bmat(:));
%                 new_svec = min(smat(:)):0.1:max(smat(:));                        
%             case {'ksde'}
%                 new_bvec = min(bmat(:)):1:max(bmat(:));
%                 new_svec = min(smat(:)):0.5:max(smat(:));  
%                 smoo_text = sprintf('Smoothing bandwidth (mm)');
%                 grid_text = sprintf('Bin size (mm)');                  
%             otherwise
%                 keyboard 
%         end                          
%         [Xq,Yq] = meshgrid(new_bvec,new_svec); 
% 
%         mat1 = mean(dme,3,'double','omitnan');    
%         Vq1 = interp2(bmat,smat,mat1,Xq,Yq,interp_method);
% 
%         surf(Xq,Yq,Vq1,'EdgeColor','none','FaceAlpha',transp); hold on;
%         colormap(gca,error_colormap) 
%         set(gca,'ColorScale','log')
%         grid off
%         caxis(clims) 
% 
%         % find the 'best' solution (closest to minimising everything)
%         load([scratch_space '\' ename '_' mname '_sigma' num2str(sigma_now) '_duration' num2str(time_now) '_gamultiobj.mat'],'balance_solution2','balance_solution1','minerr_solution','x')    
% 
%         % plot literature values  
%         jitter = 2;
%         x = dat_ksde(:,2) + normrnd(0,jitter,size(dat_ksde(:,2)));   
%         y = dat_ksde(:,3) + normrnd(0,jitter,size(dat_ksde(:,2)));  
%         pc2 = plot3(x,y,repmat(1000,size(dat_ksde(:,2))),'k*','MarkerSize',10,'LineStyle','none');          
% 
%         pb1 = plot3(balance_solution1(1),balance_solution1(2),1000,'ko','MarkerSize',8,'LineStyle','none','MarkerFaceColor','k');            
%         pc1 = plot3(minerr_solution(1),minerr_solution(2),1000,'kv','MarkerSize',8,'LineStyle','none','MarkerFaceColor','k');             
% 
% %         [~,leg] = legendflex([pb1 pc1 pc2]...
% %             ,{'Pareto-optimal','Minimum error','Reported values'}...
% %             ,'anchor',{'se','se'},'ncol',1,'box','off','buffer',[-20,-100],'xscale',.5,'fontsize',9);  
% 
%         axis ij
%         ax1.XScale = 'log';
%         ax1.YScale = 'log';  
%         ax1.XLim = [min(new_bvec) max(new_bvec)];
%         ax1.YLim = [min(new_svec) max(new_svec)];             
%         switch mname
%             case {'histogram'}
%                 ax1.YLim = [1 max(new_svec)];    
%                 ax1.XTick = [1 2 5 10 20 40 80 160 320 640];
%                 ax1.YTick = [1 2 5 10 20 40 80 160 320 640];
%             case {'ash'}                     
%                 ax1.XTick = [1 2 5 10 20 40 80 160 320 640];
%                 ax1.YTick = [0.5 1 2 4 8 16 32 64];     
%             case {'kadaptive'}
%                 ax1.XTick = [1 2 5 10 20 40 80 160 320 640];
%                 ax1.YTick = [128 256 512 1024 2048 4096 8192 16000 32000];  
%                 smoo_text = sprintf('Smoothing %c',945);
%                 grid_text = sprintf('Bin size (mm)');                 
%             case {'kyadaptive'}
%                 ax1.XTick = [1 2 5 10 20 40 80 160 320 640];
%                 ax1.YTick = [0.5 1 2 4 8];                        
%             case {'ksde'}     
%                 ax1.XTick = [2.5 5 10 20 40 80 160 320 640];
%                 ax1.YTick = [10 20 40 80 160 320]; 
%                 smoo_text = sprintf('Smoothing bandwidth (mm)');
%                 grid_text = sprintf('Bin size (mm)');                  
%             case {'fyhn'}
%                 ax1.XTick = [5 10 20 40 80 160 320 640];
%                 ax1.YTick = [10 20 40 80 160 320 640];                        
%             otherwise
%                 keyboard
%         end               
%         view(0,90)         
%         ylabel(smoo_text);
%         xlabel(grid_text);
% 
%         axc = axes('Units','pixels','Position',[ax1.Position(1)+ax1.Position(3)+30 ax1.Position(2)+0 18 ax1.Position(4)*0.9]); 
%             mat = linspace(ax1.CLim(1),ax1.CLim(2),100)';
%             imagesc(ones(size(mat)),mat,mat,'alphadata',ones(size(mat)).*transp);
%             colormap(axc,error_colormap);
%             axis xy
% 
%             axc.YTick = linspace(ax1.CLim(1),ax1.CLim(2),4);
%             axc.XTick = [];
%             axc.YAxisLocation = 'right';
%             text(0.5,1.1,sprintf('MISE'),'HorizontalAl','left','Units','normalized','FontSize',8);          
% 
% %% #################### Combine figure parts
%         g2 = getframe(fig1);
%         close(fig1);    
%         g3 = [g1.cdata(80:end-220,1:end-520,:); g2.cdata(160:end-570,1:end-520,:)];
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% ################################################################# %% Save the figure
%     % Save the figure    
%     disp(sprintf('\tSaving figure...'))    
%     figname = [fig_dir2 '\fig_lit_rev.png'];
%     [~,~,~] = mkdir(fig_dir2);    
%     %exportgraphics(gcf,figname,'BackgroundColor','w','ContentType','image','Resolution',250);  
%     imwrite(g3,figname);
%     close(gcf) 
% 
% 
% 
    
    
    
    
    
    
    
    
    
    
    
    
    






