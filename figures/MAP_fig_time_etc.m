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


    ename = 'arena120cm'; 
    override_pareto = 0;
    override_maps = 0;
    rmap_colormap = 'turbo';
    tfonts = [20 12];
    fsiz = 9;

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

    % method settings
    xnow = 70;
    ynow = 730;
    mnames = {'histogram','ash','kadaptive','ksde','fyhn','kyadaptive'};
    dnames = {'Histogram','Averaged Shifted Histogram (ASH)','Adaptive Smoothing','Kernel Smoothed Density Estimate (KSDE)','Temporal KSDE','Adaptive Binning'}; 
    ls = {'A','B','C','D','E','F','G'};
    mapidx = 4;
    xsiz = 200;
    xbuff = 240;
    xvec = xnow : (xsiz+xbuff) : 1000;
    ysiz = 200;
    ybuff = 120;
    yvec = ynow : -(ysiz+ybuff) : 0;
    [xx,yy] = meshgrid(xvec,yvec);
    xx = xx';
    yy = yy';

    % plot settings
    n_color_levels = 64;
    error_colormap = flipud(inferno(n_color_levels));    
    clims = [1*10^-12 1*10^-11];   
    interp_method = 'nearest';  
  

    for mm = 1:length(mnames) % for every method
        mname = mnames{mm};
        sigma_now = '4000';
        time_now = '16'; 

        % main panel showing place field error 
        ax1 = axes('Units','pixels','Position',[xx(mm) yy(mm) xsiz ysiz],'FontSize',fsiz);  
            load([scratch_space '\' ename '_' mname '_sigma' num2str(sigma_now) '_duration' num2str(time_now) '_datamats.mat'],'data_matrix_cputime','data_matrix_fields','data_matrix_error','data_matrix_missing','data_matrix_params')  
            bmat = data_matrix_params(:,:,1); % binsizes used when generating maps
            smat = data_matrix_params(:,:,2); % smoothing used when generating maps 

            switch mname
                case {'histogram'}
                    ah = add_panel_title(ls{mm},sprintf(dnames{mm}),'yoffset',28,'xoffset',0,'width',400,'fontsize',tfonts);
                otherwise
                    ah = add_panel_title(ls{mm},sprintf(dnames{mm}),'yoffset',5,'xoffset',0,'width',400,'fontsize',tfonts);
            end
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
            switch mname
                case {'histogram'}
                    colplot = [0 50; 0 1000; 0 0.25; 0 1]; % time, memory, missing, field error
                    colname = {'0ms','50ms';'0kb',sprintf('1{\\times}10^{3}kb');'0%','25%';'0','1'};
                case {'kyadaptive'}
                    colplot = [0 2000; 0 1000; 0 0.01; 0 0.5];
                    colname = {'0ms','2000ms';'0kb',sprintf('1{\\times}10^{3}kb');'0%','1%';'0','0.5'};  
                case {'kadaptive'}
                    colplot = [0 1000; 0 1000; 0 0.01; 0 0.5]; % time, memory, missing, field error
                    colname = {'0ms','1000ms';'0kb',sprintf('1{\\times}10^{3}kb');'0%','1%';'0 fields','0.5 fields'};                
                case {'ash'}
                    colplot = [0 50; 0 1000; 0 0.25; 0 0.5];
                    colname = {'0ms','50ms';'0kb',sprintf('1{\\times}10^{3}kb');'0%','25%';'0','0.5'};
                case {'ksde'}
                    colplot = [0 2000; 0 1000; 0 0.25; 0 0.5];
                    colname = {'0 ms','2000 ms';'0 kb',sprintf('1{\\times}10^{3}kb');'0%','25%';'0','0.5'};   
                case {'fyhn'}
                    colplot = [0 2000; 0 1000; 0 0.25; 0 0.5]; % time, memory, missing, field error
                    colname = {'0ms','2000ms';'0kb',sprintf('1{\\times}10^{3}kb');'0%','25%';'0','0.5'};                    
                otherwise
                    keyboard
            end         
            plotname = {'Computation time (ms) ';'Disk space (kb) ';'Proportion of map empty ';'Place field detection error (fields) '};
            cmapnow = parula(256);             
            [Xq,Yq] = meshgrid(new_bvec,new_svec);   

            dat = data_matrix_fields(:,:,:);
            mat1 = mean(dat,3,'double','omitnan');    
            bmat = data_matrix_params(:,:,1); % binsizes used when generating maps
            smat = data_matrix_params(:,:,2); % smoothing used when generating maps                

            dme = data_matrix_error(:,:,:,4);
            mat_e = mean(dme,3,'double','omitnan');  
            [r,p] = corr(mat1(:),mat_e(:),'rows','pairwise','type','Pearson');
            disp(sprintf('\tcorr between error,fields %s: %.3f, p = %.4f',mname,r,p))

            Vq1 = interp2(bmat,smat,mat1,Xq,Yq,interp_method);
            im = surf(Xq,Yq,Vq1,'EdgeColor','none'); hold on;
            colormap(gca,cmapnow) 
            set(gca,'ColorScale','log')

            axis ij
            ax1.XScale = 'log';
            ax1.YScale = 'log';  
            ax1.XLim = [min(new_bvec) max(new_bvec)];
            ax1.YLim = [min(new_svec) max(new_svec)]; 
            ax1.FontSize = fsiz;
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

            caxis(colplot(4,:))
            switch mname
                case {'histogram'}
                    text(0.01,1.18,plotname{4},'HorizontalAl','left','Units','normalized','FontSize',fsiz);        
                otherwise
                    text(0.01,1.05,plotname{4},'HorizontalAl','left','Units','normalized','FontSize',fsiz);        
            end           
            ylabel(smoo_text);
            xlabel(grid_text);

            % need to run MAP_fig_1_multi for this next file to exist
            load([scratch_space '\' ename '_' mname '_sigma' num2str(sigma_now) '_duration' num2str(time_now) '_gamultiobj.mat'],'balance_solution2','balance_solution1','minerr_solution','x')
            pb = plot3(balance_solution1(1),balance_solution1(2),10000,'ko','MarkerSize',8,'LineStyle','none','MarkerFaceColor','k');            
            pc = plot3(minerr_solution(1),minerr_solution(2),10000,'kv','MarkerSize',8,'LineStyle','none','MarkerFaceColor','k');          

            if mm==4
                [~,leg] = legendflex([pb pc]...
                    ,{'Balanced solution','Minimum error solution'}...
                    ,'anchor',{'n','n'},'ncol',1,'box','off','buffer',[-30,-260],'xscale',.5,'fontsize',fsiz);  
            end

            if strcmp(mname,'histogram')
                ax1b = axes('Units','pixels','Position',[ax1.Position(1) ax1.Position(2)+ax1.Position(4)+6 ax1.Position(3) 15]);  
                    imagesc(mat1(1,:)); hold on;
                    colormap(gca,cmapnow) 
                    set(gca,'ColorScale','log')
    
                    axis ij
                    caxis(ax1.CLim);
                    ax1b.YTick = 1;
                    ax1b.YTickLabel = {'none'};
                    ax1b.XTick = [];
                    axis off
                    text(-0.03,0.5,sprintf('none'),'FontSize',fsiz,'HorizontalAl','right','Units','normalized')      
            end     

            axc = axes('Units','pixels','Position',[ax1.Position(1)+ax1.Position(3)+10 ax1.Position(2) 10 ax1.Position(4)]); 
                mat = linspace(ax1.CLim(1),ax1.CLim(2),100);
                imagesc('XData',ones(size(mat)),'YData',mat,'CData',mat');
                colormap(axc,cmapnow);
                axis xy

                axc.XTick = [];
                axc.YTickLabelRotation = 0; 
                axc.YLim = ax1.CLim;
                axc.YAxisLocation = 'right';                    

%%%%%%%%%%%%%% computation time
            ax2 = axes('Units','pixels','Position',[ax1.Position(1)+ax1.Position(3)+50 (ax1.Position(2)+ax1.Position(4))-(ax1.Position(4).*0.4) ax1.Position(3).*0.4 ax1.Position(4).*0.4],'FontSize',fsiz);   
                switch mname
                    case {'fyhn'}
                        dat = data_matrix_cputime(:,:,1:10).*1000;
                    otherwise
                        dat = data_matrix_cputime(:,:,:).*1000;
                end 
                mat1 = mean(dat,3,'double','omitnan');    
                bmat = data_matrix_params(:,:,1); % binsizes used when generating maps
                smat = data_matrix_params(:,:,2); % smoothing used when generating maps                
    
                Vq1 = interp2(bmat,smat,mat1,Xq,Yq,interp_method);
                im = surf(Xq,Yq,Vq1,'EdgeColor','none'); hold on;
                colormap(gca,cmapnow) 
                set(gca,'ColorScale','log')
    
                axis ij
                ax2.XScale = 'log';
                ax2.YScale = 'log';  
                ax2.XLim = [min(new_bvec) max(new_bvec)];
                ax2.YLim = [min(new_svec) max(new_svec)];             
                ax2.YTick = [];    
                ax2.XTick = [];               
                view(0,90)
    
                caxis(colplot(1,:))
                switch mname
                    case {'histogram'}
                        text(0.01,1.40,plotname{1},'HorizontalAl','left','Units','normalized','FontSize',8);        
                    otherwise
                        text(0.01,1.15,plotname{1},'HorizontalAl','left','Units','normalized','FontSize',8);        
                end           
    
                % need to run MAP_fig_1_multi for this next file to exist
                load([scratch_space '\' ename '_' mname '_sigma' num2str(sigma_now) '_duration' num2str(time_now) '_gamultiobj.mat'],'balance_solution2','balance_solution1','minerr_solution','x')
                pb100 = plot3(balance_solution1(1),balance_solution1(2),100000,'ko','MarkerSize',8,'LineStyle','none','MarkerFaceColor','k');            
                pc00 = plot3(minerr_solution(1),minerr_solution(2),100000,'kv','MarkerSize',8,'LineStyle','none','MarkerFaceColor','k');          
    
                if strcmp(mname,'histogram')
                    ax2b = axes('Units','pixels','Position',[ax2.Position(1) ax2.Position(2)+ax2.Position(4)+6 ax2.Position(3) 15]);  
                        imagesc(mat1(1,:)); hold on;
                        colormap(gca,cmapnow) 
                        set(gca,'ColorScale','log')
        
                        axis ij
                        caxis(ax2.CLim);
                        ax2b.YTick = 1;
                        ax2b.YTickLabel = {'none'};
                        ax2b.XTick = [];
                        axis off
                end   
    
            axc = axes('Units','pixels','Position',[ax2.Position(1)+ax2.Position(3)+10 ax2.Position(2) 10 ax2.Position(4)]); 
                mat = linspace(ax2.CLim(1),ax2.CLim(2),100);
                imagesc('XData',ones(size(mat)),'YData',mat,'CData',mat');
                colormap(axc,cmapnow);
                axis xy

                axc.XTick = [];
                axc.YTickLabelRotation = 0; 
                axc.YLim = ax2.CLim;
                axc.YAxisLocation = 'right';

%%%%%%%%%%%%%%%%%%%% missing pixels
            ax3 = axes('Units','pixels','Position',[ax1.Position(1)+ax1.Position(3)+50 ax1.Position(2)-10 ax1.Position(3).*0.4 ax1.Position(4).*0.4],'FontSize',fsiz);    
                dat = data_matrix_missing(:,:,:,1);
                mat1 = mean(dat,3,'double','omitnan');    
                bmat = data_matrix_params(:,:,1); % binsizes used when generating maps
                smat = data_matrix_params(:,:,2); % smoothing used when generating maps                
    
                Vq1 = interp2(bmat,smat,mat1,Xq,Yq,interp_method);
                im = surf(Xq,Yq,Vq1,'EdgeColor','none'); hold on;
                colormap(gca,cmapnow) 
                set(gca,'ColorScale','log')
    
                axis ij
                ax3.XScale = 'log';
                ax3.YScale = 'log';  
                ax3.XLim = [min(new_bvec) max(new_bvec)];
                ax3.YLim = [min(new_svec) max(new_svec)];                   
                ax3.YTick = [];    
                ax3.XTick = [];
                view(0,90)
    
                caxis(colplot(3,:))
                switch mname
                    case {'histogram'}
                        text(0.01,1.40,plotname{3},'HorizontalAl','left','Units','normalized','FontSize',8);        
                    otherwise
                        text(0.01,1.15,plotname{3},'HorizontalAl','left','Units','normalized','FontSize',8);        
                end           
    
                % need to run MAP_fig_1_multi for this next file to exist
                load([scratch_space '\' ename '_' mname '_sigma' num2str(sigma_now) '_duration' num2str(time_now) '_gamultiobj.mat'],'balance_solution2','balance_solution1','minerr_solution','x')
                pb100 = plot3(balance_solution1(1),balance_solution1(2),100000,'ko','MarkerSize',8,'LineStyle','none','MarkerFaceColor','k');            
                pc00 = plot3(minerr_solution(1),minerr_solution(2),100000,'kv','MarkerSize',8,'LineStyle','none','MarkerFaceColor','k');          
    
                if strcmp(mname,'histogram')
                    ax3b = axes('Units','pixels','Position',[ax3.Position(1) ax3.Position(2)+ax3.Position(4)+6 ax3.Position(3) 15]);  
                        imagesc(mat1(1,:)); hold on;
                        colormap(gca,cmapnow) 
                        set(gca,'ColorScale','log')
        
                        axis ij
                        caxis(ax3.CLim);
                        ax3b.YTick = 1;
                        ax3b.YTickLabel = {'none'};
                        ax3b.XTick = [];
                        axis off
                end   

            axc = axes('Units','pixels','Position',[ax3.Position(1)+ax3.Position(3)+10 ax3.Position(2) 10 ax3.Position(4)]); 
                mat = linspace(ax3.CLim(1),ax3.CLim(2),100);
                imagesc('XData',ones(size(mat)),'YData',mat,'CData',mat');
                colormap(axc,cmapnow);
                axis xy

                axc.XTick = [];
                axc.YTickLabelRotation = 0; 
                axc.YLim = ax3.CLim;
                axc.YAxisLocation = 'right';

    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ################################################################# %% Save the figure
    % Save the figure    
    disp(sprintf('\tSaving figure...'))    
    figname = [fig_dir2 '\Fig 8.png'];
    [~,~,~] = mkdir(fig_dir2);    
    exportgraphics(fig1,figname,'BackgroundColor','w','ContentType','image','Resolution',350);  
    close(fig1)    
    
    


  
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    