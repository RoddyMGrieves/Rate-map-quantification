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
    ynow = 780;
    mnames = {'histogram','ash','kadaptive','ksde','fyhn','kyadaptive'};
    dnames = {'Histogram','ASH','Adaptive smoothing','KSDE','tKSDE','Adaptive binning'}; 
    ls = {'a','b','c','d','e','f','g'};
    mapidx = 4;
    xsiz = 200;
    xbuff = 160;
    xvec = xnow : (xsiz+xbuff) : 1000;
    ysiz = 200;
    ybuff = 100;
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

        % main panel showing 16 minute error results
        ax1 = axes('Units','pixels','Position',[xx(mm) yy(mm) xsiz ysiz]);  
            load([scratch_space '\' ename '_' mname '_sigma' num2str(sigma_now) '_duration' num2str(time_now) '_datamats.mat'],'data_matrix_error','data_matrix_params')  
            bmat = data_matrix_params(:,:,1); % binsizes used when generating maps
            smat = data_matrix_params(:,:,2); % smoothing used when generating maps 

            switch mname
                case {'histogram'}
                    ah = add_panel_title(ls{mm},sprintf(dnames{mm}),'yoffset',24,'xoffset',0,'width',400);
                otherwise
                    ah = add_panel_title(ls{mm},sprintf(dnames{mm}),'yoffset',5,'xoffset',0,'width',400);
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
            switch mname
                case {'histogram'}
                    text(0,1.15,'16 minutes coverage','HorizontalAl','left','Units','normalized','FontSize',8);                
                otherwise
                    text(0,1.04,'16 minutes coverage','HorizontalAl','left','Units','normalized','FontSize',8);
            end            
            ylabel(smoo_text);
            xlabel(grid_text);

            % need to run MAP_fig_1_multi for this next file to exist
            load([scratch_space '\' ename '_' mname '_sigma' num2str(sigma_now) '_duration' num2str(time_now) '_gamultiobj.mat'],'balance_solution2','balance_solution1','minerr_solution','x')
            pb100 = plot3(balance_solution1(1),balance_solution1(2),1000,'ko','MarkerSize',8,'LineStyle','none','MarkerFaceColor','k');            
            pc00 = plot3(minerr_solution(1),minerr_solution(2),1000,'kv','MarkerSize',8,'LineStyle','none','MarkerFaceColor','k');          

            if mm==5
                [~,leg] = legendflex([pb100 pc00]...
                    ,{'Balanced solution','Minimum error solution'}...
                    ,'anchor',{'s','s'},'ncol',1,'box','off','buffer',[-130,-110],'xscale',.5,'fontsize',9);  
            end

            if mm == 4 
                axc = axes('Units','pixels','Position',[ax1.Position(1)-0 ax1.Position(2)-90 ax1.Position(3) 10]); 
                    mat = linspace(ax1.CLim(1),ax1.CLim(2),100);
                    imagesc('YData',ones(size(mat)),'XData',mat,'CData',mat);
                    colormap(axc,error_colormap);
                    axis xy
    
                    axc.YTick = [];
                    axc.XTickLabelRotation = 0; 
                    axc.XLim = ax1.CLim;
                    text(0.5,1.7,'Map error (MISE)','FontSize',8,'HorizontalAl','center','Units','normalized')                    
            end

            if strcmp(mname,'histogram')
                ax1b = axes('Units','pixels','Position',[ax1.Position(1) ax1.Position(2)+ax1.Position(4)+6 ax1.Position(3) 15]);  
                    imagesc(mat1(1,:)); hold on;
                    colormap(gca,error_colormap) 
                    set(gca,'ColorScale','log')
    
                    axis ij
                    caxis(ax1.CLim);
                    ax1b.YTick = 1;
                    ax1b.YTickLabel = {'none'};
                    ax1b.XTick = [];
                    axis off
                    text(-0.03,0.5,sprintf('none'),'FontSize',8,'HorizontalAl','right','Units','normalized')      
            end     
    
            % inset showing 4 minute error results
            ax2 = axes('Units','pixels','Position',[ax1.Position(1)+ax1.Position(3)+10 (ax1.Position(2)+ax1.Position(4))-(ax1.Position(4).*0.4) ax1.Position(3).*0.4 ax1.Position(4).*0.4]);   
                time_now = '4';                     
                load([scratch_space '\' ename '_' mname '_sigma' num2str(sigma_now) '_duration' num2str(time_now) '_datamats.mat'],'data_matrix_error','data_matrix_params')    
                dme = data_matrix_error(:,:,:,mapidx);
                bmat = data_matrix_params(:,:,1); % binsizes used when generating maps
                smat = data_matrix_params(:,:,2); % smoothing used when generating maps
    
                mat2 = mean(dme,3,'double','omitnan');    
                Vq2 = interp2(bmat,smat,mat2,Xq,Yq,interp_method);
    
                im = surf(Xq,Yq,Vq2,'EdgeColor','none'); hold on;
                colormap(gca,error_colormap) 
                set(gca,'ColorScale','log')
                caxis(ax1.CLim)            
                axis ij off
                ax2.XScale = 'log';
                ax2.YScale = 'log'; 
                ax2.XLim = ax1.XLim;
                ax2.YLim = ax1.YLim;          
                view(0,90)
                switch mname
                    case {'histogram'}
                        text(0.01,1.25,'4 minutes','HorizontalAl','left','Units','normalized','FontSize',8);                
                    otherwise
                        text(0.01,1.08,'4 minutes','HorizontalAl','left','Units','normalized','FontSize',8);
                end     

                % need to run MAP_fig_1_multi for this next file to exist
                load([scratch_space '\' ename '_' mname '_sigma' num2str(sigma_now) '_duration' num2str(time_now) '_gamultiobj.mat'],'balance_solution2','balance_solution1','minerr_solution','x')
                pb1 = plot3(balance_solution1(1),balance_solution1(2),1000,'ko','MarkerSize',5,'LineStyle','none','MarkerFaceColor','k');            
                pc = plot3(minerr_solution(1),minerr_solution(2),1000,'kv','MarkerSize',5,'LineStyle','none','MarkerFaceColor','k'); 

                if strcmp(mname,'histogram')
                    ax2b = axes('Units','pixels','Position',[ax2.Position(1) ax2.Position(2)+ax2.Position(4)+4 ax2.Position(3) 10]);  
                        imagesc(mat2(1,:)); hold on;
                        colormap(gca,error_colormap) 
                        set(gca,'ColorScale','log')
    
                        axis ij
                        caxis(ax2.CLim);
                        ax2b.YTick = 1;
                        ax2b.YTickLabel = {'none'};
                        ax2b.XTick = [];
                        axis off     
                end
    
            % inset showing 64 minute error results
            ax3 = axes('Units','pixels','Position',[ax1.Position(1)+ax1.Position(3)+10 ax1.Position(2) ax1.Position(3).*0.4 ax1.Position(4).*0.4]);    
                switch mname
                    case {'ksde'}
                        time_now = '24';                     
                    otherwise
                        time_now = '64';                     
                end              
                load([scratch_space '\' ename '_' mname '_sigma' num2str(sigma_now) '_duration' num2str(time_now) '_datamats.mat'],'data_matrix_error','data_matrix_params') 
                dme = data_matrix_error(:,:,:,mapidx);            
                bmat = data_matrix_params(:,:,1); % binsizes used when generating maps
                smat = data_matrix_params(:,:,2); % smoothing used when generating maps
    
                mat3 = mean(dme,3,'double','omitnan');    
                Vq3 = interp2(bmat,smat,mat3,Xq,Yq,interp_method);
    
                im = surf(Xq,Yq,Vq3,'EdgeColor','none'); hold on;
                colormap(gca,error_colormap) 
                set(gca,'ColorScale','log')
                caxis(ax1.CLim)
                axis ij off
                ax3.XScale = 'log';
                ax3.YScale = 'log';
                ax3.XLim = ax1.XLim;
                ax3.YLim = ax1.YLim;          
                view(0,90)
                switch mname
                    case {'histogram'}
                        text(0.01,1.25,'64 minutes','HorizontalAl','left','Units','normalized','FontSize',8);   
                    case {'ksde'}
                        text(0.01,1.08,'24 minutes','HorizontalAl','left','Units','normalized','FontSize',8);                    
                    otherwise
                        text(0.01,1.08,'64 minutes','HorizontalAl','left','Units','normalized','FontSize',8);
                end 

                % need to run MAP_fig_1_multi for this next file to exist
                load([scratch_space '\' ename '_' mname '_sigma' num2str(sigma_now) '_duration' num2str(time_now) '_gamultiobj.mat'],'balance_solution2','balance_solution1','minerr_solution','x')
                pb1 = plot3(balance_solution1(1),balance_solution1(2),1000,'ko','MarkerSize',5,'LineStyle','none','MarkerFaceColor','k');            
                pc = plot3(minerr_solution(1),minerr_solution(2),1000,'kv','MarkerSize',5,'LineStyle','none','MarkerFaceColor','k'); 

                if strcmp(mname,'histogram')
                    ax3b = axes('Units','pixels','Position',[ax3.Position(1) ax3.Position(2)+ax3.Position(4)+4 ax3.Position(3) 10]);  
                        imagesc(mat3(1,:)); hold on;
                        colormap(gca,error_colormap) 
                        set(gca,'ColorScale','log')
    
                        axis ij
                        caxis(ax3.CLim);
                        ax3b.YTick = 1;
                        ax3b.YTickLabel = {'none'};
                        ax3b.XTick = [];
                        axis off            
                end     
    end
% return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ################################################################# %% Save the figure
    % Save the figure    
    disp(sprintf('\tSaving figure...'))    
    figname = [fig_dir2 '\v2_fig_error_matrix.png'];
    [~,~,~] = mkdir(fig_dir2);    
    exportgraphics(fig1,figname,'BackgroundColor','w','ContentType','image','Resolution',250);  
    close(fig1)    
    
    


  
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    