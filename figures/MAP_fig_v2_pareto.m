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
    fsiz = 12;

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
    ynow = 690;
    mnames = {'histogram','ash','kadaptive','ksde','fyhn','kyadaptive'};
    dnames = {'Histogram',sprintf('Averaged Shifted Histogram (ASH)'),'Adaptive Smoothing',sprintf('Kernel Smoothed Density Estimate (KSDE)'),'Temporal KSDE','Adaptive Binning'}; 
    ls = {'A','B','C','D','E','F','G'};
    mapidx = 4;
    xsiz = 180;
    xbuff = 220;
    xvec = xnow : (xsiz+xbuff) : 1000;
    ysiz = 200;
    ybuff = 140;
    yvec = ynow : -(ysiz+ybuff) : 0;
    [xx,yy] = meshgrid(xvec,yvec);
    xx = xx';
    yy = yy';

    pls = cell(4,1);
    for mm = 1:length(mnames) % for every method
        mname = mnames{mm};
        sigma_now = '4000';
        time_now = '16'; 

        % main panel showing 16 minute error results
        ax1 = axes('Units','pixels','Position',[xx(mm) yy(mm) xsiz ysiz],'FontSize',fsiz);  
            ah = add_panel_title(ls{mm},sprintf(dnames{mm}),'yoffset',10,'xoffset',0,'width',400,'fontsize',tfonts);

            % load data matrices and prepare them for pareto optimisation 
            cols0 = winter(4);
            shaps = {'.','.','.'};
            mins = NaN(config.npcells,3);
            all_f = cell(1,3);
            all_x = cell(1,3);

            switch mname
                case {'ksde'}
                    ds = [4 16 24];   
                    cols = cols0([1 2 3],:);
                    idxs = [1 2 3];
                otherwise
                    ds = [4 16 64];   
                    cols = cols0([1 2 4],:);         
                    idxs = [1 2 4];
            end   
            switch mname
                case {'histogram'}
                    smoo_text = sprintf('Smoothing %c (mm)',963);
                    grid_text = sprintf('Bin size (mm)');
                case {'ash'}
                    smoo_text = 'Smoothing {\itm}';
                    grid_text = sprintf('Bin size (mm)');                
                case {'ksde'}
                    smoo_text = sprintf('Smoothing bandwidth (mm)');
                    grid_text = sprintf('Bin size (mm)');                
                case {'kyadaptive'}
                    smoo_text = 'Smoothing {\itt} (s)';
                    grid_text = sprintf('Bin size (mm)');         
                case {'kadaptive'}
                    smoo_text = sprintf('Smoothing %c',945);
                    grid_text = sprintf('Bin size (mm)');  
                case {'fyhn'}
                    smoo_text = sprintf('Smoothing %c (mm)',963);
                    grid_text = sprintf('Bin size (mm)');                   
                otherwise
                    keyboard
            end

            fs = [4000 4000 4000];
            for ii = 1:3   
                if ii==1
                    disp(sprintf('\t\t%d mins ',ds(ii)))
                else
                    disp(sprintf('\b| %d mins',ds(ii)))                
                end                
                fname = [scratch_space '\' ename '_' mname '_sigma' num2str(fs(ii)) '_duration' num2str(ds(ii)) '_gamultiobj.mat'];
                if exist(fname,'file') && 1 && ~override_pareto
                    load(fname,'x','fval','balance_solution1','balance_solution2','balance_solution3','minerr_solution');

                    a = fval(:,1);
                    b = fval(:,2); 
                    all_f{ii} = [a b];
                    all_x{ii} = x;

                    % plot the pareto optimal results       
                    plot(a,b,'Color',cols(ii,:),'Marker',shaps{ii},'LineStyle','none','MarkerFaceColor',cols(ii,:)); hold on;

                    idx = convhull(a,b);    
                    chull = [a(idx) b(idx)];
                    chull_scaled = chull ./ max(chull,[],1,'omitnan');
                    dst = pdist2(chull_scaled,[max(chull_scaled(:,1)) min(chull_scaled(:,2))]); % bottom right
                    [~,cut_point2] = min(dst);                           
                    dst = pdist2(chull_scaled,[min(chull_scaled(:,1)) max(chull_scaled(:,2))]); % top left
                    [~,cut_point1] = min(dst);  
                    if cut_point1>cut_point2
                        chull = [chull(cut_point1:end,:); chull(1:cut_point2,:); NaN(1,2)];                            
                    else
                        chull = [chull(cut_point1:cut_point2,:); NaN(1,2)];                            
                    end
                    pls{idxs(ii)} = plot(chull(:,1),chull(:,2),'Color',cols(ii,:),'LineWidth',1.5,'Marker',shaps{ii},'MarkerFaceColor',cols(ii,:),'MarkerSize',10);                    

                    pb1 = plot3(balance_solution1(3),balance_solution1(4),1000,'Color','k','Marker','o','MarkerSize',10,'LineStyle','none');            

                    continue
                else
                    load([scratch_space '\' ename '_' mname '_sigma' num2str(fs(ii)) '_duration' num2str(ds(ii)) '_datamats.mat'],'data_matrix_cputime','data_matrix_fields','data_matrix_error','data_matrix_missing','data_matrix_params')   
                end
                data_matrix_error = data_matrix_error(:,:,:,mapidx);
                bmat = double(data_matrix_params(:,:,1)); % binsizes used when generating maps
                smat = double(data_matrix_params(:,:,2)); % smoothing used when generating maps

                cmap1 = double( mean(data_matrix_error,3,'omitnan') );    
                switch mname
                    case {'fyhn'}
                        tmap1 = double( mean(data_matrix_cputime(:,:,1:10),3,'omitnan') ) .* 1000;                            
                    otherwise
                        tmap1 = double( mean(data_matrix_cputime,3,'omitnan') ) .* 1000;
                end
                fmap1 = abs( double( mean(data_matrix_fields,3,'omitnan') ) );                
                emap1 = abs( double( mean(data_matrix_missing(:,:,:,1),3,'omitnan') ) );      

                % actual functions that will be used by gamultiobj               
                meth = 'nearest';
                f1 = @(xy) interp2(bmat,smat,cmap1,xy(:,1),xy(:,2),meth,inf); % MSE
                f2 = @(xy) interp2(bmat,smat,tmap1,xy(:,1),xy(:,2),meth,inf); % time
                f3 = @(xy) interp2(bmat,smat,emap1,xy(:,1),xy(:,2),meth,inf); % missing data (NaNs)
                f4 = @(xy) interp2(bmat,smat,fmap1,xy(:,1),xy(:,2),meth,inf); % field count error
                fx = @(xy) [f1(xy) f2(xy) f3(xy) f4(xy)];

                rng(999); % for reproducibility
                options = optimoptions('gamultiobj','Display','off','PopulationSize',250);

                % first find the pareto front when minimising error and time only
                % this result will form the 2D plot in the figure
                % fit curve to objective pareto front
                [x,fval] = gamultiobj(fx,2,[],[],[],[],[min(bmat(:)) min(smat(:))],[max(bmat(:)) max(smat(:))],[],options);
                a = fval(:,1);
                b = fval(:,2); 
                all_f{ii} = [a b];
                all_x{ii} = x;

                a2 = a./max(a(:));
                b2 = b./max(b(:));
                utopia_point = [min(a2) min(b2)];

                dom = false(size(a));                    
                for kk = 1:length(a)
                    log_idx = (1:length(a))'==kk;
                    dom(kk) = ~any( a(~log_idx)<a(kk) & b(~log_idx)<b(kk) );   
                end        
                doms = [x(dom,:) fval(dom,:)];

                d = pdist2(utopia_point,[a2(dom) b2(dom)],'euclidean');
                [~,midx] = min(d(:),[],1,'omitnan');                                
                balance_solution1 = doms(midx,:); 

                [~,midx] = min(a,[],1,'omitnan');                
                balance_solution2 = [x(midx,:) fval(midx,:)];                 

                [~,midx] = min(b,[],1,'omitnan');                
                balance_solution3 = [x(midx,:) fval(midx,:)];                 

                mat1 = mean(data_matrix_error,3,'double','omitnan');                    
                i = find(mat1==min(mat1(:),[],'omitnan'),1,'first');
                minerr_solution = [bmat(i),smat(i),mat1(i),tmap1(i)];                 
                save([scratch_space '\' ename '_' mname '_sigma' num2str(fs(ii)) '_duration' num2str(ds(ii)) '_gamultiobj.mat'],'x','fval','balance_solution1','balance_solution2','balance_solution3','minerr_solution');

                % plot the pareto optimal results    
                plot(a,b,'Color',cols(ii,:),'Marker',shaps{ii},'LineStyle','none','MarkerFaceColor',cols(ii,:)); hold on; 
                idx = convhull(a,b);    
                chull = [a(idx) b(idx)];
                chull_scaled = chull ./ max(chull,[],1,'omitnan');
                dst = pdist2(chull_scaled,[max(chull_scaled(:,1)) min(chull_scaled(:,2))]); % bottom right
                [~,cut_point2] = min(dst);                           
                dst = pdist2(chull_scaled,[min(chull_scaled(:,1)) max(chull_scaled(:,2))]); % top left
                [~,cut_point1] = min(dst);  
                if cut_point1>cut_point2
                    chull = [chull(cut_point1:end,:); chull(1:cut_point2,:); NaN(1,2)];                            
                else
                    chull = [chull(cut_point1:cut_point2,:); NaN(1,2)];                            
                end
                pls{idxs(ii)} = plot(chull(:,1),chull(:,2),'Color',cols(ii,:),'LineWidth',1.5);                    
                
                plot3(balance_solution1(3),balance_solution1(4),1000,'Color','k','Marker','o','MarkerSize',10,'LineStyle','none');            

            end     
            xlabel(sprintf('MISE'), 'Interpreter', 'none' );
            ylabel(sprintf('Computation time (ms)'), 'Interpreter', 'none' );
            ax = gca;    
            ax.YLim(1) = 0;
            ax.XLim(1) = 0;
            ax.XScale = 'log';
            ax.YScale = 'log';
            grid off
            box off
            axis 'auto'
            switch mname
                case {'histogram'}
                    ax.YTick = [1 2 4 8 16 32 64];   
                    ax.YLim(1) = 1;
                case {'ash'}
                    ax.YTick = [1 2 4 8 16 32 64];   
                    ax.YLim(1) = 0.5;                    
                otherwise
                    ax.YTick = [10^-8 10^-7 10^-6 10^-5 10^-4 10^-3 10^-2 10^-1 10^0 10^1 10^2 10^3 10^4 10^5 10^6 10^7 10^8 10^9 10^10];
            end
            % ax.XTick = logspace(-13, -10, 20);
            % xtickformat('%0.1g');
            ax.XTick = [0.60 1 2 4 8 16 32 64 100 200 400 800 1600 3200 6400] .* 10^-12;
            ax.XLim = [0.6.*10^-12 16.*10^-12];

            labs_now = {'4 mins Pareto front','16 mins','24 mins','64 mins'};   
            if mm==4
                [~,leg] = legendflex([pls{1} pls{2} pls{3} pls{4}] ,labs_now,'anchor',{'s','s'},'nrow',1,'box','off','buffer',[50,-100],'xscale',.5,'fontsize',9); 
                % leg(5).MarkerSize = 12;
                % leg(7).MarkerSize = 12;
                % leg(9).MarkerSize = 12;
            end
            % if mm==1
                text(0,1.02,sprintf('Pareto-front'),'Units','normalized','HorizontalAlignment','left','VerticalAlignment','bottom','FontSize',9)
            % end

        % inset showing 4 minute error results
        % MAP_fig_1_multi needs to be run for this to work
        % collect regression results
        fs = [4000 8000 16000];
        rfs = [127 179 253];
        switch mname
            case {'ksde'}
                ds = [4 16 24];
            otherwise
                ds = [4 16 64];
        end          
        dat = [];
        for f = 1:length(fs)
            for d = 1:length(ds)
                nme = [ename '_' mname '_sigma' num2str(fs(f)) '_duration' num2str(ds(d)) '_gamultiobj.mat'];
                if exist(nme,'file')
                    load(nme);
                    dat = [dat; rfs(f) ds(d) balance_solution1(1:2) minerr_solution(1:2)];
                end
            end
        end


        disp(sprintf('\tPareto-optimal regression results...'))   
        [b,s,r] = mvregress([ones(size(dat,1),1),dat(:,1:2)],dat(:,3:4)); 
        f = @(xy) [b(1,1) + xy(:,1).*b(2,1) + xy(:,2).*b(3,1), b(1,2) + xy(:,1).*b(2,2) + xy(:,2).*b(3,2)]; % [field size, duration] > [bin size, smoothing]

        colmapnow = magma(256);

        % inset showing 4 minute error results
        ax2 = axes('Units','pixels','Position',[ax1.Position(1)+ax1.Position(3)+35 (ax1.Position(2)+ax1.Position(4))-(ax1.Position(4).*0.4) ax1.Position(3).*0.4 ax1.Position(4).*0.4]);  
            disp(sprintf('\t\t...bin size plot'))

            fs = 125 : 1 : 300;
            switch mname
                case {'ksde'}
                    du = 4 : 1 : 32;
                otherwise
                    du = 4 : 1 : 64;
            end                 

            [xq,yq] = meshgrid(fs,du);
            fx = f([xq(:),yq(:)]);
            bi = reshape(fx(:,1),size(xq)); %                
            sm = reshape(fx(:,2),size(xq)); %       

            im = surf(fs,du,bi,'EdgeColor','none'); hold on;        
            colormap(gca,colmapnow) 

            axis ij
            ax2.XScale = 'log';
            ax2.YScale = 'log';
            ax2.XLim = [min(fs) max(fs)];
            ax2.YLim = [min(du) max(du)];     
            ax2.XTick = [125:50:300];
            switch mname
                case {'ksde'}
                    ax2.YTick = [4 8 16 24 32];
                otherwise
                    ax2.YTick = [4 8 16 32 64];
            end                 
            view(0,90)
            ylabel(sprintf('Duration (mins)'));
            % xlabel(sprintf('Avg. field radius (mm)'));
            % text(0,1.05,grid_text,'Units','normalized','HorizontalAl','left','FontSize',8)        
            % if mm==1
                text(0,1.1,sprintf('Regression fit'),'Units','normalized','HorizontalAlignment','left','VerticalAlignment','bottom','FontSize',9)
            % end

            % colorbar
            axc = axes('Units','pixels','Position',[ax2.Position(1)+ax2.Position(3)+10 ax2.Position(2) 12 ax2.Position(4)]);
                mat = (linspace(ax2.CLim(1),ax2.CLim(2),100))';
                imagesc('YData',mat,'XData',ones(size(mat)),'CData',mat);
                colormap(axc,colmapnow);
                axis xy
                axc.YLim = ax2.CLim;     
                axc.XTick = [];
                axc.YAxisLocation = 'right';
                ylabel(grid_text)

        % inset showing 64 minute error results
        ax3 = axes('Units','pixels','Position',[ax1.Position(1)+ax1.Position(3)+35 ax1.Position(2) ax1.Position(3).*0.4 ax1.Position(4).*0.4]);     
            disp(sprintf('\t\t...smoothing plot'))

            im = surf(fs,du,sm,'EdgeColor','none'); hold on;        
            colormap(gca,colmapnow) 

            axis ij
            ax3.XScale = 'log';
            ax3.YScale = 'log';
            ax3.XLim = [min(fs) max(fs)];
            ax3.YLim = [min(du) max(du)];     
            ax3.XTick = [125:50:300];
            ax3.YTick = [4 8 16 32 64];
            view(0,90)
            xlabel(sprintf('Avg. field radius (mm)'));
            ylabel(sprintf('Duration (mins)'));            
            % text(0,1.05,smoo_text,'Units','normalized','HorizontalAl','left','FontSize',8)        

            % colorbar
            axc = axes('Units','pixels','Position',[ax3.Position(1)+ax3.Position(3)+10 ax3.Position(2) 12 ax3.Position(4)]);
                mat = (linspace(ax3.CLim(1),ax3.CLim(2),100))';
                imagesc('YData',mat,'XData',ones(size(mat)),'CData',mat);
                colormap(axc,colmapnow);
                axis xy
                axc.YLim = ax3.CLim;     
                axc.XTick = [];
                axc.YAxisLocation = 'right';
                ylabel(smoo_text)

            %     if mm==2
            % keyboard
            %     end
    end
% return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ################################################################# %% Save the figure
    % Save the figure    
    disp(sprintf('\tSaving figure...'))    
    figname = [fig_dir2 '\Fig 9.png'];
    [~,~,~] = mkdir(fig_dir2);    
    exportgraphics(fig1,figname,'BackgroundColor','w','ContentType','image','Resolution',350);  
    close(fig1)    
    
    


    
    
    
    
  
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    