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
    fig1 = figure('visible','on','Position',[50,60,1000,920]); 
    set(gcf,'InvertHardCopy','off'); % this stops white lines being plotted black      
    set(gcf,'color','w'); % makes the background colour white
    colormap(jet(256)); % to make sure the colormap is not the horrible default one

    % method settings
    xnow = 80;
    ynow = 760;
    mnames = {'histogram','ash','kadaptive','ksde','fyhn','kyadaptive'};    
    dnames = {'Histogram','ASH','Adaptive smoothing','KSDE','tKSDE','Adaptive binning'}; 
    enames = {'MISE','Pearson','Euclidean'};
    mapidx = [4 3 2];
    xsiz = 100;
    xbuff = 20;
    ysiz = 100;
    ybuff = 30;

    xx = [xnow xnow+xsiz+xbuff xnow+xsiz*2+xbuff*2 xnow+xsiz*3+xbuff*3];
    yvec = ynow : -(ysiz+ybuff) : 0;

    % plot settings
    n_color_levels = 64;
    error_colormap = flipud(inferno(n_color_levels));    
    clims = [1*10^-12 1*10^-11; 0.75 1; 0.001 0.01; 0.1 1];
    interp_method = 'nearest';  

    sols = NaN(length(mnames),2,2,3);
  
    for mm = 1:length(mnames) % for every method
        disp(sprintf('\t%s...',mnames{mm}))
        
        for ee = 1:length(mapidx) % for every measure
            mname = mnames{mm};
            sigma_now = '16000';
            time_now = '8'; 
            disp(sprintf('\t\t...%d',mapidx(ee)))
    
            % main panel showing 16 minute error results
            ax1 = axes('Units','pixels','Position',[xx(ee) yvec(mm) xsiz ysiz]);  
                disp(sprintf('\t\t...errormap'))

                load([scratch_space '\' ename '_' mname '_sigma' num2str(sigma_now) '_duration' num2str(time_now) 'normal_datamats.mat'],'data_matrix_error','data_matrix_params')  
                bmat = data_matrix_params(:,:,1); % binsizes used when generating maps
                smat = data_matrix_params(:,:,2); % smoothing used when generating maps 
    
                if mm==1 && ee==1
                    ah = add_panel_title('',sprintf(dnames{mm}),'yoffset',-15,'xoffset',0,'width',400);
                elseif mm>1 && ee==1
                    ah = add_panel_title('',sprintf(dnames{mm}),'yoffset',-15,'xoffset',0,'width',400);
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
                        smoo_text = sprintf('Smoothing (mm)');
                        grid_text = sprintf('Bin size (mm)');                      
                        new_bvec = min(bmat(:)):1:max(bmat(:));
                        new_svec = min(smat(:)):0.5:max(smat(:));                    
                    otherwise
                        keyboard 
                end                          
                [Xq,Yq] = meshgrid(new_bvec,new_svec);     
    
                dme = data_matrix_error(:,:,:,mapidx(ee));
                mat1 = mean(dme,3,'double','omitnan');    
    
                Vq1 = interp2(bmat,smat,mat1,Xq,Yq,interp_method);
                surf(Xq,Yq,Vq1,'EdgeColor','none'); hold on;
                if ee==2 || ee==4
                    colormap(gca,flipud(error_colormap)) 
                else
                    colormap(gca,error_colormap) 
                end
                set(gca,'ColorScale','log')
                caxis(clims(ee,:)) 

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

% if ee==4
% 
% keyboard
% end

                if mm==1 && ee==1
                    text(0,1.35,enames{ee},'HorizontalAl','left','Units','normalized','FontSize',11);
                    text(-0.40,1.38,'a','FontSize',18,'Units','normalized')                    
                elseif mm==1 && ee==2 
                    text(0,1.35,enames{ee},'HorizontalAl','left','Units','normalized','FontSize',11);
                elseif mm==1 && ee==3
                    text(0,1.35,enames{ee},'HorizontalAl','left','Units','normalized','FontSize',11); 
                elseif mm==1 && ee==4
                    text(0,1.35,enames{ee},'HorizontalAl','left','Units','normalized','FontSize',11);                    
                end

                if mm==length(mnames)
                    xlabel(grid_text);    
                else
                    ax1.XTick = [];
                end
                if ee==1
                    ylabel(smoo_text);                    
                else
                    ax1.YTick = [];
                end

%% #################### calculate Pareto front
            fname = [scratch_space '\' ename '_' mname '_sigma' num2str(sigma_now) '_duration' num2str(time_now) 'normal_mapidx' num2str(mapidx(ee)) '_gamultiobj.mat'];
            if 1
                disp(sprintf('\t\t...pareto front'))

                if exist(fname,'file') && 1 && ~override_pareto
                    % do nothing
                else
                    load([scratch_space '\' ename '_' mname '_sigma' num2str(sigma_now) '_duration' num2str(time_now) 'normal_datamats.mat'],'data_matrix_cputime','data_matrix_fields','data_matrix_error','data_matrix_missing','data_matrix_params')   
                
                    disp(sprintf('\t\t\t...calculating'))
                    
                    data_matrix_error = data_matrix_error(:,:,:,mapidx(ee));
                    bmat = double(data_matrix_params(:,:,1)); % binsizes used when generating maps
                    smat = double(data_matrix_params(:,:,2)); % smoothing used when generating maps
    
                    cmap1 = double( mean(data_matrix_error,3,'omitnan') );    
                    if mapidx(ee)==3 || mapidx(ee)==1
                        cmap1 = cmap1 ./ max(cmap1(:),[],'omitnan'); % we want to maximise, not minimise MI and correlation
                        cmap1 = ones(size(cmap1))-cmap1;
                    end
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
    
                    rng(1000); % for reproducibility
                    options = optimoptions('gamultiobj','Display','off','PopulationSize',250);
    
                    % first find the pareto front when minimising error and time only
                    % this result will form the 2D plot in the figure
                    % fit curve to objective pareto front
                    [x,fval] = gamultiobj(fx,2,[],[],[],[],[min(bmat(:)) min(smat(:))],[max(bmat(:)) max(smat(:))],[],options);
                    a = fval(:,1);
                    b = fval(:,2); 
    
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
                    if mapidx(ee)==3 || mapidx(ee)==1
                        i = find(mat1==max(mat1(:),[],'omitnan'),1,'first'); % we want to maximise, not minimise MI and correlation
                    else
                        i = find(mat1==min(mat1(:),[],'omitnan'),1,'first');
                    end                    
                    minerr_solution = [bmat(i),smat(i),mat1(i),tmap1(i)];                 
                    save(fname,'x','fval','balance_solution1','balance_solution2','balance_solution3','minerr_solution');                
                end
                
            end        

%% #################### Plot Solutions onto error map 
            if 1 
                disp(sprintf('\t\t...plotting solutions'))
                
                load(fname,'balance_solution2','balance_solution1','minerr_solution','x');
                % find the 'best' solution (closest to minimising everything)
                pb1 = plot3(balance_solution1(1),balance_solution1(2),1000,'ko','MarkerSize',8,'LineStyle','none');   
                sols(mm,:,1,ee) = [balance_solution1(1),balance_solution1(2)];

                % find the minimal error solution 
                pc = plot3(minerr_solution(1),minerr_solution(2),1000,'kv','MarkerSize',8,'LineStyle','none');          
                sols(mm,:,2,ee) = [minerr_solution(1),minerr_solution(2)];  
            end

%% #################### colorbar
            if mm==length(mnames)
                axc = axes('Units','pixels','Position',[ax1.Position(1) ax1.Position(2)-65 ax1.Position(3) 10]); 
                    mat = linspace(ax1.CLim(1),ax1.CLim(2),100);
                    imagesc([ax1.CLim(1),ax1.CLim(2)],[1 1],mat);
                    colormap(axc,ax1.Colormap);
                    axis xy

                    axc.XTick = linspace(ax1.CLim(1),ax1.CLim(2),4);
                    axc.YTick = [];
                    axc.XAxisLocation = 'bottom';
                    text(0.5,1.5,sprintf(enames{ee}),'FontSize',8,'HorizontalAl','center','Units','normalized')
            end  
        end        
    end

%% #################### Differences in binsize and smoothing
    xnow = 510;
    ynow = 700;
    jit = 0.05;        
    ax1 = axes('Units','pixels','Position',[xnow ynow 160 130]);  
        ah = add_panel_title('b',sprintf('Different error measures'),'yoffset',10,'xoffset',-20,'width',300);                    

        bal_bin = squeeze(sols(:,1,1,:)); % all methods, binsize, balanced, all measures
        bal_smoo = squeeze(sols(:,2,1,:)); % all methods, smoothing, balanced, all measures

        cols = winter(length(mnames));
        bal_bin_y = [];
        for jj = 1:6
            y = (bal_bin(jj,:) - bal_bin(jj,1))./(bal_bin(jj,:) + bal_bin(jj,1));
            bal_bin_y = [bal_bin_y;y];
            plot((2:3)+normrnd(0,jit,[1 2]),y(2:3),'Color',cols(jj,:),'Marker','o'); hold on;
        end

        % ax1.YScale = 'log';
        ax1.YLim = [-1 1];
        ax1.XTick = [];
        ylabel(sprintf('Balanced bin size (mm)\n(a-b)/(a+b)'))
        ax1.XLim = [1.5 3.5];
        text(0,1.15,'Balanced solution bin size','Units','normalized','FontSize',10)
        line(ax1.XLim,[0 0],'Color',[.5 .5 .5])
        box off

        for jj = 2:3
            [H,P,CI] = ttest(bal_bin_y(:,jj),0);
            if P<=.05
                text(jj,ax1.YLim(2).*0.95,'*','FontSize',15,'HorizontalAl','center')
            else
                text(jj,ax1.YLim(2).*1,'n.s.','FontSize',8,'HorizontalAl','center')            
            end
        end
        
    ax2 = axes('Units','pixels','Position',[xnow ynow-180 160 130]);  
        bal_smoo_y = [];    
        for jj = 1:6
            y = (bal_smoo(jj,:) - bal_smoo(jj,1))./(bal_smoo(jj,:) + bal_smoo(jj,1));
            bal_smoo_y = [bal_smoo_y;y];            
            plot((2:3)+normrnd(0,jit,[1 2]),y(2:3),'Color',cols(jj,:),'Marker','o'); hold on;
        end

        ax2.YLim = ax1.YLim;
        ax2.XTick = [];
        ylabel(sprintf('Balanced smoothing (mm)\n(a-b)/(a+b)'))
        ax2.XLim = [1.5 3.5];
        text(0,1.15,'Balanced solution smoothing','Units','normalized','FontSize',10)
        line(ax1.XLim,[0 0],'Color',[.5 .5 .5])
        box off

        for jj = 2:3
            [H,P,CI] = ttest(bal_smoo_y(:,jj),0);
            if P<=.05
                text(jj,ax1.YLim(2).*0.95,'*','FontSize',15,'HorizontalAl','center')
            else
                text(jj,ax1.YLim(2).*1,'n.s.','FontSize',8,'HorizontalAl','center')            
            end
        end

    ax3 = axes('Units','pixels','Position',[xnow ynow-360 160 130]);  
        min_bin = squeeze(sols(:,1,2,:)); % all methods, binsize, balanced, all measures
        min_smoo = squeeze(sols(:,2,2,:)); % all methods, smoothing, balanced, all measures

        min_bin_y = [];

        cols = winter(length(mnames));
        for jj = 1:6
            y = (min_bin(jj,:) - min_bin(jj,1))./(min_bin(jj,:) + min_bin(jj,1));
            min_bin_y = [min_bin_y;y];                        
            plot((2:3)+normrnd(0,jit,[1 2]),y(2:3),'Color',cols(jj,:),'Marker','o'); hold on;
        end

        ax3.YLim = ax1.YLim;        
        % ax3.YScale = 'log';
        ax3.XTick = [];
        ylabel(sprintf('Minimum bin size (mm)\n(a-b)/(a+b)'))
        ax3.XLim = [1.5 3.5];
        text(0,1.15,'Minimum solution bin size','Units','normalized','FontSize',10)
        line(ax1.XLim,[0 0],'Color',[.5 .5 .5])
        box off

        for jj = 2:3
            [H,P,CI] = ttest(min_bin_y(:,jj),0);
            if P<=.05
                text(jj,ax1.YLim(2).*0.95,'*','FontSize',15,'HorizontalAl','center')
            else
                text(jj,ax1.YLim(2).*1,'n.s.','FontSize',8,'HorizontalAl','center')            
            end
        end

    ax4 = axes('Units','pixels','Position',[xnow ynow-540 160 130]);  
        min_smoo_y = [];
    
        ps = [];
        for jj = 1:6
            y = (min_smoo(jj,:) - min_smoo(jj,1))./(min_smoo(jj,:) + min_smoo(jj,1));
            min_smoo_y = [min_smoo_y;y];                                    
            ps(jj) = plot((2:3)+normrnd(0,jit,[1 2]),y(2:3),'Color',cols(jj,:),'Marker','o'); hold on;
        end

        ax4.YLim = ax1.YLim;        
        % ax4.YScale = 'log';
        ax4.XTick = 1:3;
        ax4.XTickLabel = enames;
        ylabel(sprintf('Minimum smoothing (mm)\n(a-b)/(a+b)'))
        ax4.XLim = [1.5 3.5];
        text(0,1.15,'Minimum solution smoothing','Units','normalized','FontSize',10)
        line(ax1.XLim,[0 0],'Color',[.5 .5 .5])
        box off

        for jj = 2:3
            [H,P,CI] = ttest(min_smoo_y(:,jj),0);
            if P<=.05
                text(jj,ax1.YLim(2).*0.95,'*','FontSize',15,'HorizontalAl','center')
            else
                text(jj,ax1.YLim(2).*1,'n.s.','FontSize',8,'HorizontalAl','center')            
            end
        end

        [~,leg] = legendflex(ps,dnames,'anchor',{'s','s'},'ncol',2,'box','off','buffer',[0,-120],'xscale',0.5,'fontsize',9); 

        % return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ################################################################# %% Save the figure
    % Save the figure    
    disp(sprintf('\tSaving figure...'))    
    figname = [fig_dir2 '\fig_error_measures.png'];
    [~,~,~] = mkdir(fig_dir2);    
    exportgraphics(fig1,figname,'BackgroundColor','w','ContentType','image','Resolution',250);  
    close(fig1)    
    
    


  
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    