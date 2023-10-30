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
ovals = {'overdisperse0-0','overdisperse0-8'};
nmes = {'No overdispersion','Overdispersion'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ################################################################# %% FUNCTION BODY
%% #################### Heatmap of error by binsize and smoothing
    % Figure settings
    fig1 = figure('visible','on','Position',[50,60,1000,900]); 
    set(gcf,'InvertHardCopy','off'); % this stops white lines being plotted black      
    set(gcf,'color','w'); % makes the background colour white
    colormap(jet(256)); % to make sure the colormap is not the horrible default one

    % method settings
    xnow = 520;
    ynow = 730;
    mnames = {'histogram','ash','kadaptive','ksde','fyhn','kyadaptive'};
    
    dnames = {'Histogram','ASH','Adaptive smoothing','KSDE','tKSDE','Adaptive binning'}; 
    ls = {'a','b','c','d','e','f','g'};
    mapidx = 4;
    xsiz = 100;
    xbuff = 20;
    xvec = xnow : (xsiz+xbuff) : 1000;
    ysiz = 100;
    ybuff = 30;
    yvec = ynow : -(ysiz+ybuff) : 0;
    [xx,yy] = meshgrid(xvec,yvec);
    xx = xx';
    yy = yy;

xx = [xnow xnow+xsiz+xbuff];
yvec = ynow : -(ysiz+ybuff) : 0;


    % plot settings
    n_color_levels = 64;
    error_colormap = flipud(inferno(n_color_levels));    
    clims = [1*10^-12 1*10^-11];   
    interp_method = 'nearest';  

    sols = NaN(length(mnames),2,2,length(ovals));
  
    for mm = 1:length(mnames) % for every method
        disp(sprintf('\t%s...',mnames{mm}))
        
        for oo = 1:length(ovals) % for every overdispersion value
            mname = mnames{mm};
            sigma_now = '16000';
            time_now = '4'; 
            disp(sprintf('\t\t...%s',ovals{oo}))
    
            % main panel showing 16 minute error results
            ax1 = axes('Units','pixels','Position',[xx(oo) yvec(mm) xsiz ysiz]);  
                disp(sprintf('\t\t...errormap'))

                load([scratch_space '\' ename '_' mname '_sigma' num2str(sigma_now) '_duration' num2str(time_now) ovals{oo} '_datamats.mat'],'data_matrix_error','data_matrix_params')  
                bmat = data_matrix_params(:,:,1); % binsizes used when generating maps
                smat = data_matrix_params(:,:,2); % smoothing used when generating maps 
    
                % switch mname
                %     case {'histogram'}
                %         ah = add_panel_title(ls{mm},sprintf(dnames{mm}),'yoffset',5,'xoffset',0,'width',400);
                %     otherwise
                if mm==1 && oo==1
                    ah = add_panel_title('',sprintf(dnames{mm}),'yoffset',-15,'xoffset',0,'width',400);
                elseif mm>1 && oo==1
                    ah = add_panel_title('',sprintf(dnames{mm}),'yoffset',-15,'xoffset',0,'width',400);
                end
                % end
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
                % switch mname
                %     case {'histogram'}
                %         text(0,1.15,'4 minutes coverage','HorizontalAl','left','Units','normalized','FontSize',8);                
                %     otherwise
                if mm==1 && oo==1
                    % text(0,1.04,'4 minutes coverage','HorizontalAl','left','Units','normalized','FontSize',8);
                    text(-0.15,1.35,nmes{1},'HorizontalAl','left','Units','normalized','FontSize',11);
                    text(-0.40,1.38,'b','FontSize',18,'Units','normalized')                    
                elseif mm==1 && oo==2 
                    text(0,1.35,nmes{2},'HorizontalAl','left','Units','normalized','FontSize',11);
                end
                % end  
                if mm==length(mnames)
                    xlabel(grid_text);    
                else
                    ax1.XTick = [];
                end
                if oo==1
                    ylabel(smoo_text);                    
                else
                    ax1.YTick = [];
                end

%% #################### calculate Pareto front
            if 1
                disp(sprintf('\t\t...pareto front'))

                fname = [scratch_space '\' ename '_' mname '_sigma' num2str(sigma_now) '_duration' num2str(time_now) ovals{oo} '_gamultiobj.mat'];
                if exist(fname,'file') && 1 && ~override_pareto
                    % do nothing
                else
                    load([scratch_space '\' ename '_' mname '_sigma' num2str(sigma_now) '_duration' num2str(time_now) ovals{oo} '_datamats.mat'],'data_matrix_cputime','data_matrix_fields','data_matrix_error','data_matrix_missing','data_matrix_params')   
                
                    disp(sprintf('\t\t\t...calculating'))
                    
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
                    save([scratch_space '\' ename '_' mname '_sigma' num2str(sigma_now) '_duration' num2str(time_now) ovals{oo} '_gamultiobj.mat'],'x','fval','balance_solution1','balance_solution2','balance_solution3','minerr_solution');                
                end
                
            end        

%% #################### Plot Solutions onto error map 
            if 1     
                disp(sprintf('\t\t...plotting solutions'))
                
                load([scratch_space '\' ename '_' mname '_sigma' num2str(sigma_now) '_duration' num2str(time_now) ovals{oo} '_gamultiobj.mat'],'balance_solution2','balance_solution1','minerr_solution','x');
                % find the 'best' solution (closest to minimising everything)
                pb1 = plot3(balance_solution1(1),balance_solution1(2),1000,'ko','MarkerSize',8,'LineStyle','none');   
                sols(mm,:,1,oo) = [balance_solution1(1),balance_solution1(2)];

                % find the minimal error solution 
                pc = plot3(minerr_solution(1),minerr_solution(2),1000,'kv','MarkerSize',8,'LineStyle','none');          
                sols(mm,:,2,oo) = [minerr_solution(1),minerr_solution(2)];

                % [~,leg] = legendflex([pb1 pc]...
                %     ,{'Balanced solution','Minimum error solution'}...
                %     ,'anchor',{'se','se'},'ncol',1,'box','off','buffer',[-30,-90],'xscale',.5,'fontsize',9);    
            end

%% #################### colorbar
            if mm==1 && oo==2
                axc = axes('Units','pixels','Position',[ax1.Position(1)+ax1.Position(3)+10 ax1.Position(2) 12 ax1.Position(4)]); 
                    mat = linspace(ax1.CLim(1),ax1.CLim(2),100)';
                    imagesc(ones(size(mat)),mat,mat);
                    colormap(axc,error_colormap);
                    axis xy

                    axc.YTick = linspace(ax1.CLim(1),ax1.CLim(2),4);
                    axc.XTick = [];
                    axc.YAxisLocation = 'right';
                    text(0.5,1.2,sprintf('MISE'),'FontSize',8,'HorizontalAl','center','Units','normalized')
            end  
        end
    end

%% #################### Example overdispersion
    xnow = 70;
    ynow = 680;
    ax1 = axes('Units','pixels','Position',[xnow ynow 160 160]);  
        ah = add_panel_title('a',sprintf('Bin size, smoothing & overdispersion'),'yoffset',10,'xoffset',-25,'width',300);                    

            sdir = [scratch_space '\' ename '_sigma16000_duration64' ovals{1} '_sdata.mat'];            
            disp(sprintf('\t\t...loading %s',sdir));            
            load(sdir,'sdata'); 

            % cut the position and spike data to the trial length
            % we want
            excell = 9;
            excell = 17;
            pindx = sdata.pos_index(excell);
            pos_now = all_pos{pindx,1};
            pox = pos_now(:,1);
            poy = pos_now(:,2);
            pot = pos_now(:,3); 

            spk = sdata.spk{excell};
            spx = spk(:,1);
            spy = spk(:,2);
            spt = spk(:,3);

            tcut = 4*60;
            ppox = pox(pot<tcut);
            ppoy = poy(pot<tcut);
            ppot = pot(pot<tcut);            
            pspx = spx(spt<tcut);
            pspy = spy(spt<tcut);
            pspt = spt(spt<tcut);            
            pos = [ppox ppoy]; % positions in mm
            spk = [pspx pspy]; % spikes in mm

            plot(ppox,ppoy,'k'); hold on;
            plot(pspx,pspy,'r.','MarkerSize',20);
            daspect([1 1 1])
            axis xy off
            text(0,1.1,sprintf(nmes{1}),'Units','normalized','HorizontalAlignment','left','FontSize',10)

        ax1b = axes('Units','pixels','Position',[ax1.Position(1) ax1.Position(2)-90 120 60]); % inset
            tcut = 64*60;
            ppox = pox(pot<tcut);
            ppoy = poy(pot<tcut);
            ppot = pot(pot<tcut);            
            pspx = spx(spt<tcut);
            pspy = spy(spt<tcut);
            pspt = spt(spt<tcut);            
            pos = [ppox ppoy]; % positions in mm
            spk = [pspx pspy]; % spikes in mm

            rmset = mapset; % rate mapper settings structure - passing mapset directly to ratemapper causes a memory leak
            rmset.pot = ppot;
            rmset.spt = pspt;
            [ratemap,dwellmap,spikemap,rmset,speedlift] = rate_mapper(pos,spk,rmset);

            [z,overdispersion,r] = getOVERDISPERSE(rmset.map_pos,ppot,pspt,ratemap);

            xi = -10:0.5:10;
            % f = ksdensity(z(:),xi,'bandwidth',1,'kernel','normal');
            % area(xi,f,'FaceColor','k','FaceAlpha',0.2,'EdgeColor','none');   
            f = histcounts(z(:),xi,'Normalization','pdf');
            x = movmean(xi,2,'EndPoints','discard');
            bar(x,f,1,'k');            
            hold on;
            box off

            xi = -10:0.1:10;            
            y = normpdf(xi,0,1);
            plot(xi,y,'r'); 

            xlabel('Z score')
            ylabel('Probability density')
            ax1b.XLim = [-8 8];
            text(0.1,1.15,sprintf('%c = %.2f',963,overdispersion),'Units','normalized','HorizontalAlignment','left','FontSize',10)

    ax1 = axes('Units','pixels','Position',[ax1b.Position(1)+210 ynow 160 160]);  
            sdir = [scratch_space '\' ename '_sigma16000_duration64' ovals{2} '_sdata.mat'];            
            disp(sprintf('\t\t...loading %s',sdir));            
            load(sdir,'sdata'); 

            % cut the position and spike data to the trial length
            % we want
            excell = 13;
            excell = 20;
            pindx = sdata.pos_index(excell);
            pos_now = all_pos{pindx,1};
            pox = pos_now(:,1);
            poy = pos_now(:,2);
            pot = pos_now(:,3); 

            spk = sdata.spk{excell};
            spx = spk(:,1);
            spy = spk(:,2);
            spt = spk(:,3);

            tcut = 4*60;
            ppox = pox(pot<tcut);
            ppoy = poy(pot<tcut);
            ppot = pot(pot<tcut);            
            pspx = spx(spt<tcut);
            pspy = spy(spt<tcut);
            pspt = spt(spt<tcut);            
            pos = [ppox ppoy]; % positions in mm
            spk = [pspx pspy]; % spikes in mm

            plot(ppox,ppoy,'k'); hold on;
            plot(pspx,pspy,'r.','MarkerSize',20);
            daspect([1 1 1])
            axis xy off
            text(0,1.1,sprintf(nmes{2}),'Units','normalized','HorizontalAlignment','left','FontSize',10)

        ax1b = axes('Units','pixels','Position',[ax1.Position(1) ax1.Position(2)-90 120 60]); % inset
            tcut = 64*60;
            ppox = pox(pot<tcut);
            ppoy = poy(pot<tcut);
            ppot = pot(pot<tcut);            
            pspx = spx(spt<tcut);
            pspy = spy(spt<tcut);
            pspt = spt(spt<tcut);            
            pos = [ppox ppoy]; % positions in mm
            spk = [pspx pspy]; % spikes in mm

            rmset = mapset; % rate mapper settings structure - passing mapset directly to ratemapper causes a memory leak
            rmset.pot = ppot;
            rmset.spt = pspt;
            [ratemap,dwellmap,spikemap,rmset,speedlift] = rate_mapper(pos,spk,rmset);

            [z,overdispersion,r] = getOVERDISPERSE(rmset.map_pos,ppot,pspt,ratemap);

            xi = -10:0.5:10;
            % f = ksdensity(z(:),xi,'bandwidth',1,'kernel','normal');
            % area(xi,f,'FaceColor','k','FaceAlpha',0.2,'EdgeColor','none');  
            f = histcounts(z(:),xi,'Normalization','pdf');
            x = movmean(xi,2,'EndPoints','discard');
            bar(x,f,1,'k');            
            hold on;
            box off

            xi = -10:0.1:10;            
            y = normpdf(xi,0,1);
            plot(xi,y,'r'); 

            xlabel('Z score')
            ylabel('Probability density')
            ax1b.XLim = [-8 8];
            text(0.1,1.1,sprintf('%c = %.2f',963,overdispersion),'Units','normalized','HorizontalAlignment','left','FontSize',10)
% return
%% #################### Differences in binsize and smoothing
    xnow = 80;
    ynow = ynow-310;
    ax1 = axes('Units','pixels','Position',[xnow ynow 140 130]);  
        ah = add_panel_title('c',sprintf('Bin size, smoothing & overdispersion'),'yoffset',5,'xoffset',-30,'width',300);                    

        b1 = sols(:,1,1,1); % binsize, balanced, no overdispersion
        b2 = sols(:,1,1,2); % binsize, balanced, overdispersion
        s1 = sols(:,2,1,1); % smoothing, balanced, no overdispersion
        s2 = sols(:,2,1,2); % smoothing, balanced, overdispersion

        b3 = sols(:,1,2,1); % binsize, minimised, no overdispersion
        b4 = sols(:,1,2,2); % binsize, minimised, overdispersion
        s3 = sols(:,2,2,1); % smoothing, minimised, no overdispersion
        s4 = sols(:,2,2,2); % smoothing, minimised, overdispersion

        y1 = (b2-b1) ./ (b2+b1);
        y2 = (b4-b3) ./ (b4+b3);                
        cols = winter(length(mnames));
        msiz = 50;
        malph = 1;
        jit = 0.05;
        errorbar(zeros(size(y1))+1,mean(y1),std(y1),'ko'); hold on;   
        errorbar(zeros(size(y2))+2,mean(y2),std(y2),'k^');
        m1 = scatter(zeros(size(y1))+1.5+normrnd(0,jit,size(y2)),y1,msiz,cols,'filled','o','MarkerFaceAlpha',malph); hold on;          
        m2 = scatter(zeros(size(y2))+2.5+normrnd(0,jit,size(y2)),y2,msiz,cols,'filled','^','MarkerFaceAlpha',malph);

        % ax1.XColor = 'none';
        ax1.XLim = [0.5 3];
        ylabel(sprintf('Bin size difference\n(a-b)/(a+b)'))
        box off
        line(ax1.XLim,[0 0],'Color',[.5 .5 .5]);
        % text(1.25,ax1.YLim(1),'No overisp.','FontSize',10,'HorizontalAlignment','right','rotation',25)
        % text(2.25,ax1.YLim(1),'Overdisp.','FontSize',10,'HorizontalAlignment','right','rotation',25)
        ax1.XTick = [1.25 2.25];
        ax1.XTickLabel = {'No overisp.','Overdisp.'};
        xtickangle(ax1,25);

        p1 = plot(-1000,1,'ko','MarkerFaceColor','k');
        p2 = plot(-1000,1,'k^','MarkerFaceColor','k');
        [~,leg] = legendflex([p1 p2]...
            ,{'Balanced solution','Minimum error solution'}...
            ,'anchor',{'s','s'},'ncol',1,'box','off','buffer',[-20,-80],'xscale',0.5,'fontsize',9); 

        xm = [];
        for cm = 1:length(mnames)
            xm(cm) = plot(-1000,1,'color',cols(cm,:),'marker','o','LineStyle','none','MarkerFaceColor',cols(cm,:));
        end
        [~,leg] = legendflex(xm...
            ,dnames...
            ,'anchor',{'s','s'},'ncol',2,'box','off','buffer',[190,-100],'xscale',0.5,'fontsize',9); 

        % stats
        [H,P,CI] = ttest(y1,0);
        if P<=.05
            text(1.25,ax1.YLim(2).*1.15,'*','FontSize',20,'HorizontalAl','center','VerticalAl','top')
        else
            text(1.25,ax1.YLim(2).*1.6,'n.s.','FontSize',8,'HorizontalAl','center','VerticalAl','top')            
        end
        [H,P,CI] = ttest(y2,0);
        if P<=.05
            text(2.25,ax1.YLim(2).*1.15,'*','FontSize',20,'HorizontalAl','center','VerticalAl','top')
        else
            text(2.25,ax1.YLim(2).*1.6,'n.s.','FontSize',8,'HorizontalAl','center','VerticalAl','top')                       
        end        

    ax2 = axes('Units','pixels','Position',[xnow+ax1.Position(3)+90 ynow ax1.Position(3) ax1.Position(4)]);  
        y1 = (s2-s1) ./ (s2+s1);
        y2 = (s4-s3) ./ (s4+s3);                
        errorbar(zeros(size(y1))+1,mean(y1),std(y1),'ko'); hold on;   
        errorbar(zeros(size(y2))+2,mean(y2),std(y2),'k^');
        m1 = scatter(zeros(size(y1))+1.5+normrnd(0,jit,size(y2)),y1,msiz,cols,'filled','o','MarkerFaceAlpha',malph); hold on;          
        m2 = scatter(zeros(size(y2))+2.5+normrnd(0,jit,size(y2)),y2,msiz,cols,'filled','^','MarkerFaceAlpha',malph);

        % ax2.XColor = 'none';
        ax2.XLim = [0.5 3];        
        ylabel(sprintf('Smoothing difference\n(a-b)/(a+b)'))
        box off
        line(ax2.XLim,[0 0],'Color',[.5 .5 .5]);
        % text(1.25,ax2.YLim(1)-(ax2.YLim(1).*0.1),'No overisp.','FontSize',10,'HorizontalAlignment','right','rotation',25)
        % text(2.25,ax2.YLim(1)-(ax2.YLim(1).*0.1),'Overdisp.','FontSize',10,'HorizontalAlignment','right','rotation',25)
        ax2.XTick = [1.25 2.25];
        ax2.XTickLabel = {'No overisp.','Overdisp.'};
        xtickangle(ax2,25);

        % stats
        [H,P,CI] = ttest(y1,0);
        if P<=.05
            text(1.25,ax2.YLim(2).*1.15,'*','FontSize',20,'HorizontalAl','center','VerticalAl','top')
        else
            text(1.25,ax2.YLim(2).*1.6,'n.s.','FontSize',8,'HorizontalAl','center','VerticalAl','top')
        end
        [H,P,CI] = ttest(y2,0);
        if P<=.05
            text(2.25,ax2.YLim(2).*1.15,'*','FontSize',20,'HorizontalAl','center','VerticalAl','top')
        else
            text(2.25,ax2.YLim(2).*1.6,'n.s.','FontSize',8,'HorizontalAl','center','VerticalAl','top')
        end

        % return
%% #################### error and place field error
    mapidx = 4;

    sigma_now = '16000';
    time_now = '4';   
    err_vals = NaN(config.npcells,length(mnames),2,length(ovals));
    p_vals = NaN(config.npcells,length(mnames),2,length(ovals));  
    m_vals = NaN(config.npcells,length(mnames),2,length(ovals));                
    time_vals = NaN(config.npcells,length(mnames),2,length(ovals));  
    bsols = NaN(length(mnames),2,length(ovals));
    msols = NaN(length(mnames),2,length(ovals));        
    for mm = 1:length(mnames) % for every mapping method
        mname = mnames{mm};
        for oo = 1:length(ovals) % for every overdispersion value

            % get error matrix
            load([scratch_space '\' ename '_' mname '_sigma' num2str(sigma_now) '_duration' num2str(time_now) ovals{oo} '_datamats.mat'],'data_matrix_error','data_matrix_params','data_matrix_fields','data_matrix_cputime','data_matrix_missing')  
            bmat = data_matrix_params(:,:,1); % binsizes used when generating maps
            smat = data_matrix_params(:,:,2); % smoothing used when generating maps
            dme = data_matrix_error(:,:,:,mapidx);
            dmt = data_matrix_cputime;
            dmf = data_matrix_fields;
            dmm = data_matrix_missing;

            % get balanced solution values
            fname = [scratch_space '\' ename '_' mname '_sigma' num2str(sigma_now) '_duration' num2str(time_now) ovals{oo} '_gamultiobj.mat'];
            load(fname,'x','fval','balance_solution1','balance_solution2','balance_solution3','minerr_solution');
            bsols(mm,:,oo) = balance_solution1(1:2);
            msols(mm,:,oo) = minerr_solution(1:2);

            % get all error values for this solution
            for ii = 1:size(dme,3)
                err_vals(ii,mm,1,oo) = interp2(bmat,smat,dme(:,:,ii),balance_solution1(1),balance_solution1(2));
                err_vals(ii,mm,2,oo) = interp2(bmat,smat,dme(:,:,ii),minerr_solution(1),minerr_solution(2)); 

                p_vals(ii,mm,1,oo) = interp2(bmat,smat,dmf(:,:,ii),balance_solution1(1),balance_solution1(2),'nearest');
                p_vals(ii,mm,2,oo) = interp2(bmat,smat,dmf(:,:,ii),minerr_solution(1),minerr_solution(2),'nearest'); 

                m_vals(ii,mm,1,oo) = interp2(bmat,smat,dmm(:,:,ii),balance_solution1(1),balance_solution1(2));
                m_vals(ii,mm,2,oo) = interp2(bmat,smat,dmm(:,:,ii),minerr_solution(1),minerr_solution(2));                 
                if strcmp(mname,'fyhn')
                    if ii>10
                        continue
                    end
                end
                time_vals(ii,mm,1,oo) = interp2(bmat,smat,dmt(:,:,ii),balance_solution1(1),balance_solution1(2));            
                time_vals(ii,mm,2,oo) = interp2(bmat,smat,dmt(:,:,ii),minerr_solution(1),minerr_solution(2));            
            end
        end
    end

    ynow = ynow-270;

    % error results
    ax3 = axes('Units','pixels','Position',[ax1.Position(1) ynow ax1.Position(3) ax1.Position(4)]); 
        ah = add_panel_title('d',sprintf('MISE, field detection & dispersion'),'yoffset',5,'xoffset',-30,'width',300);                    
   
        % plot results
        dispersion = [1 2];
        cols = winter(length(mnames));
        means = NaN(length(mnames),length(ovals)); % means
        sems = NaN(length(mnames),length(ovals)); % SEMs
        for mm = 1:length(mnames) % for every mapping method
            enow = squeeze(err_vals(:,mm,1,:)); % balanced
            means(mm,:) = mean(enow,1,'omitnan');
            sems(mm,:) = nansem(enow,1);

            errorbar(dispersion,means(mm,:),sems(mm,:),'Color',cols(mm,:),'Marker','o'); hold on;
        end

        ax = gca;
        ax.XLim = [0.5 2.5];
        ax.XTick = dispersion;
        ax.XTickLabel = nmes;        
        ax.YTick = (ones(1,20).*10).^-[20:-1:1];        
        % ax.YLim = [-1*10^-12 16*10^-12]; 
        % ax.YTick = 0:1:5;
        % ax.XTickLabel = [];
        ax.YScale = 'log'; 
        ax.YLim(1) = 10^-12;
        ax.YLim(2) = 10^-11;        
        ylabel(sprintf('MISE')) 
        % xlabel('Dispersion')
        % text(0.99,1.1,'N = 256 cells','Units','normalized','HorizontalAlignment','right','VerticalAlignment','top')
        box off

        [result1,result2] = plotsigbrackets([means(:,1); means(:,2)],[ones(size(means,1),1); ones(size(means,1),1).*2],'plot_omnibus',0);

    % field detection results
    ax4 = axes('Units','pixels','Position',[ax2.Position(1) ynow ax3.Position(3) ax3.Position(4)]); 
        % ah = add_panel_title('c',sprintf('Place field accuracy across mapping methods'),'yoffset',5,'xoffset',-30,'width',300);                    
   
        % plot results
        cols = winter(length(mnames));
        vals = NaN(length(mnames),length(ovals)); % 
        for mm = 1:length(mnames) % for every mapping method
            enow = squeeze(p_vals(:,mm,1,:)); % balanced
            vals(mm,:) = sum(enow==0,1,'omitnan') ./ size(enow,1);

            plot(dispersion+normrnd(0,0.05,size(vals(mm,:))),vals(mm,:),'Color',cols(mm,:),'Marker','o'); hold on;
        end

        ax = gca;
        ax.XLim = [0.5 2.5];
        ax.XTick = dispersion;
        ax.XTickLabel = nmes;
        ax.YLim = [0.94 1]; 
        ylabel(sprintf('Prop. correct place field count'))  
        % xlabel('Dispersion (a_{mod})')        
        % text(0.99,1.1,'N = 256 cells','Units','normalized','HorizontalAlignment','right','VerticalAlignment','top')
        box off

        [result1,result2] = plotsigbrackets([vals(:,1); vals(:,2)],[ones(size(vals,1),1); ones(size(vals,1),1).*2],'plot_omnibus',0);
        % text(1.5,ax.YLim(2).*1.1,'n.s.','FontSize',10,'HorizontalAl','center')            
        % 
        % keyboard
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ################################################################# %% Save the figure
    % Save the figure    
    disp(sprintf('\tSaving figure...'))    
    figname = [fig_dir2 '\v2_fig_overdispersion.png'];
    [~,~,~] = mkdir(fig_dir2);    
    exportgraphics(fig1,figname,'BackgroundColor','w','ContentType','image','Resolution',250);  
    close(fig1)    
    
    


  
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    