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
    if 1
        % Figure settings
        fig1 = figure('visible','on','Position',[50,50,1300,1050]); 
        set(gcf,'InvertHardCopy','off'); % this stops white lines being plotted black      
        set(gcf,'color','w'); % makes the background colour white
        colormap(jet(256)); % to make sure the colormap is not the horrible default one


%% #################### Literature values over time
        xnow = 100;
        ynow = 650;

        disp(sprintf('\tOver time...'))

        % main panel showing 16 minute error results
        ax1 = axes('Units','pixels','Position',[xnow+50 ynow+100 550 200]);    
            ah = add_panel_title('a',sprintf('Reported bin sizes over time'),'yoffset',10,'xoffset',-50,'width',400);
            zidx = dat(:,3)>0;        
            msiz = 80;

            rng(999);            
            jitter = 0.2;
            alph = 0.7;
            x = dat(:,1) + normrnd(0,jitter*10,size(dat(:,1)));
            y = dat(:,2) + normrnd(0,jitter,size(dat(:,1)));
            s2 = scatter(x(~zidx),y(~zidx),msiz,rgb('DarkSlateGray'),'filled','Marker','o','MarkerFaceAlpha',alph); hold on;
            s1 = scatter(x(zidx),y(zidx),msiz*1.2,rgb('Orange'),'filled','Marker','s','MarkerFaceAlpha',alph); hold on;

            xb = dat_box(:,1) + normrnd(0,jitter*10,size(dat_box(:,1)));
            yb = dat_box(:,2) + normrnd(0,jitter,size(dat_box(:,1)));
            s5 = scatter(xb,yb,50,rgb('MediumBlue'),'filled','Marker','^','MarkerFaceAlpha',alph); hold on;

            xlabel('Year of publication')
            ylabel('Bin size (mm)')
            s3 = line(ax1.XLim,[balance_solution1(1) balance_solution1(1)],'Color',[.5 .5 .5],'LineStyle','--','LineWidth',1.5);
            p = polyfit([dat(:,1);dat_box(:,1)], [dat(:,2);dat_box(:,2)], 1);
            [r,pr] = corr([dat(:,1);dat_box(:,1)], [dat(:,2);dat_box(:,2)],'Type','Spearman','rows','pairwise');
            s4 = refline(p(1),p(2));
            set(s4,'Color','k','LineStyle','-','LineWidth',1.5)

            ax1.XLim = [1985 2024];

            lbl = {'No smoothing','Boxcar smoothing','Gaussian smoothing','Pareto-optimal bin size',sprintf('Linear fit (r = %.1f,p=%.3f)',r,pr)};
            [~,leg] = legendflex([s2 s5 s1 s3 s4],lbl,'anchor',{'ne','ne'},'ncol',1,'box','off','buffer',[50,40],'xscale',.5,'fontsize',9);  

    % figure
    % xi = 1985:0.5:2025;
    % bw = 1.5;
    % y1 = ksdensity(dat(~zidx,1),xi,'Kernel','normal','Bandwidth',bw);
    % y2 = ksdensity(dat(zidx,1),xi,'Kernel','normal','Bandwidth',bw);
    % y3 = ksdensity(dat_box(:,1),xi,'Kernel','normal','Bandwidth',bw);
    %      
    % plot(xi,y1,'b'); hold on;
    % plot(xi,y2,'g');
    % plot(xi,y3,'r');

%% #################### Example maps
        ynow = ynow-130;
        xnow = xnow;
        if 1
            disp(sprintf('\tExample rate maps...'))        
            sdir = [scratch_space '\' ename '_sigma8000_duration64_sdata.mat'];            
            disp(sprintf('\t\t...loading %s',sdir));            
            load(sdir,'sdata'); 

            psiz = 140;
            hbuff = 10;
            xvec = [xnow xnow+psiz+hbuff xnow+(psiz*2)+(hbuff*2) xnow+(psiz*3)+(hbuff*3) xnow+(psiz*4)+(hbuff*4) xnow+(psiz*5)+(hbuff*5) xnow+(psiz*6)+(hbuff*6) xnow+(psiz*7)+(hbuff*7)]+10;

            excell = [1 1 1];
            maps_for_plotting = cell(1,length(excell));
            bvec = [39 30 26];
            svec = [0 3 35];
            smeth = [0 4 1];
            tvec = {sprintf('Unsmothed\nbin size = 39 mm'),sprintf('3 \\times 3 boxcar\nbin size = 30 mm'),sprintf('Gaussian smoothed\nbin size = 26 mm\n%c = 35 mm',963)};
            tvecy = {'1985-1995','2000-2010','2015-2023'};

            override_maps = 0;
            fname = [fig_dir mname '_over_time_maps.mat'];
            if exist(fname,'file') && ~override_maps
                mem = 1;
                disp(sprintf('\t\t...loading old maps'))
                load(fname,'maps_for_plotting')
            else
                disp(sprintf('\t\t...generating maps'))            
                mem = 0;
            end

            pindx = sdata.pos_index(excell);
            pos_now = all_pos{pindx,1};
            pox = pos_now(:,1);
            poy = pos_now(:,2);
            pot = pos_now(:,3); 
        
            spk = sdata.spk{excell};
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
        
            ax = axes('Units','pixels','Position',[xvec(1) ynow psiz psiz]);    
                plot(ppox,ppoy,'Color',[.5 .5 .5 .5]); hold on;
                scatter(pspx,pspy,20,'r','filled','MarkerEdgeColor','none','MarkerFaceAlpha',0.5);
                daspect([1 1 1])
                axis xy off  

            loopout = looper(length(excell));     
            for xx = 1:3
                ax = axes('Units','pixels','Position',[xvec(xx+1) ynow psiz psiz]);    
                    if mem
                        ratemap = maps_for_plotting{1,xx};
                    else                    
                        rmset = mapset; % rate mapper settings structure - passing mapset directly to ratemapper causes a memory leak
                        rmset.method = 'histogram';
                        rmset.binsize = bvec(xx);
                        rmset.ssigma = svec(xx);
                        rmset.smethod = smeth(xx);
                        rmset.pot = ppot;
                        rmset.spt = pspt;
                        rmset.maxdist = 320;
                        rmset.mindist = 50;
                        rmset.maplims = [-max(abs(epoly(:,1))) -max(abs(epoly(:,2))) max(abs(epoly(:,1))) max(abs(epoly(:,2)))];                     
                        [ratemap,dwellmap,spikemap,rmset,speedlift] = rate_mapper(pos,spk,rmset);
                        maps_for_plotting{1,xx} = ratemap;
                    end

                    im = imagesc(ratemap,'alphadata',~isnan(ratemap)); hold on;
                    daspect([1 1 1])
                    axis xy off tight
                    colormap(gca,rmap_colormap)    
                    loopout = looper(loopout);  
                    text(0.05,-0.05,tvec{xx},'Units','normalized','VerticalAl','top')
                    text(0.5,1.01,tvecy{xx},'Units','normalized','VerticalAl','bottom','HorizontalAl','center')
            end

            axc = axes('Units','pixels','Position',[ax.Position(1)+ax.Position(3)+25 ax.Position(2) 12 ax.Position(4)]); 
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
        g1 = getframe(fig1);
        % close(fig1);
    end
    

%% #################### Literature values on MISE heatmap
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

    ynow = 150;
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

% return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ################################################################# %% Save the figure
    % Save the figure    
    disp(sprintf('\tSaving figure...'))    
    figname = [fig_dir2 '\v2_fig_lit_review.png'];
    [~,~,~] = mkdir(fig_dir2);    
    exportgraphics(fig1,figname,'BackgroundColor','w','ContentType','image','Resolution',250);  
    close(fig1)    
   
    
    
    
    
    
    
    






