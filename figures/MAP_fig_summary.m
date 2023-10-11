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

    if ~exist('sdata','var')
        sdir = [scratch_space '\' config.environs{1} '_sigma16000_duration64_sdata.mat'];            
        disp(sprintf('\t\t...loading %s',sdir));            
        load(sdir,'sdata'); 
    end

%% ##################################### Heading 2
%% #################### Heading 3
%% Heading 4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ################################################################# %% FUNCTION BODY
%% #################### Heatmap of error by binsize and smoothing
    % Figure settings
    fig1 = figure('visible','on','Position',[50,60,1550,1050]); 
    set(gcf,'InvertHardCopy','off'); % this stops white lines being plotted black      
    set(gcf,'color','w'); % makes the background colour white
    colormap(jet(256)); % to make sure the colormap is not the horrible default one

    xnow = 80;
    ynow = 800;

    mapidx = 4;
    disp(sprintf('\tMISE map...'))

    % main panel showing 16 minute error results
    ax1 = axes('Units','pixels','Position',[xnow ynow-5 300 170]); 
        ah = add_panel_title('a',sprintf('Error across mapping methods'),'yoffset',5,'xoffset',-30,'width',300);                    
   
        mnames = {'fyhn','histogram','ash','kadaptive','ksde','kyadaptive'};
        dnames = {'tKSDE','Histogram','ASH','Adaptive\nsmoothing','KSDE','Adaptive\nbinning'};
        sigma_now = '4000';
        time_now = '16';   
        err_vals = NaN(256,length(mnames),2);
        p_vals = NaN(256,length(mnames),2);  
        m_vals = NaN(256,length(mnames),2);                
        time_vals = NaN(256,length(mnames),2);  
        bsols = NaN(length(mnames),2);
        msols = NaN(length(mnames),2);        
        for mm = 1:length(mnames) % for every mapping method
            mname = mnames{mm};
     
            % get error matrix
            load([scratch_space '\' ename '_' mname '_sigma' num2str(sigma_now) '_duration' num2str(time_now) '_datamats.mat'],'data_matrix_error','data_matrix_params','data_matrix_fields','data_matrix_cputime','data_matrix_missing')  
            bmat = data_matrix_params(:,:,1); % binsizes used when generating maps
            smat = data_matrix_params(:,:,2); % smoothing used when generating maps
            dme = data_matrix_error(:,:,:,mapidx);
            dmt = data_matrix_cputime;
            dmf = data_matrix_fields;
            dmm = data_matrix_missing;

            % get balanced solution values
            fs = 4000;
            ds = 16;
            fname = [scratch_space '\' ename '_' mname '_sigma' num2str(fs) '_duration' num2str(ds) '_gamultiobj.mat'];
            load(fname,'x','fval','balance_solution1','balance_solution2','balance_solution3','minerr_solution');
            bsols(mm,:) = balance_solution1(1:2);
            msols(mm,:) = minerr_solution(1:2);

            % get all error values for this solution
            for ii = 1:size(dme,3)
                err_vals(ii,mm,1) = interp2(bmat,smat,dme(:,:,ii),balance_solution1(1),balance_solution1(2));
                err_vals(ii,mm,2) = interp2(bmat,smat,dme(:,:,ii),minerr_solution(1),minerr_solution(2)); 

                p_vals(ii,mm,1) = interp2(bmat,smat,dmf(:,:,ii),balance_solution1(1),balance_solution1(2),'nearest');
                p_vals(ii,mm,2) = interp2(bmat,smat,dmf(:,:,ii),minerr_solution(1),minerr_solution(2),'nearest'); 

                m_vals(ii,mm,1) = interp2(bmat,smat,dmm(:,:,ii),balance_solution1(1),balance_solution1(2));
                m_vals(ii,mm,2) = interp2(bmat,smat,dmm(:,:,ii),minerr_solution(1),minerr_solution(2));                 
                if strcmp(mname,'fyhn')
                    if ii>10
                        continue
                    end
                end
                time_vals(ii,mm,1) = interp2(bmat,smat,dmt(:,:,ii),balance_solution1(1),balance_solution1(2));            
                time_vals(ii,mm,2) = interp2(bmat,smat,dmt(:,:,ii),minerr_solution(1),minerr_solution(2));            
            end
        end

        [g,~] = meshgrid(1:size(err_vals,2),1:size(err_vals,1));
        gs = g(:);
        ds = reshape(err_vals(:,:,1),[],1);
        cols = winter(length(mnames));
        cmat = mat2cell(cols,ones(size(cols,1),1));
        meanplot(ds,gs-0.2,'linecolor',{'none'},'meancolor',{'k'},'meanmarker',{'o'},'meansize',6,'dots',1,'dotsize',8,'dotalpha',0.5,'dotsigma',0.05,'dotcolor',cmat,'meanlinewidth',1.5); hold on;

        ds = reshape(err_vals(:,:,2),[],1);
        cols = winter(length(mnames));
        cmat = mat2cell(cols,ones(size(cols,1),1));
        meanplot(ds,gs+0.2,'linecolor',{'none'},'meancolor',{'k'},'meanmarker',{'v'},'meansize',6,'dots',1,'dotsize',8,'dotalpha',0.5,'dotsigma',0.05,'dotcolor',cmat,'meanlinewidth',1.5); hold on;

        x = [(mean(g,1,'omitnan')-0.2)' (mean(g,1,'omitnan')+0.2)'];
        y = [mean(err_vals(:,:,1),1,'omitnan')' mean(err_vals(:,:,2),1,'omitnan')'];
        plot(x',y','k');

        ax = gca;
        ax.XLim = [.5 length(mnames)+.5];
        ax.XTick = [1:length(mnames)];
        ax.YTick = (ones(1,20).*10).^-[20:-1:1];        
        % ax.YLim = [-1*10^-12 16*10^-12]; 
        % ax.YTick = 0:1:5;
        ax.XTickLabel = [];
        ax.YScale = 'log'; 
        ax.YLim(1) = 10^-13;
        ylabel(sprintf('MISE'))  
        text(0.99,1.1,'N = 256 cells','Units','normalized','HorizontalAlignment','right','VerticalAlignment','top')

        for ii = 1:length(mnames)
            text(ii,ax.YLim(1),sprintf(dnames{ii}),'VerticalAlignment','top','HorizontalAlignment','center','FontSize',8)
        end

        axis manual
        p1 = plot(-10,1,'ko','LineWidth',1.5);
        p2 = plot(-10,1,'kv','LineWidth',1.5);
        [~,leg] = legendflex([p1 p2],{'Balanced solution','Minimum error solution'},'anchor',{'nw','nw'},'ncol',1,'box','off','buffer',[0,20],'xscale',.5,'fontsize',9); 
        

%% #################### computation time
    xnow = xnow+400;
    sep = 10;
    ax1 = axes('Units','pixels','Position',[xnow ynow+sep 300 170-sep],'Color','none','Clipping','off'); 
        ah = add_panel_title('b',sprintf('Computation time across mapping methods'),'yoffset',0,'xoffset',-30,'width',300);                    
    
        [g,~] = meshgrid(1:size(time_vals,2),1:size(time_vals,1));
        gs = g(:);
        ds1 = reshape(time_vals(:,:,1),[],1);
        cols = winter(length(mnames));
        cmat = mat2cell(cols,ones(size(cols,1),1));
        meanplot(ds1,gs-0.2,'linecolor',{'none'},'meancolor',{'k'},'meanmarker',{'o'},'meansize',6,'dots',1,'dotsize',8,'dotalpha',0.5,'dotsigma',0.05,'dotcolor',cmat,'meanlinewidth',1.5); hold on;

        ds2 = reshape(time_vals(:,:,2),[],1);
        cols = winter(length(mnames));
        cmat = mat2cell(cols,ones(size(cols,1),1));
        meanplot(ds2,gs+0.2,'linecolor',{'none'},'meancolor',{'k'},'meanmarker',{'v'},'meansize',6,'dots',1,'dotsize',8,'dotalpha',0.5,'dotsigma',0.05,'dotcolor',cmat,'meanlinewidth',1.5); hold on;

        x = [(mean(g,1,'omitnan')-0.2)' (mean(g,1,'omitnan')+0.2)'];
        y = [mean(time_vals(:,:,1),1,'omitnan')' mean(time_vals(:,:,2),1,'omitnan')'];
        plot(x',y','k');

        ax = gca;
        ax.XLim = [.5 length(mnames)+.5];
        % ax.YLim(1) = 0;        
        ax.XTick = [];
        % ax.YTick = (ones(1,10).*10).^[1:10];
        ax.XTickLabel = [];
        ax.XColor = 'none';
        ax.YScale = 'log'; 
        ylabel(sprintf('Computation time (ms)'))  
        text(0.99,1.1,'N = 256 cells','Units','normalized','HorizontalAlignment','right','VerticalAlignment','top')

    ax2 = axes('Units','pixels','Position',ax1.Position-[0 sep 0 0],'Color','none'); 
        ax2.YColor = 'none';
        ax2.YTick = [];
        ax2.XLim = [.5 length(mnames)+.5];
        ax2.XTick = [];        

        % axis manual
        % p3 = plot(-10,1,'ko','LineWidth',1.5); hold on;
        % p4 = plot(-10,1,'kv','LineWidth',1.5);
        % [~,leg] = legendflex([p3 p4],{'Balanced solution','Minimum error solution'},'anchor',{'se','se'},'ncol',1,'box','off','buffer',[0,0],'xscale',.5,'fontsize',9); 

        for ii = 1:length(mnames)
            text(ii,ax.YLim(1),sprintf(dnames{ii}),'VerticalAlignment','top','HorizontalAlignment','center','FontSize',8)
        end


%% #################### field error
    xnow = 80;
    ynow = ynow-250;
    ax1 = axes('Units','pixels','Position',[xnow ynow 300 170],'Color','none','Clipping','off'); 
        ah = add_panel_title('c',sprintf('Place field accuracy across mapping methods'),'yoffset',0,'xoffset',-30,'width',300);                    

        x = [(mean(g,1,'omitnan')-0.2)' (mean(g,1,'omitnan')+0.2)'];
        y = [sum(p_vals(:,:,1)==0,1,'omitnan')' sum(p_vals(:,:,2)==0,1,'omitnan')'] ./ size(p_vals,1);
        cols = winter(length(mnames));
        cols = [cols;cols];
        for i = 1:length(x(:))
            if i>6
                m = 'v';
            else
                m = 'o';
            end
            stem(x(i),y(i),'filled','Color',cols(i,:),'Marker',m); hold on;
        end

        ax = gca;
        ax.XLim = [.5 length(mnames)+.5];
        ax.YLim = [0.95 1]; 
        ax.XTick = [];
        ylabel(sprintf('Prop. correct place field count'))  
        text(0.99,1.1,'N = 256 cells','Units','normalized','HorizontalAlignment','right','VerticalAlignment','top')
        box off

        for ii = 1:length(mnames)
            text(ii,ax.YLim(1),sprintf(dnames{ii}),'VerticalAlignment','top','HorizontalAlignment','center','FontSize',8)
        end

%% #################### empty bins
    xnow = xnow+400;
    ax1 = axes('Units','pixels','Position',[xnow ynow+sep 300 170-sep],'Color','none','Clipping','off'); 
        ah = add_panel_title('d',sprintf('Empty bins across mapping methods'),'yoffset',0,'xoffset',-30,'width',300);                    
    
        [g,~] = meshgrid(1:size(m_vals,2),1:size(m_vals,1));
        gs = g(:);
        ds1 = reshape(m_vals(:,:,1),[],1);
        cols = winter(length(mnames));
        cmat = mat2cell(cols,ones(size(cols,1),1));
        meanplot(ds1,gs-0.2,'linecolor',{'none'},'meancolor',{'k'},'meanmarker',{'o'},'meansize',6,'dots',1,'dotsize',8,'dotalpha',0.5,'dotsigma',0.05,'dotcolor',cmat,'meanlinewidth',1.5); hold on;

        ds2 = reshape(m_vals(:,:,2),[],1);
        cols = winter(length(mnames));
        cmat = mat2cell(cols,ones(size(cols,1),1));
        meanplot(ds2,gs+0.2,'linecolor',{'none'},'meancolor',{'k'},'meanmarker',{'v'},'meansize',6,'dots',1,'dotsize',8,'dotalpha',0.5,'dotsigma',0.05,'dotcolor',cmat,'meanlinewidth',1.5); hold on;

        x = [(mean(g,1,'omitnan')-0.2)' (mean(g,1,'omitnan')+0.2)'];
        y = [mean(m_vals(:,:,1),1,'omitnan')' mean(m_vals(:,:,2),1,'omitnan')'];
        plot(x',y','k');

        ax = gca;
        ax.XLim = [.5 length(mnames)+.5];
        % ax.YLim(1) = 0;        
        ax.XTick = [];
        % ax.YTick = (ones(1,10).*10).^[1:10];
        ax.XTickLabel = [];
        ax.XColor = 'none';
        % ax.YScale = 'log'; 
        ylabel(sprintf('Prop. empty bins'))  
        text(0.99,1.1,'N = 256 cells','Units','normalized','HorizontalAlignment','right','VerticalAlignment','top')

    ax2 = axes('Units','pixels','Position',ax1.Position-[0 sep 0 0],'Color','none'); 
        ax2.YColor = 'none';
        ax2.YTick = [];
        ax2.XLim = [.5 length(mnames)+.5];
        ax2.XTick = [];        

        for ii = 1:length(mnames)
            text(ii,ax.YLim(1),sprintf(dnames{ii}),'VerticalAlignment','top','HorizontalAlignment','center','FontSize',8)
        end

%% #################### Heatmap of error by binsize and smoothing
    dnames = {'tKSDE','Histogram','ASH','Adaptive smoothing','KSDE','Adaptive binning'};        
    excell = 140;
    excell = 20;

    xnow = 50;
    ynow = ynow-210;
    siz = 110; 
    spac = 10;
    xvec = xnow : siz+spac : 1000;
    yvec = ones(size(xvec))*ynow;

    ax1 = axes('Units','pixels','Position',[xvec(1) yvec(1) siz siz],'Color','none','Clipping','off');         
        pmap = sdata.pmap_1mm{excell,1};
        im = imagesc(pmap); hold on;
        daspect([1 1 1])
        axis xy off tight
        colormap(gca,rmap_colormap)

        text(0,1.05,sprintf('Spike probability map'),'HorizontalAl','left','Units','normalized','FontSize',8,'rotation',0,'Color','k','VerticalAl','bottom')  
        ah = add_panel_title('e',sprintf('Minimum error and balanced rate maps'),'yoffset',10,'xoffset',0,'width',300);    

    % cut the position and spike data to the trial length
    % we want
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

    ax1 = axes('Units','pixels','Position',[xvec(2) yvec(1) siz siz],'Color','none','Clipping','off');         
        plot(ppox,ppoy,'Color',[.5 .5 .5 .5]); hold on;
        scatter(pspx,pspy,5,'r','filled','MarkerEdgeColor','none','MarkerFaceAlpha',0.5);
        daspect([1 1 1])
        axis xy off  

        text(0,1.05,sprintf('Positions and spikes'),'HorizontalAl','left','Units','normalized','FontSize',8,'rotation',0,'Color','k','VerticalAl','bottom')  

    ynow = ynow-150;
    yvec = ones(size(xvec))*ynow;

    maps_for_plotting = cell(length(mnames),2);
    fname = [fig_dir '_summary_maps.mat'];
    override_maps = 0;
    if exist(fname,'file') && ~override_maps
        mem = 1;
        disp(sprintf('\t\t...loading old maps'))
        load(fname,'maps_for_plotting')
    else
        disp(sprintf('\t\t...generating maps'))            
        mem = 0;
    end
    for mm = 1:length(mnames) % for every mapping method
        mname = mnames{mm};

        % balanced map
        params = bsols(mm,:); % balanced solution
        if mem
            ratemap = maps_for_plotting{mm,1};
        else                    
            rmset = mapset; % rate mapper settings structure - passing mapset directly to ratemapper causes a memory leak
            rmset.method = mname;
            rmset.binsize = params(1);
            rmset.ssigma = params(2);
            rmset.ash = params(2);            
            rmset.maplims = [-max(abs(epoly(:,1))) -max(abs(epoly(:,2))) max(abs(epoly(:,1))) max(abs(epoly(:,2)))]; 
            rmset.pot = ppot;
            rmset.spt = pspt;
            if strcmp(mname,'fyhn')
                rmset.maxdist = 320;
                rmset.mindist = 50;                                
            end
            [ratemap,dwellmap,spikemap,rmset,speedlift] = rate_mapper(pos,spk,rmset);
            maps_for_plotting{mm,1} = ratemap;
        end

        ax1 = axes('Units','pixels','Position',[xvec(mm) yvec(mm) siz siz],'Color','none','Clipping','off');         
            im = imagesc(ratemap,'alphadata',~isnan(ratemap)); hold on;
            daspect([1 1 1])
            axis xy off tight
            colormap(gca,'turbo'); 
            clim([0 prctile(ratemap(:),99)])

            text(0,1.1,sprintf('%s',sprintf(dnames{mm})),'HorizontalAl','left','Units','normalized','FontSize',8,'rotation',0)                             
            text(0.5,0,sprintf('%.1f, %.1f',params(1),params(2)),'HorizontalAl','center','Units','normalized','FontSize',8,'rotation',0,'Color','w','VerticalAl','bottom')                             
            if mm==1
                text(-0.015,0.5,sprintf('Balanced solution'),'HorizontalAl','center','Units','normalized','FontSize',8,'rotation',90,'Color','k','VerticalAl','bottom')  
                % ah = add_panel_title('c',sprintf('Example rate maps'),'yoffset',5,'xoffset',0,'width',300);                                    
            end

        % minimum error map   
        params = msols(mm,:); % balanced solution
        if mem
            ratemap = maps_for_plotting{mm,2};
        else                    
            rmset = mapset; % rate mapper settings structure - passing mapset directly to ratemapper causes a memory leak
            rmset.method = mname;
            rmset.binsize = params(1);
            rmset.ssigma = params(2);
            rmset.maplims = [-max(abs(epoly(:,1))) -max(abs(epoly(:,2))) max(abs(epoly(:,1))) max(abs(epoly(:,2)))]; 
            rmset.pot = ppot;
            rmset.spt = pspt;
            if strcmp(mname,'fyhn')
                rmset.maxdist = 320;
                rmset.mindist = 50;                                
            end
            [ratemap,dwellmap,spikemap,rmset,speedlift] = rate_mapper(pos,spk,rmset);
            maps_for_plotting{mm,2} = ratemap;
        end

        ax1 = axes('Units','pixels','Position',[xvec(mm) yvec(mm)-siz-spac siz siz],'Color','none','Clipping','off');         
            im = imagesc(ratemap,'alphadata',~isnan(ratemap)); hold on;
            daspect([1 1 1])
            axis xy off tight
            colormap(gca,'turbo'); 
            clim([0 prctile(ratemap(:),99)])
            text(0.5,0,sprintf('%.1f, %.1f',params(1),params(2)),'HorizontalAl','center','Units','normalized','FontSize',8,'rotation',0,'Color','w','VerticalAl','bottom')                             
            if mm==1
                text(-0.015,0.5,sprintf('Minimum error solution'),'HorizontalAl','center','Units','normalized','FontSize',8,'rotation',90,'Color','k','VerticalAl','bottom')                             
            end
    end
    if ~mem
        save(fname,'maps_for_plotting');
    end
         
    % return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ################################################################# %% Save the figure
    % Save the figure    
    disp(sprintf('\tSaving figure...'))    
    figname = [fig_dir2 '\method_summary.png'];
    [~,~,~] = mkdir(fig_dir2);    
    exportgraphics(fig1,figname,'BackgroundColor','w','ContentType','image','Resolution',250);  
    close(fig1)    
    
    
    
    
    
  
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    