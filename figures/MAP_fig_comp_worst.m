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

    sigma_now = '4000';  
    time_now = '16'; 
    mapidx = 4;
    
    % Figure settings
    fig1 = figure('visible','off','Position',[50,0,1200,1050]); 
    set(gcf,'InvertHardCopy','off'); % this stops white lines being plotted black      
    set(gcf,'color','w'); % makes the background colour white
    colormap(jet(256)); % to make sure the colormap is not the horrible default one

%% #################### 
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

    ynow = 645;
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
    xnow = 80;
    xvecp = [xnow xnow xnow];
    yvecp = [ynow ynow-290 ynow-580];

    ms = {'histogram','kadaptive','ksde'};
    ds = {'Histogram','Adaptive smoothing','KSDE'};
    ps = {[9 4.5;3 30;29 230;46 14],[28 128;38 1000;43 5000;67 12500],[31 21;16 64;75 60]};
    rng(999);
    maps_for_plotting = cell(length(ms),10);

    for ii = 1:3
        % main panel showing 16 minute error results
        mname = ms{ii};
        ax1 = axes('Units','pixels','Position',[xvecp(ii) yvecp(ii) xsiz ysiz]);    
            % if ii==1
            %     ah = add_panel_title('b',sprintf('Reported bin size and smoothing combinations'),'yoffset',35,'xoffset',-30,'width',400);
            % end
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

            % % compare smoothing & no smoothing
            % x1 = mat1(1,:);
            % x1a = dme(1,:,:);
            % [~,mindx] = min(x1,[],'all','omitmissing');
            % [i,j] = ind2sub(size(x1),mindx);
            % v1 = reshape(x1a(i,j,:),[],1);
            % 
            % x2 = mat1(2:end,:);
            % x2a = dme(2:end,:,:);
            % [~,mindx] = min(x2,[],'all','omitmissing');
            % [i,j] = ind2sub(size(x2),mindx);
            % v2 = reshape(x2a(i,j,:),[],1);
            % 
            % [~,p,~,s] = ttest2(v1(:),v2(:))
            % [mean(v1(:),'all','omitmissing') mean(v2(:),'all','omitmissing')]
            % keyboard

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
                    text(0,1.14,sprintf(ds{ii}),'Units','normalized','HorizontalAlignment','left','VerticalAlignment','bottom','FontSize',12)                    
                case {'kadaptive'}
                    x = dat_adapt(:,2) + normrnd(0,jitter,size(dat_adapt(:,2)));   
                    y = dat_adapt(:,3) + normrnd(0,jitter,size(dat_adapt(:,2)));  
                    pc2 = plot3(x,y,repmat(1000,size(dat_adapt(:,2))),'k*','MarkerSize',10,'LineStyle','none');   
                    text(0,1.02,sprintf(ds{ii}),'Units','normalized','HorizontalAlignment','left','VerticalAlignment','bottom','FontSize',12)                    
                case {'ksde'}
                    x = dat_ksde(:,2) + normrnd(0,jitter,size(dat_ksde(:,2)));   
                    y = dat_ksde(:,3) + normrnd(0,jitter,size(dat_ksde(:,2)));  
                    pc2 = plot3(x,y,repmat(1000,size(dat_ksde(:,2))),'k*','MarkerSize',10,'LineStyle','none');    
                    text(0,1.02,sprintf(ds{ii}),'Units','normalized','HorizontalAlignment','left','VerticalAlignment','bottom','FontSize',12)                    
            end
            pb = plot3(balance_solution1(1),balance_solution1(2),1000,'ko','MarkerSize',8,'LineStyle','none','MarkerFaceColor','k');            
            pc = plot3(minerr_solution(1),minerr_solution(2),1000,'kv','MarkerSize',8,'LineStyle','none','MarkerFaceColor','k');  

            % if ii==2
            %     [~,leg] = legendflex([pb pc pc2]...
            %         ,{'Balanced solution','Minimum error solution','Literature reported values'}...
            %         ,'anchor',{'n','n'},'ncol',1,'box','off','buffer',[30,80],'xscale',.5,'fontsize',9);  
            % end
    
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
    
            tnow = ps{ii};
            for j = 1:size(tnow,1)
                text(tnow(j,1),tnow(j,2),sprintf('%d',j),'FontSize',10,'HorizontalAlignment','center','VerticalAlignment','middle');
                plot(tnow(j,1),tnow(j,2),'ko','MarkerSize',14,'LineStyle','none')
            end

            if strcmp(mname,'histogram')
                ax1b = axes('Units','pixels','Position',[ax1.Position(1) ax1.Position(2)+ax1.Position(4)+6 ax1.Position(3) 20],'Clipping','off'); 
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

                    text(30,-0.5,sprintf('%d',5),'FontSize',10,'HorizontalAlignment','center','VerticalAlignment','middle','Clipping','off');
                    plot(30,-0.5,'ko','MarkerSize',14,'LineStyle','none','Clipping','off')                    
            end
% if ii==3
%     % keyboard
%     return
% end
%% #################### get the error, time, place field error values for the numbered values
            dmt = data_matrix_cputime;
            dmf = data_matrix_fields;
            dmm = data_matrix_missing; 
            pnow = ps{ii};
            actual_params = NaN(size(pnow)+[3 0]);
            switch mname
                case {'histogram'}
                    dindx = dat;   
                    zidx = dat(:,3)>0;
                case {'kadaptive'}
                    dindx = dat_adapt;     
                    zidx = true(size(dindx,1),1);
                case {'ksde'}     
                    dindx = dat_ksde;   
                    zidx = true(size(dindx,1),1);                    
                otherwise
                    keyboard
            end               

            % get the error, time, place field error values for this worst value
            e = NaN(256,size(pnow,1)+3);
            % p = NaN(256,size(pnow,1)+3);  
            % m = NaN(256,size(pnow,1)+3);                
            % t = NaN(256,size(pnow,1)+3);              
            for jj = 1:size(dme,3) % for every cell
                e(jj,2) = interp2(bmat,smat,dme(:,:,jj),balance_solution1(1),balance_solution1(2));
                e(jj,1) = interp2(bmat,smat,dme(:,:,jj),minerr_solution(1),minerr_solution(2)); 
        
                % p(jj,2) = interp2(bmat,smat,dmf(:,:,jj),balance_solution1(1),balance_solution1(2),'nearest');
                % p(jj,1) = interp2(bmat,smat,dmf(:,:,jj),minerr_solution(1),minerr_solution(2),'nearest'); 
                % 
                % m(jj,2) = interp2(bmat,smat,dmm(:,:,jj),balance_solution1(1),balance_solution1(2));
                % m(jj,1) = interp2(bmat,smat,dmm(:,:,jj),minerr_solution(1),minerr_solution(2));                 
                % 
                % t(jj,2) = interp2(bmat,smat,dmt(:,:,jj),balance_solution1(1),balance_solution1(2));            
                % t(jj,1) = interp2(bmat,smat,dmt(:,:,jj),minerr_solution(1),minerr_solution(2));      

                actual_params(1,:) = [minerr_solution(1),minerr_solution(2)];
                actual_params(2,:) = [balance_solution1(1),balance_solution1(2)];
                dindx2 = dindx(zidx,:);
                for pp = 1:size(pnow,1)
                    c_idx = knnsearch(dindx2(:,2:3),pnow(pp,:));
                    a_point = dindx2(c_idx,2:3);
                    actual_params(pp+2,:) = a_point;
                    e(jj,pp+2) = interp2(bmat,smat,dme(:,:,jj),a_point(1),a_point(2));
                    % p(jj,pp+2) = interp2(bmat,smat,dmf(:,:,jj),a_point(1),a_point(2),'nearest');
                    % m(jj,pp+2) = interp2(bmat,smat,dmm(:,:,jj),a_point(1),a_point(2));
                    % t(jj,pp+2) = interp2(bmat,smat,dmt(:,:,jj),a_point(1),a_point(2));            
                end
                if strcmp(mname,'histogram') % add zero smoothing value (histogram only)
                    c_idx = knnsearch(dindx(:,2:3),[30 0]);
                    a_point = dindx(c_idx,2:3);
                    actual_params(end,:) = a_point;                    
                    e(jj,end) = interp2(bmat,smat,dme(:,:,jj),a_point(1),a_point(2));
                    % p(jj,end) = interp2(bmat,smat,dmf(:,:,jj),a_point(1),a_point(2),'nearest');
                    % m(jj,end) = interp2(bmat,smat,dmm(:,:,jj),a_point(1),a_point(2));
                    % t(jj,end) = interp2(bmat,smat,dmt(:,:,jj),a_point(1),a_point(2)); 
                end                
            end


            % plot the error results
            xnow = ax1.Position(1)+270;
            ax = axes('Units','pixels','Position',[xnow ax1.Position(2) ax1.Position(3) 180],'Clipping','off'); 
                switch mname
                    case {'histogram'}
                        [ds1,gs1] = vectorDATAGROUP([],e(:,1),e(:,2),e(:,3),e(:,4),e(:,5),e(:,6),e(:,7));
                        cols = winter(7);
                        ax.XLim = [.5 7.5];
                        ax.XTick = [1:7];
                        ax.XTickLabel = {'Minimum','Balanced','1','2','3','4','5'};
                    case {'kadaptive'}
                        [ds1,gs1] = vectorDATAGROUP([],e(:,1),e(:,2),e(:,3),e(:,4),e(:,5),e(:,6));
                        cols = winter(7);
                        ax.XLim = [.5 6.5];
                        ax.XTick = [1:6];
                        ax.XTickLabel = {'Minimum','Balanced','1','2','3','4'};
                    case {'ksde'}     
                        [ds1,gs1] = vectorDATAGROUP([],e(:,1),e(:,2),e(:,3),e(:,4),e(:,5));
                        cols = winter(7);
                        ax.XLim = [.5 5.5];
                        ax.XTick = [1:5];
                        ax.XTickLabel = {'Minimum','Balanced','1','2','3'};                   
                    otherwise
                        keyboard
                end      
                cmat = mat2cell(cols,ones(size(cols,1),1));                
                meanplot(ds1,gs1,'linecolor',{'none'},'meancolor',{'k'},'meanmarker',{'o'},'meansize',6,'dots',1,'dotsize',8,'dotalpha',0.5,'dotsigma',0.05,'dotcolor',cmat,'meanlinewidth',1.5); hold on;
        % keyboard
                ax = gca;
                ax.YScale = 'log';
                ylabel(sprintf('MISE'))  
                text(0.99,1.1,'N = 256 cells','Units','normalized','HorizontalAlignment','right','VerticalAlignment','top')            

        ax_sig = axes('Units','pixels','Position',ax.Position,'Color','none'); 
            ax_sig.XLim = ax.XLim;
            ax_sig.YLim = ax.YLim;            
            axis off
            [result1,result2] = plotsigbrackets(ds1,gs1,'bracket_text_fsize',0,'bracket_y_gap_coeff',55);

%% #################### plot the example ratemaps
            % cut the position and spike data to the trial length
            % we want
            excell = 1;
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

            siz = 100;
            xnow = ax.Position(1)+240;
            xvec = xnow : siz+10 : 10000;
            xvec = [xvec(1:2) xvec(1:5)];
            ynow = ax.Position(2);
            yvec = [ones([1 2])*ynow+siz+10 ones([1 5])*ynow-10];
        
            if ii==1
                ax1 = axes('Units','pixels','Position',[xvec(1) yvec(1) siz siz],'Color','none','Clipping','off');         
                    p1 = plot(ppox,ppoy,'Color',[.5 .5 .5 .5]); hold on;
                    s1 = scatter(pspx,pspy,5,'r','filled','MarkerEdgeColor','none','MarkerFaceAlpha',0.5);
                    daspect([1 1 1])
                    axis xy off  
            
                    text(0,1.05,sprintf('Positions and spikes'),'HorizontalAl','left','Units','normalized','FontSize',8,'rotation',0,'Color','k','VerticalAl','bottom')  
            
                    [~,leg] = legendflex([p1 s1]...
                        ,{'trajectory','spikes'}...
                        ,'anchor',{'e','e'},'ncol',1,'box','off','buffer',[90,0],'xscale',0.5,'fontsize',9); 

                axc = axes('Units','pixels','Position',[ax1.Position(1)+450 ax1.Position(2)+20 12 ax1.Position(4)-20]); 
                    mat = linspace(0,100,100)';
                    imagesc(ones(size(mat)),mat,mat);
                    colormap(axc,'turbo');
                    axis xy
        
                    axc.YTick = [0 100];
                    axc.YTickLabel = {'0 Hz','Max'};
                    axc.XTick = [];
                    axc.YAxisLocation = 'right';
                    text(0,1.25,sprintf('Firing\nrate (Hz)'),'FontSize',8,'HorizontalAl','left','Units','normalized')
                    % add univisted box
                    axc2 = axes('Units','pixels','Position',[axc.Position(1) axc.Position(2)-20 axc.Position(3) 10]); 
                    axis on; box on; axc2.XTick = []; axc2.YTick = 0.5; axc2.YTickLabel = {'Unvisited'}; axc2.YAxisLocation = 'right';


                xvec = xnow : siz+10 : 10000;
                xvec = [xvec(3:4) xvec(1:5)];
            end

            fname = [fig_dir '_comp_worst_maps.mat'];
            override_maps = 0;
            if exist(fname,'file') && ~override_maps
                mem = 1;
                disp(sprintf('\t\t...loading old maps'))
                load(fname,'maps_for_plotting')
            else
                disp(sprintf('\t\t...generating maps'))            
                mem = 0;
            end

            titles = {'Minimum','Balanced','1','2','3','4','5'};
            for mm = 1:size(actual_params,1) % for every set of map settings        
                params = actual_params(mm,:);
                if all(isnan(params(:)))
                    continue
                end
                if mem
                    ratemap = maps_for_plotting{ii,mm};
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
                    maps_for_plotting{ii,mm} = ratemap;
                end
        
                ax1 = axes('Units','pixels','Position',[xvec(mm) yvec(mm) siz siz],'Color','none','Clipping','off');         
                    im = imagesc(ratemap,'alphadata',~isnan(ratemap)); hold on;
                    daspect([1 1 1])
                    axis xy off tight
                    colormap(gca,'turbo'); 
                    clim([0 prctile(ratemap(:),99)])
        
                    text(0,1,sprintf('%s',titles{mm}),'HorizontalAl','left','Units','normalized','FontSize',10,'rotation',0,'Color','k','VerticalAl','bottom')                             
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
    figname = [fig_dir2 '\S10_Fig.png'];
    [~,~,~] = mkdir(fig_dir2);    
    exportgraphics(fig1,figname,'BackgroundColor','w','ContentType','image','Resolution',350);  
    close(fig1)    
   
    
    
    
    
    
    
    






