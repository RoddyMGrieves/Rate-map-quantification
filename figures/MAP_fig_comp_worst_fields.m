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
    fig1 = figure('visible','on','Position',[50,0,1200,1050]); 
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

    ynow = 635;
    xsiz = 200;
    ysiz = 200;

    disp(sprintf('\tMISE map...'))
    n_color_levels = 64;
    error_colormap = parula(n_color_levels);    
    clims = [0 0.5];
    interp_method = 'nearest';  
    mapidx = 4;
    transp = 0.5;
    xbuff = xsiz+70;
    xnow = 80;
    xvecp = [xnow xnow xnow];
    yvecp = [ynow ynow-280 ynow-560];

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
            load([scratch_space '\' ename '_' mname '_sigma' num2str(sigma_now) '_duration' num2str(time_now) '_datamats.mat'],'data_matrix_error','data_matrix_fields','data_matrix_cputime','data_matrix_params','data_matrix_missing')    
            % dme = data_matrix_error(:,:,:,mapidx);
            dme = data_matrix_fields;
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
            % e = NaN(256,size(pnow,1)+3);
            p = NaN(256,size(pnow,1)+3);  
            m = NaN(256,size(pnow,1)+3);                
            t = NaN(256,size(pnow,1)+3);              
            for jj = 1:size(dme,3) % for every cell
                % e(jj,2) = interp2(bmat,smat,dme(:,:,jj),balance_solution1(1),balance_solution1(2));
                % e(jj,1) = interp2(bmat,smat,dme(:,:,jj),minerr_solution(1),minerr_solution(2)); 
        
                p(jj,2) = interp2(bmat,smat,dmf(:,:,jj),balance_solution1(1),balance_solution1(2),'nearest');
                p(jj,1) = interp2(bmat,smat,dmf(:,:,jj),minerr_solution(1),minerr_solution(2),'nearest'); 

                m(jj,2) = interp2(bmat,smat,dmm(:,:,jj),balance_solution1(1),balance_solution1(2));
                m(jj,1) = interp2(bmat,smat,dmm(:,:,jj),minerr_solution(1),minerr_solution(2));                 

                t(jj,2) = interp2(bmat,smat,dmt(:,:,jj),balance_solution1(1),balance_solution1(2));            
                t(jj,1) = interp2(bmat,smat,dmt(:,:,jj),minerr_solution(1),minerr_solution(2));      

                actual_params(1,:) = [minerr_solution(1),minerr_solution(2)];
                actual_params(2,:) = [balance_solution1(1),balance_solution1(2)];
                dindx2 = dindx(zidx,:);
                for pp = 1:size(pnow,1)
                    c_idx = knnsearch(dindx2(:,2:3),pnow(pp,:));
                    a_point = dindx2(c_idx,2:3);
                    actual_params(pp+2,:) = a_point;
                    % e(jj,pp+2) = interp2(bmat,smat,dme(:,:,jj),a_point(1),a_point(2));
                    p(jj,pp+2) = interp2(bmat,smat,dmf(:,:,jj),a_point(1),a_point(2),'nearest');
                    m(jj,pp+2) = interp2(bmat,smat,dmm(:,:,jj),a_point(1),a_point(2));
                    t(jj,pp+2) = interp2(bmat,smat,dmt(:,:,jj),a_point(1),a_point(2));            
                end
                if strcmp(mname,'histogram') % add zero smoothing value (histogram only)
                    c_idx = knnsearch(dindx(:,2:3),[30 0]);
                    a_point = dindx(c_idx,2:3);
                    actual_params(end,:) = a_point;                    
                    % e(jj,end) = interp2(bmat,smat,dme(:,:,jj),a_point(1),a_point(2));
                    p(jj,end) = interp2(bmat,smat,dmf(:,:,jj),a_point(1),a_point(2),'nearest');
                    m(jj,end) = interp2(bmat,smat,dmm(:,:,jj),a_point(1),a_point(2));
                    t(jj,end) = interp2(bmat,smat,dmt(:,:,jj),a_point(1),a_point(2)); 
                end                
            end


%% #################### plot the field error results
            xnow = ax1.Position(1)+270;
            ax = axes('Units','pixels','Position',[xnow ax1.Position(2) ax1.Position(3) 180],'Clipping','off'); 
                switch mname
                    case {'histogram'}
                        [ds1,gs1] = vectorDATAGROUP([],p(:,1),p(:,2),p(:,3),p(:,4),p(:,5),p(:,6),p(:,7));
                        cols = winter(7);
                        ax.XLim = [.5 7.5];
                        ax.XTick = [1:7];
                        ax.XTickLabel = {'Minimum','Balanced','1','2','3','4','5'};
                    case {'kadaptive'}
                        [ds1,gs1] = vectorDATAGROUP([],p(:,1),p(:,2),p(:,3),p(:,4),p(:,5),p(:,6));
                        cols = winter(7);
                        ax.XLim = [.5 6.5];
                        ax.XTick = [1:6];
                        ax.XTickLabel = {'Minimum','Balanced','1','2','3','4'};
                    case {'ksde'}     
                        [ds1,gs1] = vectorDATAGROUP([],p(:,1),p(:,2),p(:,3),p(:,4),p(:,5));
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
                % ax.YScale = 'log';
                ylabel(sprintf('Field error (fields)'))  
                text(0.99,0,'N = 256 cells','Units','normalized','HorizontalAlignment','right','VerticalAlignment','bottom')            

            ax_sig = axes('Units','pixels','Position',ax.Position,'Color','none'); 
                ax_sig.XLim = ax.XLim;
                ax_sig.YLim = ax.YLim;            
                axis off
                [result1,result2] = plotsigbrackets(ds1,gs1,'bracket_text_fsize',0,'bracket_y_gap_coeff',55,'test','kw');


%% #################### plot the comp time results
            xnow = ax.Position(1)+270;
            ax = axes('Units','pixels','Position',[xnow ax1.Position(2) ax1.Position(3) 180],'Clipping','off'); 
                switch mname
                    case {'histogram'}
                        [ds1,gs1] = vectorDATAGROUP([],t(:,1),t(:,2),t(:,3),t(:,4),t(:,5),t(:,6),t(:,7));
                        cols = winter(7);
                        ax.XLim = [.5 7.5];
                        ax.XTick = [1:7];
                        ax.XTickLabel = {'Minimum','Balanced','1','2','3','4','5'};
                    case {'kadaptive'}
                        [ds1,gs1] = vectorDATAGROUP([],t(:,1),t(:,2),t(:,3),t(:,4),t(:,5),t(:,6));
                        cols = winter(7);
                        ax.XLim = [.5 6.5];
                        ax.XTick = [1:6];
                        ax.XTickLabel = {'Minimum','Balanced','1','2','3','4'};
                    case {'ksde'}     
                        [ds1,gs1] = vectorDATAGROUP([],t(:,1),t(:,2),t(:,3),t(:,4),t(:,5));
                        cols = winter(7);
                        ax.XLim = [.5 5.5];
                        ax.XTick = [1:5];
                        ax.XTickLabel = {'Minimum','Balanced','1','2','3'};                   
                    otherwise
                        keyboard
                end      
                cmat = mat2cell(cols,ones(size(cols,1),1));                
                meanplot(ds1,gs1,'linecolor',{'none'},'meancolor',{'k'},'meanmarker',{'o'},'meansize',6,'dots',1,'dotsize',8,'dotalpha',0.5,'dotsigma',0.05,'dotcolor',cmat,'meanlinewidth',1.5); hold on;
        
                ax = gca;
                % ax.YScale = 'log';
                ylabel(sprintf('Computation time (ms)'))  

            ax_sig = axes('Units','pixels','Position',ax.Position,'Color','none'); 
                ax_sig.XLim = ax.XLim;
                ax_sig.YLim = ax.YLim;            
                axis off
                [result1,result2] = plotsigbrackets(ds1,gs1,'bracket_text_fsize',0,'bracket_y_gap_coeff',55,'test','kw');

%% #################### plot the missing values results
            xnow = ax.Position(1)+270;
            ax = axes('Units','pixels','Position',[xnow ax1.Position(2) ax1.Position(3) 180],'Clipping','off'); 
                switch mname
                    case {'histogram'}
                        [ds1,gs1] = vectorDATAGROUP([],m(:,1),m(:,2),m(:,3),m(:,4),m(:,5),m(:,6),m(:,7));
                        cols = winter(7);
                        ax.XLim = [.5 7.5];
                        ax.XTick = [1:7];
                        ax.XTickLabel = {'Minimum','Balanced','1','2','3','4','5'};
                    case {'kadaptive'}
                        [ds1,gs1] = vectorDATAGROUP([],m(:,1),m(:,2),m(:,3),m(:,4),m(:,5),m(:,6));
                        cols = winter(7);
                        ax.XLim = [.5 6.5];
                        ax.XTick = [1:6];
                        ax.XTickLabel = {'Minimum','Balanced','1','2','3','4'};
                    case {'ksde'}     
                        [ds1,gs1] = vectorDATAGROUP([],m(:,1),m(:,2),m(:,3),m(:,4),m(:,5));
                        cols = winter(7);
                        ax.XLim = [.5 5.5];
                        ax.XTick = [1:5];
                        ax.XTickLabel = {'Minimum','Balanced','1','2','3'};                   
                    otherwise
                        keyboard
                end      
                cmat = mat2cell(cols,ones(size(cols,1),1));                
                meanplot(ds1,gs1,'linecolor',{'none'},'meancolor',{'k'},'meanmarker',{'o'},'meansize',6,'dots',1,'dotsize',8,'dotalpha',0.5,'dotsigma',0.05,'dotcolor',cmat,'meanlinewidth',1.5); hold on;
        
                ax = gca;
                % ax.YScale = 'log';
                ylabel(sprintf('Empty bins (prop.)'))  

            ax_sig = axes('Units','pixels','Position',ax.Position,'Color','none'); 
                ax_sig.XLim = ax.XLim;
                ax_sig.YLim = ax.YLim;            
                axis off
                [result1,result2] = plotsigbrackets(ds1,gs1,'bracket_text_fsize',0,'bracket_y_gap_coeff',55,'test','kw');


% keyboard
    end

    % return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ################################################################# %% Save the figure
    % Save the figure    
    disp(sprintf('\tSaving figure...'))    
    figname = [fig_dir2 '\S11_Fig.png'];
    [~,~,~] = mkdir(fig_dir2);    
    exportgraphics(fig1,figname,'BackgroundColor','w','ContentType','image','Resolution',350);  
    close(fig1)    
   
    
    
    
    
    
    
    






