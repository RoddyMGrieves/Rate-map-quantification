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
dispersion = [0 0.8];

if 0
%% #################### Example fields

            sdir = [scratch_space '\' ename '_sigma16000_duration64' ovals{3} '_sdata.mat'];            
            disp(sprintf('\t\t...loading %s',sdir));            
            load(sdir,'sdata'); 

            % cut the position and spike data to the trial length
            % we want
            excell = 5;
            pindx = sdata.pos_index(excell);
            pos_now = all_pos{pindx,1};
            pox = pos_now(:,1);
            poy = pos_now(:,2);
            pot = pos_now(:,3); 

            spk = sdata.spk{excell};
            spx = spk(:,1);
            spy = spk(:,2);
            spt = spk(:,3);

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
            rmset.method = 'histogram';
            rmset.binsize = 20;
            rmset.ssigma = 40;
            rmset.maplims = [-max(abs(epoly(:,1))) -max(abs(epoly(:,2))) max(abs(epoly(:,1))) max(abs(epoly(:,2)))]; 
            rmset.pot = ppot;
            rmset.spt = pspt;
            [ratemap,dwellmap,spikemap,rmset,speedlift] = rate_mapper(pos,spk,rmset);

            [z,overdispersion,r] = getOVERDISPERSE(rmset.map_pos,ppot,pspt,ratemap);

figure
xi = -10:0.5:10;
% f = ksdensity(z,xi,'bandwidth',1,'kernel','normal');
% area(xi,f,'FaceColor','k','FaceAlpha',0.2,'EdgeColor','none'); 
f = histcounts(z(:),xi,'Normalization','probability');
x = movmean(xi,2,'EndPoints','discard');
bar(x,f,1,'k');
hold on;

xi = -10:0.1:10;
y = normpdf(-10:0.1:10,0,1);
y = y ./ nanmax(y(:)) .* nanmax(f(:));
plot(xi,y,'r'); 

% return



            % keyboard
    poxmap = rmset.map_pos(:,1);
    poymap = rmset.map_pos(:,2);
    inds = sub2ind(size(ratemap),round(poymap),round(poxmap));
    exp_frate = ratemap(inds);

    log_field = exp_frate > (max(ratemap(:)).*0.2);
    n_pass = bwlabel(log_field);

    dv = 50;
    m = NaN(max(n_pass),dv);

    inst_spikes = histcounts(spt,min(ppot(:))-0.01:0.02:max(ppot(:))+0.01);
    z = NaN(max(n_pass),1);
    v = NaN(max(n_pass),2);
    for pp = 1:max(n_pass)
        Nobs = sum(inst_spikes(n_pass==pp));
        Nexp = sum(exp_frate(n_pass==pp).*(1/50));
        z(pp) = (Nobs-Nexp)/sqrt(Nexp);
        v(pp,:) = [Nobs Nexp];

        i = inst_spikes(n_pass==pp);
        if length(i)>dv
            m(pp,:) = i(1:dv);
        else
            m(pp,1:length(i)) = i;

        end
    end


figure
xi = -10:0.5:10;
% f = ksdensity(z,xi,'bandwidth',1,'kernel','normal');
% area(xi,f,'FaceColor','k','FaceAlpha',0.2,'EdgeColor','none'); 
f = histcounts(z(:),xi,'Normalization','probability');
x = movmean(xi,2,'EndPoints','discard');
bar(x,f,1,'k');
hold on;

xi = -10:0.1:10;
y = normpdf(-10:0.1:10,0,1);
y = y ./ nanmax(y(:)) .* nanmax(f(:));
plot(xi,y,'r'); 

% return
% figure
% imagesc(m)
% 
% 
% keyboard


figure
subplot(2,2,1)
plot(ppox,ppoy,'k'); hold on;
plot(pspx,pspy,'r.','MarkerSize',20);
daspect([1 1 1])

subplot(2,2,2)
imagesc(ratemap)
daspect([1 1 1])
axis xy


return
subplot(2,2,3)
xi = -10:0.5:10;
% f = ksdensity(z,xi,'bandwidth',1,'kernel','normal');
% area(xi,f,'FaceColor','k','FaceAlpha',0.2,'EdgeColor','none'); 
f = histcounts(z(:),xi,'Normalization','probability');
x = movmean(xi,2,'EndPoints','discard');
bar(x,f,1,'k');
hold on;

xi = -10:0.1:10;
y = normpdf(-10:0.1:10,0,1);
y = y ./ nanmax(y(:)) .* nanmax(f(:));
plot(xi,y,'r'); 






















return
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ################################################################# %% FUNCTION BODY
%% #################### Heatmap of error by binsize and smoothing
    % Figure settings
    fig1 = figure('visible','on','Position',[50,60,1000,800]); 
    set(gcf,'InvertHardCopy','off'); % this stops white lines being plotted black      
    set(gcf,'color','w'); % makes the background colour white
    colormap(jet(256)); % to make sure the colormap is not the horrible default one

    % method settings
    xnow = 70;
    ynow = 600;
    mnames = {'histogram','ash','kadaptive','ksde','fyhn','kyadaptive'};
mnames = {'histogram','ash'}
    
    dnames = {'Histogram','ASH','Adaptive smoothing','KSDE','tKSDE','Adaptive binning'}; 
    ls = {'a','b','c','d','e','f','g'};
    mapidx = 4;
    xsiz = 150;
    xbuff = 50;
    xvec = xnow : (xsiz+xbuff) : 1000;
    ysiz = 150;
    ybuff = 70;
    yvec = ynow : -(ysiz+ybuff) : 0;
    [xx,yy] = meshgrid(xvec,yvec);
    xx = xx';
    yy = yy;

    % plot settings
    n_color_levels = 64;
    error_colormap = flipud(inferno(n_color_levels));    
    clims = [1*10^-12 1*10^-11];   
    interp_method = 'nearest';  
  
    for mm = 1:length(mnames) % for every method
        disp(sprintf('\t%s...',mnames{mm}))
        
        for oo = 1:length(ovals) % for every overdispersion value
            mname = mnames{mm};
            sigma_now = '16000';
            time_now = '4'; 
            disp(sprintf('\t\t...%s',ovals{oo}))
    
            % main panel showing 16 minute error results
            ax1 = axes('Units','pixels','Position',[xx(oo) yy(mm) xsiz ysiz]);  
                disp(sprintf('\t\t...errormap'))

                load([scratch_space '\' ename '_' mname '_sigma' num2str(sigma_now) '_duration' num2str(time_now) ovals{oo} '_datamats.mat'],'data_matrix_error','data_matrix_params')  
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
                        text(0,1.15,'4 minutes coverage','HorizontalAl','left','Units','normalized','FontSize',8);                
                    otherwise
                        text(0,1.04,'4 minutes coverage','HorizontalAl','left','Units','normalized','FontSize',8);
                end            
                ylabel(smoo_text);
                xlabel(grid_text);

%% #################### calculate Pareto front
            if 1
                disp(sprintf('\t\t...pareto front'))

                fname = [scratch_space '\' ename '_' mname '_sigma' num2str(sigma_now) '_duration' num2str(time_now) ovals{oo} '_gamultiobj.mat'];
                if exist(fname,'file') && 1 && ~override_pareto
                    continue
                else
                    load([scratch_space '\' ename '_' mname '_sigma' num2str(sigma_now) '_duration' num2str(time_now) ovals{oo} '_datamats.mat'],'data_matrix_cputime','data_matrix_fields','data_matrix_error','data_matrix_missing','data_matrix_params')   
                end
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

%% #################### Plot Solutions onto error map 
            if 1     
                disp(sprintf('\t\t...plotting solutions'))
                
                load([scratch_space '\' ename '_' mname '_sigma' num2str(sigma_now) '_duration' num2str(time_now) ovals{oo} '_gamultiobj.mat'],'balance_solution2','balance_solution1','minerr_solution','x');
                % find the 'best' solution (closest to minimising everything)
                pb1 = plot3(balance_solution1(1),balance_solution1(2),1000,'ko','MarkerSize',10,'LineStyle','none');            

                % find the minimal error solution 
                pc = plot3(minerr_solution(1),minerr_solution(2),1000,'kv','MarkerSize',10,'LineStyle','none');          

                % [~,leg] = legendflex([pb1 pc]...
                %     ,{'Balanced solution','Minimum error solution'}...
                %     ,'anchor',{'se','se'},'ncol',1,'box','off','buffer',[-30,-90],'xscale',.5,'fontsize',9);    
            end
        end
    end

% keyboard
   





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ################################################################# %% FUNCTION BODY
%% #################### Heatmap of error by binsize and smoothing
    % Figure settings
    fig1 = figure('visible','on','Position',[50,60,900,800]); 
    set(gcf,'InvertHardCopy','off'); % this stops white lines being plotted black      
    set(gcf,'color','w'); % makes the background colour white
    colormap(jet(256)); % to make sure the colormap is not the horrible default one
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

    xnow = 80;
    ynow = 550;

    % error results
    ax1 = axes('Units','pixels','Position',[xnow ynow-5 300 170]); 
        ah = add_panel_title('a',sprintf('Error across mapping methods'),'yoffset',5,'xoffset',-30,'width',300);                    
   
        % plot results
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
        ax.XLim = [0 1];
        ax.XTick = dispersion;
        ax.YTick = (ones(1,20).*10).^-[20:-1:1];        
        % ax.YLim = [-1*10^-12 16*10^-12]; 
        % ax.YTick = 0:1:5;
        % ax.XTickLabel = [];
        ax.YScale = 'log'; 
        ax.YLim(1) = 10^-13;
        ax.YLim(2) = 10^-11;        
        ylabel(sprintf('MISE')) 
        xlabel('Dispersion (a_{mod})')
        text(0.99,1.1,'N = 256 cells','Units','normalized','HorizontalAlignment','right','VerticalAlignment','top')

    % field detection results
    ax1 = axes('Units','pixels','Position',[xnow+400 ynow+10 300 170-10]); 
        ah = add_panel_title('b',sprintf('Computation time across mapping methods'),'yoffset',0,'xoffset',-30,'width',300);                    
   
        % plot results
        cols = winter(length(mnames));
        vals = NaN(length(mnames),length(ovals)); % 
        for mm = 1:length(mnames) % for every mapping method
            enow = squeeze(p_vals(:,mm,1,:)); % balanced
            vals(mm,:) = sum(enow==0,1,'omitnan') ./ size(enow,1);

            plot(dispersion,vals(mm,:),'Color',cols(mm,:),'Marker','o'); hold on;
        end

        ax = gca;
        ax.XLim = [0 1];
        ax.XTick = dispersion;
        ax.YLim = [0.8 1]; 
        ylabel(sprintf('Prop. correct place field count'))  
        xlabel('Dispersion (a_{mod})')        
        text(0.99,1.1,'N = 256 cells','Units','normalized','HorizontalAlignment','right','VerticalAlignment','top')
        box off



% keyboard


        % 
        % [g,~] = meshgrid(1:size(err_vals,2),1:size(err_vals,1),1:size(err_vals,4));
        % 
        % 
        % keyboard
        % gs = g(:);
        % ds = reshape(err_vals(:,:,1,:),[],1);
        % cols = winter(length(mnames));
        % cmat = mat2cell(cols,ones(size(cols,1),1));
        % meanplot(ds,gs-0.2,'linecolor',{'none'},'meancolor',{'k'},'meanmarker',{'o'},'meansize',6,'dots',1,'dotsize',8,'dotalpha',0.5,'dotsigma',0.05,'dotcolor',cmat,'meanlinewidth',1.5); hold on;
        % 
        % ds = reshape(err_vals(:,:,2),[],1);
        % cols = winter(length(mnames));
        % cmat = mat2cell(cols,ones(size(cols,1),1));
        % meanplot(ds,gs+0.2,'linecolor',{'none'},'meancolor',{'k'},'meanmarker',{'v'},'meansize',6,'dots',1,'dotsize',8,'dotalpha',0.5,'dotsigma',0.05,'dotcolor',cmat,'meanlinewidth',1.5); hold on;
        % 
        % x = [(mean(g,1,'omitnan')-0.2)' (mean(g,1,'omitnan')+0.2)'];
        % y = [mean(err_vals(:,:,1),1,'omitnan')' mean(err_vals(:,:,2),1,'omitnan')'];
        % plot(x',y','k');
        % 
        % ax = gca;
        % ax.XLim = [.5 length(mnames)+.5];
        % ax.XTick = [1:length(mnames)];
        % ax.YTick = (ones(1,20).*10).^-[20:-1:1];        
        % % ax.YLim = [-1*10^-12 16*10^-12]; 
        % % ax.YTick = 0:1:5;
        % ax.XTickLabel = [];
        % ax.YScale = 'log'; 
        % ax.YLim(1) = 10^-13;
        % ylabel(sprintf('MISE'))  
        % text(0.99,1.1,'N = 256 cells','Units','normalized','HorizontalAlignment','right','VerticalAlignment','top')
        % 
        % for ii = 1:length(mnames)
        %     text(ii,ax.YLim(1),sprintf(dnames{ii}),'VerticalAlignment','top','HorizontalAlignment','center','FontSize',8)
        % end
        % 
        % axis manual
        % p1 = plot(-10,1,'ko','LineWidth',1.5);
        % p2 = plot(-10,1,'kv','LineWidth',1.5);
        % [~,leg] = legendflex([p1 p2],{'Balanced solution','Minimum error solution'},'anchor',{'nw','nw'},'ncol',1,'box','off','buffer',[0,20],'xscale',.5,'fontsize',9); 
        % 












return
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% ################################################################# %% Save the figure
%     % Save the figure    
%     disp(sprintf('\tSaving figure...'))    
%     figname = [fig_dir2 '\v2_fig_error_matrix.png'];
%     [~,~,~] = mkdir(fig_dir2);    
%     exportgraphics(fig1,figname,'BackgroundColor','w','ContentType','image','Resolution',250);  
%     close(fig1)    
    
    


  
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    