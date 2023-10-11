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

%     if 0
%         pos = all_pos{1,1}(:,[1 2]);
%         N = numel(pos(:,1));
%         
%         % sqrt N rule   
%         nbins = floor(sqrt(N));
%         sqrt_n_binsize = 1200 / nbins
% 
%         % Sturges rule    
%         nbins = ceil(log2(N) +1);
%         sturges_binsize = 1200 / nbins
% 
%         % Freedman–Diaconis rule
%         IQRs = iqr(pos);              % Get the IQ ranges across each dimension
%         Hs = 2* IQRs* N^(-1/3);     % Get the bandwidths across each dimension
%         Ranges = range(pos);          % Get the range of values across each dimension
%         % Get the suggested number of bins along each dimension
%         nbins_dim1 = ceil(Ranges(1)/Hs(1)); % 12 here
%         nbins_dim2 = ceil(Ranges(2)/Hs(2)); % 15 here
%         % Get the maximum of the two
%         nbins = max( [nbins_dim1, nbins_dim2]); 
%         freedman1_binsize = 1200 / nbins
%         
%         % Freedman–Diaconis on the norms
%         Norms = sqrt(sum(pos.^2,2));        % Get the norm of each point in th 2-D sample
%         H_norms = 2* iqr(Norms)* N^(-1/3);% Get the "norm" bandwidth
%         nbins = ceil(range(Norms)/ H_norms);   % Get number of bins         
%         freedman2_binsize = 1200 / nbins
% 
%         pos = all_pos{1,1};        
%         pos2 = pos(pos(:,3)<16*60,[1 2]);
%         N = numel(pos2(:,1));
%         
%         % sqrt N rule   
%         nbins = floor(sqrt(N));
%         sqrt_n_binsize = 1200 / nbins
% 
%         % Sturges rule    
%         nbins = ceil(log2(N) +1);
%         sturges_binsize = 1200 / nbins
% 
%         % Freedman–Diaconis rule
%         IQRs = iqr(pos2);              % Get the IQ ranges across each dimension
%         Hs = 2* IQRs* N^(-1/3);     % Get the bandwidths across each dimension
%         Ranges = range(pos2);          % Get the range of values across each dimension
%         % Get the suggested number of bins along each dimension
%         nbins_dim1 = ceil(Ranges(1)/Hs(1)); % 12 here
%         nbins_dim2 = ceil(Ranges(2)/Hs(2)); % 15 here
%         % Get the maximum of the two
%         nbins = max( [nbins_dim1, nbins_dim2]); 
%         freedman1_binsize = 1200 / nbins
%         
%         % Freedman–Diaconis on the norms
%         Norms = sqrt(sum(pos2.^2,2));        % Get the norm of each point in th 2-D sample
%         H_norms = 2* iqr(Norms)* N^(-1/3);% Get the "norm" bandwidth
%         nbins = ceil(range(Norms)/ H_norms);   % Get number of bins         
%         freedman2_binsize = 1200 / nbins        
%         return
%     end

    ename = 'arena120cm'; 
    override_pareto = 0;
    override_maps = 0;
    rmap_colormap = 'turbo';

%% ##################################### Heading 2
%% #################### Heading 3
%% Heading 4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ################################################################# %% FUNCTION BODY
    if 1
        % Figure settings
        fig1 = figure('visible','on','Position',[50,60,1550,900]); 
        set(gcf,'InvertHardCopy','off'); % this stops white lines being plotted black      
        set(gcf,'color','w'); % makes the background colour white
        colormap(jet(256)); % to make sure the colormap is not the horrible default one

%% #################### Example cells
        ynow = 750;
        xnow = 80;
        if 1
            sdir = [scratch_space '\' ename '_sigma8000_duration64_sdata.mat'];            
            disp(sprintf('\t\t...loading %s',sdir));            
            load(sdir,'sdata'); 

            disp(sprintf('\tExample rate maps...'))

            psiz = 42;
            buff = 47;
            xvec = xnow:buff:10000;
            yvec = ynow:-buff:-1000;
            switch mname
                case {'histogram'}
                    bvnow = [1 10 20 40 80 160 320];
                    svnow = [1 10 20 40 80 160 320];
                    smoo_text = sprintf('Smoothing %c (mm)',963);
                    grid_text = sprintf('Bin size (mm)');
                case {'ash'}
                    bvnow = [20 40 60 80 160 320 640];
                    svnow = [2 4 8 16 32 48 64]; 
                    smoo_text = 'Smoothing {\itm}';
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
            excell = 140;

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

            maps_for_plotting = cell(length(bvnow),length(svnow));
            fname = [fig_dir mname '_panelA_maps.mat'];
            if exist(fname,'file') && ~override_maps
                mem = 1;
                disp(sprintf('\t\t...loading old maps'))
                load(fname,'maps_for_plotting')
            else
                disp(sprintf('\t\t...generating maps'))            
                mem = 0;
            end
            loopout = looper(length(bvnow)*length(svnow));        
            for xx = 1:length(bvnow)
                for yy = 1:length(svnow)
                    ax = axes('Units','pixels','Position',[xvec(xx) yvec(yy) psiz psiz]);    
                        if xx==1 && yy==1
                            ah = add_panel_title('a',sprintf('Effect of bin size and smoothing'),'yoffset',5,'xoffset',-30,'width',300);                    
                        end

                        if mem
                            ratemap = maps_for_plotting{xx,yy};
                        else                    
                            rmset = mapset; % rate mapper settings structure - passing mapset directly to ratemapper causes a memory leak
                            rmset.method = mname;
                            rmset.binsize = bvnow(xx);
                            rmset.ssigma = svnow(yy);
                            rmset.maplims = [-max(abs(epoly(:,1))) -max(abs(epoly(:,2))) max(abs(epoly(:,1))) max(abs(epoly(:,2)))]; 
                            rmset.pot = ppot;
                            rmset.spt = pspt;
                            if strcmp(mname,'fyhn')
                                rmset.maxdist = 320;
                                rmset.mindist = 50;                                
                            end
                            [ratemap,dwellmap,spikemap,rmset,speedlift] = rate_mapper(pos,spk,rmset);
                            maps_for_plotting{xx,yy} = ratemap;
                        end

                        im = imagesc(ratemap,'alphadata',~isnan(ratemap)); hold on;
                        daspect([1 1 1])
                        axis xy off tight
                        colormap(gca,rmap_colormap)     
                        if xx==1 && yy==floor(length(svnow)/2)
                            text(-1.0,-0.5,smoo_text,'HorizontalAl','center','Units','normalized','FontSize',10,'rotation',90)
                        elseif xx==floor(length(bvnow)/2) && yy==length(svnow)
                            text(1.5,-0.48,grid_text,'HorizontalAl','center','Units','normalized','FontSize',10)                    
                        end
                        if yy==length(svnow)
                            text(0.5,-0.14,sprintf('%s',num2str(bvnow(xx))),'HorizontalAl','center','Units','normalized','FontSize',8,'rotation',0)   
                        end
                        if xx==1
                            text(-0.06,0.5,sprintf('%s',num2str(svnow(yy))),'HorizontalAl','right','Units','normalized','FontSize',8,'rotation',0)                             
                        end                    
                        loopout = looper(loopout);        
                end
            end
            if ~mem
                save(fname,'maps_for_plotting');
            end

            axc = axes('Units','pixels','Position',[ax.Position(1)+70 ax.Position(2)+90 12 180]);
                mat = (linspace(0,100,100))';
                imagesc(mat,ones(size(mat)),mat);
                colormap(axc,rmap_colormap);
                axis xy

                axc.YTick = [];
                axc.XTick = [];
                text(0.5,1.25,sprintf('Firing\nrate (Hz)'),'FontSize',8,'HorizontalAl','center','Units','normalized')
                text(0.5,1.1,sprintf('Max'),'FontSize',8,'HorizontalAl','center','Units','normalized')
                text(0.5,-0.1,sprintf('0 Hz'),'FontSize',8,'HorizontalAl','center','Units','normalized')      
        end

%% #################### Heatmap of error by binsize and smoothing
        mapidx = 4;

        xnow = xnow+(psiz*7)+200;
        ynow = ynow-(psiz*7)+20;

        n_color_levels = 64;
        error_colormap = flipud(inferno(n_color_levels));    
        clims = [1*10^-12 1*10^-11];
        %clims = 'auto';    
        interp_method = 'nearest';   

        if 1
            disp(sprintf('\tMISE map...'))

            % main panel showing 16 minute error results
            ax1 = axes('Units','pixels','Position',[xnow ynow-5 300 300]);    
                ah = add_panel_title('b',sprintf('Mean integrated squared error (MISE)'),'yoffset',20,'xoffset',-20,'width',400);
                sigma_now = '4000';
                time_now = '16';              
                load([scratch_space '\' ename '_' mname '_sigma' num2str(sigma_now) '_duration' num2str(time_now) '_datamats.mat'],'data_matrix_error','data_matrix_params')  
                bmat = data_matrix_params(:,:,1); % binsizes used when generating maps
                smat = data_matrix_params(:,:,2); % smoothing used when generating maps
                switch mname
                    case {'histogram','fyhn'}
                        new_bvec = min(bmat(:)):1:max(bmat(:));
                        new_svec = min(smat(:)):1:max(smat(:));  
                    case {'ash'}
                        new_bvec = min(bmat(:)):1:max(bmat(:));
                        new_svec = min(smat(:)):0.1:max(smat(:));  
                    case {'kadaptive'}
                        new_bvec = min(bmat(:)):1:max(bmat(:));
                        new_svec = min(smat(:)):100:max(smat(:)); 
                    case {'kyadaptive'}
                        new_bvec = min(bmat(:)):1:max(bmat(:));
                        new_svec = min(smat(:)):0.1:max(smat(:));                        
                    case {'ksde'}
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
                        text(0,1.1,'16 minutes coverage','HorizontalAl','left','Units','normalized','FontSize',8);                
                    otherwise
                        text(0,1.04,'16 minutes coverage','HorizontalAl','left','Units','normalized','FontSize',8);
                end            
                ylabel(smoo_text);
                xlabel(grid_text);

            axc = axes('Units','pixels','Position',[ax1.Position(1)+ax1.Position(3)+160 ax1.Position(2)+0 12 ax1.Position(4)*0.9]); 
                mat = linspace(ax1.CLim(1),ax1.CLim(2),100)';
                imagesc(ones(size(mat)),mat,mat);
                colormap(axc,error_colormap);
                axis xy

                axc.YTick = linspace(ax1.CLim(1),ax1.CLim(2),4);
                axc.XTick = [];
                axc.YAxisLocation = 'right';
                text(0.5,1.1,sprintf('MISE'),'HorizontalAl','left','Units','normalized','FontSize',8);                

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
            ax2 = axes('Units','pixels','Position',[ax1.Position(1)+ax1.Position(3)+10 ax1.Position(2)+170 130 130]);   
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
                        text(0.01,1.2,'4 minutes','HorizontalAl','left','Units','normalized','FontSize',8);                
                    otherwise
                        text(0.01,1.05,'4 minutes','HorizontalAl','left','Units','normalized','FontSize',8);
                end          
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
            ax3 = axes('Units','pixels','Position',[ax1.Position(1)+ax1.Position(3)+10 ax1.Position(2) 130 130]);    
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
                        text(0.01,1.2,'64 minutes','HorizontalAl','left','Units','normalized','FontSize',8);   
                    case {'ksde'}
                        text(0.01,1.05,'24 minutes','HorizontalAl','left','Units','normalized','FontSize',8);                    
                    otherwise
                        text(0.01,1.05,'64 minutes','HorizontalAl','left','Units','normalized','FontSize',8);
                end 
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

%% #################### Heatmap of time and field error etc
        if 1
            disp(sprintf('\tMapping factors...'))

            xnow = 70;
            ynow = ynow-250;        

            buff = 230;
            vbuff = 180;
            siz = 140;
            posplot = [xnow ynow siz siz; xnow+buff ynow siz siz; xnow ynow-vbuff siz siz; xnow+buff ynow-vbuff siz siz];
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
                    colplot = [0 2000; 0 1000; 0 0.25; 0 1]; % time, memory, missing, field error
                    colname = {'0ms','2000ms';'0kb',sprintf('1{\\times}10^{3}kb');'0%','25%';'0','1'};                    
                otherwise
                    keyboard
            end         
            plotname = {'Computation time (ms) ';'Disk space (kb) ';'Proportion of map empty ';'Place field detection error (fields) '};
            cmapnow = parula(256);        

            time_now = '16'; 
            sigma_now = '4000';
            load([scratch_space '\' ename '_' mname '_sigma' num2str(sigma_now) '_duration' num2str(time_now) '_datamats.mat'],'data_matrix_params','data_matrix_fields','data_matrix_cputime','data_matrix_missing')    

            switch mname
                case {'fyhn'}
                    dat = {data_matrix_cputime(:,:,1:10).*1000,data_matrix_missing(:,:,:,3).*8.*0.001,data_matrix_missing(:,:,:,1),abs( data_matrix_fields(:,:,:) )};
                otherwise
                    dat = {data_matrix_cputime(:,:,:).*1000,data_matrix_missing(:,:,:,3).*8.*0.001,data_matrix_missing(:,:,:,1),abs( data_matrix_fields(:,:,:) )};
            end   
            
            for pp = 1:4
                if pp==1
                    disp(sprintf('\t\t%s ',plotname{pp}))
                else
                    disp(sprintf('\b| %s',plotname{pp}))                
                end

                axp = axes('Units','pixels','Position',posplot(pp,:));  
                    if pp == 1
                        switch mname
                            case {'histogram'}
                                ah = add_panel_title('c',sprintf('Additional mapping factors (16-minute values)'),'yoffset',20,'xoffset',-20,'width',400);                            
                            otherwise
                                ah = add_panel_title('c',sprintf('Additional mapping factors (16-minute values)'),'yoffset',5,'xoffset',-20,'width',400);                            
                        end                     
                    end
                    dme = dat{pp};
                    mat1 = mean(dme,3,'double','omitnan');   
                    bmat = data_matrix_params(:,:,1); % binsizes used when generating maps
                    smat = data_matrix_params(:,:,2); % smoothing used when generating maps                

                    Vq1 = interp2(bmat,smat,mat1,Xq,Yq,interp_method);
                    im = surf(Xq,Yq,Vq1,'EdgeColor','none'); hold on;
                    colormap(gca,cmapnow) 
                    set(gca,'ColorScale','log')

                    axis ij on
                    box off
                    axp.XScale = 'log';
                    axp.YScale = 'log';
                    axp.XLim = ax1.XLim;
                    axp.YLim = ax1.YLim;    
                    axp.YTick = ax1.YTick;
                    axp.XTick = ax1.XTick;             

                    if pp==1 || pp==3
                        ylabel(smoo_text);                        
                    end
                    if pp>2
                        xlabel(grid_text); 
                    end
                    if pp<3
                        axp.XTick = [];                    
                    end

                    view(0,90)
                    caxis(colplot(pp,:))
                    switch mname
                        case {'histogram'}
                            text(0.01,1.18,plotname{pp},'HorizontalAl','left','Units','normalized','FontSize',8);        
                        otherwise
                            text(0.01,1.05,plotname{pp},'HorizontalAl','left','Units','normalized','FontSize',8);        
                    end 

                axc = axes('Units','pixels','Position',[axp.Position(1)+axp.Position(3)+15 axp.Position(2)+10 12 siz*0.8]);
                    mat = linspace(axp.CLim(1),axp.CLim(2),100)';
                    imagesc(mat,ones(size(mat)),mat);
                    colormap(axc,cmapnow);
                    axis xy
                    set(gca,'ColorScale','log')
                    caxis(colplot(pp,:))

                    axc.YTick = [];
                    axc.XTick = [];
                    text(1.1,1,colname{pp,2},'FontSize',8,'HorizontalAl','left','Units','normalized')
                    text(1.1,0,colname{pp,1},'FontSize',8,'HorizontalAl','left','Units','normalized')  

                if strcmp(mname,'histogram')
                    axpb = axes('Units','pixels','Position',[axp.Position(1) axp.Position(2)+axp.Position(4)+4 axp.Position(3) 10]);  
                        imagesc(mat1(1,:)); hold on;
                        colormap(gca,cmapnow) 
                        set(gca,'ColorScale','log')

                        axis ij
                        caxis(axp.CLim);
                        axpb.YTick = 1;
                        axpb.YTickLabel = {'none'};
                        axpb.XTick = [];
                        axis off            
                end                

            end
        end

%% #################### Pareto front
%% plot pareto front of time vs error   
        if 1
            disp(sprintf('\tPareto fronts...'))

            xnow = xnow+520;
            ynow = ynow-160;
            ax = axes('Units','pixels','Position',[xnow ynow 270 300]);  
                ah = add_panel_title('d',sprintf('Pareto-optimisation of binsize & smoothing'),'yoffset',0,'xoffset',-20,'width',400);

                % load data matrices and prepare them for pareto optimisation 
                cols = winter(3);
                shaps = {'.','.','.'};
                mins = NaN(config.npcells,3);
                all_f = cell(1,3);
                all_x = cell(1,3);
                pls = cell(1,3);

                switch mname
                    case {'ksde'}
                        ds = [4 16 24];            
                    otherwise
                        ds = [4 16 64];            
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
                        pls{ii,1} = plot(a,b,'Color',cols(ii,:),'Marker',shaps{ii},'LineStyle','none','MarkerFaceColor',cols(ii,:)); hold on;

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
                        plot(chull(:,1),chull(:,2),'Color',cols(ii,:),'LineWidth',1.5);                    

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
                    pls{ii,1} = plot(a,b,'Color',cols(ii,:),'Marker',shaps{ii},'LineStyle','none','MarkerFaceColor',cols(ii,:)); hold on; 
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
                    plot(chull(:,1),chull(:,2),'Color',cols(ii,:),'LineWidth',1.5);                    
                    
                    plot3(balance_solution1(3),balance_solution1(4),1000,'Color','k','Marker','o','MarkerSize',10,'LineStyle','none');            

                end     
                xlabel(sprintf('MISE'), 'Interpreter', 'none' );
                ylabel(sprintf('Time (ms)'), 'Interpreter', 'none' );
                ax = gca;    
                ax.YLim(1) = 0;
                ax.XLim(1) = 0;
                ax.XScale = 'log';
                ax.YScale = 'log';
                grid off
                box off
                switch mname
                    case {'ksde'}
                        labs_now = {'4 mins Pareto front','16 mins Pareto front','24 mins Pareto front'};
                    otherwise
                        labs_now = {'4 mins Pareto front','16 mins Pareto front','64 mins Pareto front'};
                end              
                [~,leg] = legendflex([pls{1,1} pls{2,1} pls{3,1}] ,labs_now,'anchor',{'e','e'},'ncol',1,'box','off','buffer',[150,-150],'xscale',.5,'fontsize',9); 
                leg(5).MarkerSize = 12;
                leg(7).MarkerSize = 12;
                leg(9).MarkerSize = 12;

        end

%% #################### Solutions on an error map 
        xnow = xnow+290;

        if 1
            disp(sprintf('\t\t...inset plot'))

            ax2 = axes('Units','pixels','Position',[xnow ynow+130 175 175]);    
                sigma_now = '4000';  
                time_now = '16';                        
                load([scratch_space '\' ename '_' mname '_sigma' num2str(sigma_now) '_duration' num2str(time_now) '_datamats.mat'],'data_matrix_error','data_matrix_fields','data_matrix_cputime','data_matrix_params')    
                dme = data_matrix_error(:,:,:,mapidx);
                bmat = double(data_matrix_params(:,:,1)); % binsizes used when generating maps
                smat = double(data_matrix_params(:,:,2)); % smoothing used when generating maps         

                mat1 = mean(dme,3,'double','omitnan');    
                Vq1 = interp2(bmat,smat,mat1,Xq,Yq,interp_method);

                surf(Xq,Yq,Vq1,'EdgeColor','none','FaceAlpha',0.7); hold on;
                colormap(gca,error_colormap) 
                set(gca,'ColorScale','log')
                grid off

                axis ij
                ax2.XScale = 'log';
                ax2.YScale = 'log';
                ax2.XLim = ax1.XLim;
                ax2.YLim = ax1.YLim;     
                ax2.XTick = ax1.XTick;
                ax2.YTick = ax1.YTick;
                view(0,90)
                ylabel(smoo_text);
                xlabel(grid_text);
                caxis(ax1.CLim);

                fval = all_f{2};
                x = all_x{2};
                hold on;

                % find the 'best' solution (closest to minimising everything)
                load([scratch_space '\' ename '_' mname '_sigma' num2str(sigma_now) '_duration' num2str(time_now) '_gamultiobj.mat'],'balance_solution2','balance_solution1','minerr_solution','x')
                pb1 = plot3(balance_solution1(1),balance_solution1(2),1000,'ko','MarkerSize',10,'LineStyle','none');            

                % find the minimal error solution 
                pc = plot3(minerr_solution(1),minerr_solution(2),1000,'kv','MarkerSize',10,'LineStyle','none');          

                [~,leg] = legendflex([pb1 pc]...
                    ,{'Balanced solution','Minimum error solution'}...
                    ,'anchor',{'se','se'},'ncol',1,'box','off','buffer',[-30,-90],'xscale',.5,'fontsize',9);    

        end
        
        g1 = getframe(fig1);
        close(fig1);
    end

%% #################### Prepare data for Field size + duration = bin size + smoothing maps 
    % Figure settings
    fig1 = figure('visible','on','Position',[50,60,1550,900]); 
    set(gcf,'InvertHardCopy','off'); % this stops white lines being plotted black      
    set(gcf,'color','w'); % makes the background colour white
    colormap(jet(256)); % to make sure the colormap is not the horrible default one
    
    ynow = 660;
    xnow = 40;

%% #################### Effect of field size on error
    ynowt = ynow+80;

    siz = 90;
    buff = siz+15;    
    posplot = [xnow ynowt siz siz; xnow ynowt-buff siz siz; xnow ynowt-2*buff siz siz; xnow ynowt-3*buff siz siz];  
    
    if 1
        disp(sprintf('\tExample fields...'))
        
        for i = 1:3
            ax1b = axes('Units','pixels','Position',posplot(i,:));  
                if i==1
                    ah = add_panel_title('e',sprintf('Field size and error'),'yoffset',0,'xoffset',10,'width',400);  
                end            
                sigma_now = config.analysis_pairs{i,1};
                xvec = 0 : 1 : ( max(abs(epoly(:,1)),[],'omitnan') + mapset.padding*10 ); 
                yvec = 0 : 1 : ( max(abs(epoly(:,2)),[],'omitnan') + mapset.padding*10 );                     
                xvec = unique(sort([-xvec xvec],'ascend')); % mirror these vectors and then sort them in ascending order
                yvec = unique(sort([-yvec yvec],'ascend'));     
                [X1,Y1] = meshgrid(movmean(xvec,2,'EndPoints','discard'),movmean(yvec,2,'EndPoints','discard'));                
                Xn = [X1(:) Y1(:)];                
                
                p = mvnpdf(Xn, [0 0],[sigma_now 0;0 sigma_now]);                     
                p = reshape(p,size(X1));

                imagesc(p);
                daspect([1 1 1])
                axis xy off tight
                colormap(gca,rmap_colormap)   
                
                sg = [sigma_now 0;0 sigma_now];                
                e = 2.*sqrt( eig(sg) );
                fsiz = double( geomean(e(:)) );                
                text(0.5,1.1,sprintf('%.1f mm',fsiz),'HorizontalAl','center','Units','normalized','FontSize',10)             
 
            ax1c = axes('Units','pixels','Position',posplot(i,:)+[buff 0 0 0]);  
                sigma_now = config.analysis_pairs{i,1};
                time_now = 16;            
                load([scratch_space '\' ename '_' mname '_sigma' num2str(sigma_now) '_duration' num2str(time_now) '_datamats.mat'],'data_matrix_error','data_matrix_fields','data_matrix_cputime')    
                dme = data_matrix_error(:,:,:,mapidx);
                mat1 = mean(dme,3,'double','omitnan');    
                Vq1 = interp2(bmat,smat,mat1,Xq,Yq,interp_method);

                im = surf(Xq,Yq,Vq1,'EdgeColor','none'); hold on;
                colormap(gca,error_colormap) 
                set(gca,'ColorScale','log')
                switch mname
                    case {'histogram','ash','ksde','kadaptive','kyadaptive','fyhn'}
                        clims = [1*10^-13 1*10^-11]; 
                    otherwise
                        keyboard
                end                 
                caxis(clims)

                axis ij
                ax1c.XScale = 'log';
                ax1c.YScale = 'log';
                ax1c.XLim = [min(new_bvec) max(new_bvec)];
                ax1c.YLim = [min(new_svec) max(new_svec)];     
                ax1c.YTick = [];
                ax1c.XTick = [];
                view(0,90) 
                
            if i==2
                axc = axes('Units','pixels','Position',[ax1c.Position(1)+ax1c.Position(3)+10 ax1c.Position(2)-10 12 ax1c.Position(4)*1.3]); 
                    mat = linspace(ax1c.CLim(1),ax1c.CLim(2),100)';
                    imagesc(ones(size(mat)),mat,mat);
                    colormap(axc,error_colormap);
                    axis xy

                    axc.YTick = linspace(ax1c.CLim(1),ax1c.CLim(2),4);
                    axc.XTick = [];
                    axc.YAxisLocation = 'right';
                    text(0.5,1.2,sprintf('MISE'),'FontSize',8,'HorizontalAl','center','Units','normalized')
            end                  
        end 
    end
        
%% #################### Prepare data for Field size + duration = bin size + smoothing maps 
    xnow = xnow+300;
    if 1
        disp(sprintf('\tPreparing Pareto data...'))
        
        fig2 = figure('visible','off','Position',[50,60,900,900]); 
        set(gcf,'InvertHardCopy','off'); % this stops white lines being plotted black      
        set(gcf,'color','w'); % makes the background colour white
        colormap(jet(256)); % to make sure the colormap is not the horrible default one
            
        % load data matrices and prepare them for pareto optimisation 
        switch mname
            case {'ksde'}
                ds = [4 4 4 16 16 16 24 24 24];
            otherwise
                ds = [4 4 4 16 16 16 64 64 64];
        end          
        fs = [4000 8000 16000 4000 8000 16000 4000 8000 16000];
        for ii = 1:length(ds)
            cols = winter(length(ds));
            
            if pp==1
                disp(sprintf('\t\tradius %d ',fs(ii)))
            else
                disp(sprintf('\b | radius %d',fs(ii)))                
            end
            
            fname = [scratch_space '\' ename '_' mname '_sigma' num2str(fs(ii)) '_duration' num2str(ds(ii)) '_gamultiobj.mat'];
            if exist(fname,'file') && 1 && ~override_pareto
                load([scratch_space '\' ename '_' mname '_sigma' num2str(fs(ii)) '_duration' num2str(ds(ii)) '_gamultiobj.mat'],'x','fval','balance_solution1','minerr_solution','a','b','utopia_point');
                
                set(0,'currentfigure',fig2);                     
                subplot(3,3,ii)
                    plot(a,b,'Color','b','Marker','.','LineStyle','none','MarkerFaceColor','b'); hold on;
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
                    plot(chull(:,1),chull(:,2),'Color',cols(ii,:),'LineWidth',1.5);                      

                    plot(balance_solution1(3),balance_solution1(4),'Color','k','Marker','o','MarkerSize',10,'LineStyle','none'); hold on;
                    plot(utopia_point(1).*max(a(:)),utopia_point(2).*max(b(:)),'Color','b','Marker','x','MarkerSize',20,'LineStyle','none'); hold on;                

                    text(0,1.1,sprintf('Duration %d, field size %d',ds(ii),fs(ii)),'Units','normalized','HorizontalAl','left','FontSize',8)

                    xlabel(sprintf('MISE'), 'Interpreter', 'none' );
                    ylabel(sprintf('Time (ms)'), 'Interpreter', 'none' );
                    ax = gca;
                    ax.YLim(1) = 0;
                    ax.XLim(1) = 0;
                    ax.XScale = 'log';
                    ax.YScale = 'log';
                    grid off
                    box off                 

                ax2 = axes('Units','pixels','Position',[ax.Position(1)+ax.Position(3)/2 ax.Position(2)+ax.Position(4)/2 ax.Position(3)/2 ax.Position(4)/2]);                           
                    Vq1 = interp2(bmat,smat,mat1,Xq,Yq,interp_method);
                    surf(Xq,Yq,Vq1,'EdgeColor','none','FaceAlpha',0.7); hold on;
                    colormap(gca,error_colormap) 
                    set(gca,'ColorScale','log')
                    grid off

                    axis ij
                    ax2.XScale = 'log';
                    ax2.YScale = 'log';
                    ax2.XLim = [min(new_bvec) max(new_bvec)];
                    ax2.YLim = [min(new_svec) max(new_svec)];     
                    ax2.XTick = [];
                    ax2.YTick = [];
                    view(0,90)
                    caxis(clims) 
                    plot(balance_solution1(1),balance_solution1(2),'Color','k','Marker','o','MarkerSize',10,'LineStyle','none'); hold on;
                
                continue
            else
                load([scratch_space '\' ename '_' mname '_sigma' num2str(fs(ii)) '_duration' num2str(ds(ii)) '_datamats.mat'],'data_matrix_params','data_matrix_cputime','data_matrix_fields','data_matrix_error','data_matrix_missing')   
            end
            data_matrix_error = data_matrix_error(:,:,:,mapidx);
            mat1 = mean(data_matrix_error,3,'double','omitnan');    
            bmat = double(data_matrix_params(:,:,1)); % binsizes used when generating maps
            smat = double(data_matrix_params(:,:,2)); % smoothing used when generating maps
            
            cmap1 = double( mean(data_matrix_error,3,'omitnan') );                
            tmap1 = double( mean(data_matrix_cputime,3,'omitnan') ) .* 1000;
            fmap1 = abs( double( mean(data_matrix_fields,3,'omitnan') ) );                
            emap1 = abs( double( mean(data_matrix_missing(:,:,:,1),3,'omitnan') ) );     

            % actual functions that will be used by gamultiobj               
            meth = 'nearest';
            f1 = @(xy) interp2(bmat,smat,cmap1,xy(:,1),xy(:,2),meth,inf); % MSE
            f2 = @(xy) interp2(bmat,smat,tmap1,xy(:,1),xy(:,2),meth,inf); % time
            f3 = @(xy) interp2(bmat,smat,emap1,xy(:,1),xy(:,2),meth,inf); % missing data (NaNs)
            f4 = @(xy) interp2(bmat,smat,fmap1,xy(:,1),xy(:,2),meth,inf); % field count error
            fx = @(xy) [f1(xy) f2(xy) f3(xy) f4(xy)];

            % first find the pareto front when minimising error and time only
            % this result will form the 2D plot in the figure
            % fit curve to objective pareto front
            rng(999);
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

            [v,midx] = min(a,[],1,'omitnan');                
            balance_solution2 = [x(midx,:) fval(midx,:)];                    

            [v,midx] = min(b,[],1,'omitnan');                
            balance_solution3 = [x(midx,:) fval(midx,:)];                  

            i = find(mat1==min(mat1(:),[],'omitnan'),1,'first');
            minerr_solution = [bmat(i),smat(i),cmap1(i),tmap1(i)];                 
            save([scratch_space '\' ename '_' mname '_sigma' num2str(fs(ii)) '_duration' num2str(ds(ii)) '_gamultiobj.mat'],'x','fval','balance_solution1','balance_solution2','balance_solution3','minerr_solution','a','b','utopia_point');
        
            set(0,'currentfigure',fig2);           
            ax = subplot(3,3,ii);
                plot(a,b,'Color','b','Marker','.','LineStyle','none','MarkerFaceColor','b'); hold on;
                idx = convhull(a,b);    
                chull = [a(idx) b(idx)];
                [~,cut_point2] = max(chull(:,2));      
                if cut_point2==1
                    chull = [a(idx(1:end-1)) b(idx(1:end-1))];                        
                    [~,cut_point2] = max(chull(:,2));                              
                end
                chull = [chull(1:cut_point2-1,:) ; NaN(1,2) ; chull(cut_point2:end,:)];
                plot(chull(:,1),chull(:,2),'Color','b','LineWidth',1.5);                    

                plot(balance_solution1(3),balance_solution1(4),'Color','k','Marker','o','MarkerSize',10,'LineStyle','none'); hold on;
                plot(utopia_point(1).*max(a(:)),utopia_point(2).*max(b(:)),'Color','b','Marker','x','MarkerSize',20,'LineStyle','none'); hold on;

                text(0,1.1,sprintf('Duration %d, field size %d',ds(ii),fs(ii)),'Units','normalized','HorizontalAl','left','FontSize',8)
                    
                xlabel(sprintf('MISE'), 'Interpreter', 'none' );
                ylabel(sprintf('Time (ms)'), 'Interpreter', 'none' );
                ax.YLim(1) = 0;
                ax.XLim(1) = 0;
                ax.XScale = 'log';
                ax.YScale = 'log';
                grid off
                box off                 
                
            ax2 = axes('Units','normalized','Position',[ax.Position(1)+ax.Position(3)/2 ax.Position(2)+ax.Position(4)/2 ax.Position(3)/2 ax.Position(4)/2]);    
                Vq1 = interp2(bmat,smat,mat1,Xq,Yq,interp_method);
                surf(Xq,Yq,Vq1,'EdgeColor','none','FaceAlpha',0.7); hold on;
                colormap(gca,error_colormap) 
                set(gca,'ColorScale','log')
                grid off

                axis ij
                ax2.XScale = 'log';
                ax2.YScale = 'log';
                ax2.XLim = [min(new_bvec) max(new_bvec)];
                ax2.YLim = [min(new_svec) max(new_svec)];     
                ax2.XTick = [];
                ax2.YTick = [];
                view(0,90)
                caxis(clims) 
                plot(balance_solution1(1),balance_solution1(2),'Color','k','Marker','o','MarkerSize',10,'LineStyle','none'); hold on;

        end   
        figname = [fig_dir2 '\' mname '_fig_all_pareto_fronts.png'];
        [~,~,~] = mkdir(fig_dir2);    
        exportgraphics(fig2,figname,'BackgroundColor','w','ContentType','image','Resolution',200);  
        close(fig2)
    end
    set(0,'currentfigure',fig1);                            

%% #################### Field size + duration = bin size + smoothing maps
    if 1
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
            
        for rr = 1:2
            switch mname
                case {'histogram','ksde'}
                    if rr==1 
                        xnowt = xnow;
                        disp(sprintf('\tMinimum regression results...'))
                        [b,s,r] = mvregress([ones(size(dat,1),1),dat(:,1:2)],dat(:,6)); 
                        f = @(xy) b(1,1) + xy(:,1).*b(2,1) + xy(:,2).*b(3,1); % [field size, duration] > [bin size, smoothing]
                        f_string = sprintf('f(r,d) = %.4f + %.4fr + %.4fd',b(1,1),b(2,1),b(3,1));
                    else
                        xnowt = xnow+410;
                        disp(sprintf('\tPareto-optimal regression results...'))   
                        [b,s,r] = mvregress([ones(size(dat,1),1),dat(:,1:2)],dat(:,3:4)); 
                        f = @(xy) [b(1,1) + xy(:,1).*b(2,1) + xy(:,2).*b(3,1), b(1,2) + xy(:,1).*b(2,2) + xy(:,2).*b(3,2)]; % [field size, duration] > [bin size, smoothing]
                        f_string = sprintf('f(r,d) = (%.4f + %.4fr + %.4fd, %.4f + %.4fr + %.4fd)',b(1,1),b(2,1),b(3,1),b(1,2),b(2,2),b(3,2));                
                    end                    
                otherwise
                    if rr==1 
                        xnowt = xnow;
                        disp(sprintf('\tMinimum regression results...'))
                        [b,s,r] = mvregress([ones(size(dat,1),1),dat(:,1:2)],dat(:,5:6)); 
                        f = @(xy) [b(1,1) + xy(:,1).*b(2,1) + xy(:,2).*b(3,1), b(1,2) + xy(:,1).*b(2,2) + xy(:,2).*b(3,2)]; % [field size, duration] > [bin size, smoothing]
                        f_string = sprintf('f(r,d) = (%.4f + %.4fr + %.4fd, %.4f + %.4fr + %.4fd)',b(1,1),b(2,1),b(3,1),b(1,2),b(2,2),b(3,2)); 
                    else
                        xnowt = xnow+410;
                        disp(sprintf('\tPareto-optimal regression results...'))   
                        [b,s,r] = mvregress([ones(size(dat,1),1),dat(:,1:2)],dat(:,3:4)); 
                        f = @(xy) [b(1,1) + xy(:,1).*b(2,1) + xy(:,2).*b(3,1), b(1,2) + xy(:,1).*b(2,2) + xy(:,2).*b(3,2)]; % [field size, duration] > [bin size, smoothing]
                        f_string = sprintf('f(r,d) = (%.4f + %.4fr + %.4fd, %.4f + %.4fr + %.4fd)',b(1,1),b(2,1),b(3,1),b(1,2),b(2,2),b(3,2)); 
                    end                    
            end  
% f_string
% keyboard
            colmapnow = magma(256);

            ax1 = axes('Units','pixels','Position',[xnowt ynow-15 160 160]); 
                disp(sprintf('\t\t...bin size plot'))

                if rr==1
                    ah = add_panel_title('f',sprintf('Minimum error solution bin size and smoothing strength'),'yoffset',25,'xoffset',-10,'width',400); 
                else
                    ah = add_panel_title('g',sprintf('Balanced solution bin size and smoothing strength'),'yoffset',25,'xoffset',0,'width',400);                     
                end

                fs = 125 : 1 : 300;
                switch mname
                    case {'ksde'}
                        du = 4 : 1 : 32;
                    otherwise
                        du = 4 : 1 : 64;
                end                 

                [xq,yq] = meshgrid(fs,du);
                fx = f([xq(:),yq(:)]);
                switch mname
                    case {'histogram','ksde'}
                        if rr==1
                            bi = ones(size(xq)); % minimum error bin size is always 1mm
                            sm = reshape(fx(:,1),size(xq)); % only smoothing is included                  
                        else
                            bi = reshape(fx(:,1),size(xq)); % minimum error bin size                 
                            sm = reshape(fx(:,2),size(xq)); % minimum error smoothing
                        end
                    otherwise
                        bi = reshape(fx(:,1),size(xq)); % minimum error bin size                 
                        sm = reshape(fx(:,2),size(xq)); % minimum error smoothing
                end                 


                im = surf(fs,du,bi,'EdgeColor','none'); hold on;        
                colormap(gca,colmapnow) 

                axis ij
                ax1.XScale = 'log';
                ax1.YScale = 'log';
                ax1.XLim = [min(fs) max(fs)];
                ax1.YLim = [min(du) max(du)];     
                ax1.XTick = [125:50:300];
                switch mname
                    case {'ksde'}
                        ax1.YTick = [4 8 16 24 32];
                    otherwise
                        ax1.YTick = [4 8 16 32 64];
                end                 
                view(0,90)
                ylabel(sprintf('Duration (mins)'));
                xlabel(sprintf('Avg. field radius (mm)'));
                text(0,1.17,f_string,'Units','normalized','HorizontalAl','left','FontSize',8)        
                text(0,1.05,grid_text,'Units','normalized','HorizontalAl','left','FontSize',8)        

            axc = axes('Units','pixels','Position',[ax1.Position(1) ax1.Position(2)-75 ax1.Position(3) 10]); 
                mat = linspace(min(bi(:)),max(bi(:)),100);
                imagesc(mat,ones(size(mat)),mat);
                colormap(axc,colmapnow);
                axis xy

                axc.YTick = [];
                axc.XTickLabelRotation = 0;    
                if rr==1
                    switch mname
                        case {'histogram'}
                            axc.XTick = 1; 
                            axc.XLim = [0.5 1.5];
                        case {'ash'}   
                            axc.XTick = 0:5:100;    
                        case {'ksde'}
                            axc.XTick = 0:20:100;   
                        case {'kyadaptive'}
                            axc.XTick = 0:0.25:100;  
                        case {'kadaptive'}
                            axc.XTick = 1:0.25:100; 
                        case {'fyhn'}
                            axc.XTick = 0:1:100;                            
                        otherwise
                            keyboard
                    end                     
                    text(0,1.7,['Minimum error ' grid_text],'FontSize',8,'HorizontalAl','left','Units','normalized') 
                else
                    switch mname
                        case {'histogram'}
                            axc.XTick = 0:10:100;                            
                        case {'ash'}   
                            axc.XTick = 0:10:1000;   
                        case {'ksde'}
                            axc.XTick = 0:5:100;   
                        case {'kyadaptive'}
                            axc.XTick = 0:5:1000;  
                        case {'kadaptive'}
                            axc.XTick = 0:2:1000; 
                        case {'fyhn'}
                            axc.XTick = 0:5:100;                            
                        otherwise
                            keyboard
                    end                    
                    text(0,1.7,['Pareto-optimal ' grid_text],'FontSize',8,'HorizontalAl','left','Units','normalized')
                end
                
            ax2 = axes('Units','pixels','Position',ax1.Position+[ax1.Position(3)+20 0 0 0]);  
                disp(sprintf('\t\t...smoothing plot'))

                im = surf(fs,du,sm,'EdgeColor','none'); hold on;        
                colormap(gca,colmapnow) 

                axis ij
                ax2.XScale = 'log';
                ax2.YScale = 'log';
                ax2.XLim = [min(fs) max(fs)];
                ax2.YLim = [min(du) max(du)];     
                ax2.XTick = [125:50:300];
                ax2.YTick = [4 8 16 32 64];
                view(0,90)
                xlabel(sprintf('Avg. field radius (mm)'));
                text(0,1.05,smoo_text,'Units','normalized','HorizontalAl','left','FontSize',8)        

            axc = axes('Units','pixels','Position',[ax2.Position(1) ax2.Position(2)-75 ax2.Position(3) 10]); 
                mat = linspace(min(sm(:)),max(sm(:)),100);
                imagesc(mat,ones(size(mat)),mat);
                colormap(axc,colmapnow);
                axis xy

                axc.YTick = [];
                axc.XTickLabelRotation = 0;                
                if rr==1
                    switch mname
                        case {'histogram'}
                            axc.XTick = 0:5:100;                            
                        case {'ash'}   
                            axc.XTick = 0:10:100;
                        case {'ksde'}
                            axc.XTick = 0:10:50;  
                        case {'kyadaptive'}
                            axc.XTick = 0:2:100; 
                        case {'kadaptive'}
                            axc.XTick = 0:10000:50000; 
                            axc.XTickLabelRotation = 45;
                        case {'fyhn'}
                            axc.XTick = 0:5:100;                            
                        otherwise
                            keyboard
                    end                     
                    text(0,1.7,['Minimum error ' smoo_text],'FontSize',8,'HorizontalAl','left','Units','normalized')                    
                else
                    switch mname
                        case {'histogram'}
                            axc.XTick = 0:10:100;  
                        case {'ash'}   
                            axc.XTick = 0:0.25:100;
                        case {'ksde'}
                            axc.XTick = 0:5:50;   
                        case {'kyadaptive'}
                            axc.XTick = 0:2:100; 
                        case {'kadaptive'}
                            axc.XTick = 0:2500:50000;  
                            axc.XTickLabelRotation = 45;  
                        case {'fyhn'}
                            axc.XTick = 0:10:100;  
                        otherwise
                            keyboard
                    end                     
                    text(0,1.7,['Pareto-optimal ' smoo_text],'FontSize',8,'HorizontalAl','left','Units','normalized')
                end
                
        end
    end
    g2 = getframe(fig1);
    close(fig1);
    
    imall = [g1.cdata; g2.cdata];
    imall = imall(140:3175,20:2750,:);

%     figure
%     imshow(imall)
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ################################################################# %% Save the figure
    % Save the figure    
    disp(sprintf('\tSaving figure...'))    
    figname1 = [fig_dir2 '\' mname '_fig_summary.png'];
    imwrite(imall,figname1,'png')
        
    
    
    
    
    
  
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    