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
    if ~exist('sdata','var') || 0
        sdir = [scratch_space '\' ename '_sigma16000_duration64_sdata.mat'];            
        disp(sprintf('\t\t...loading %s',sdir));            
        load(sdir,'sdata'); 
    end

    mname = 'histogram';
    
    fs = [4000 8000 16000];
    fs2 = NaN(size(fs));    
    for i = 1:length(fs)
        sg = [fs(i) 0;0 fs(i)];                
        e = 2.*sqrt( eig(sg) );
        fs2(1,i) = double( geomean(e(:)) ); % field size in mm
    end
    cols = copper(max(fs2)+1);    
    
%% ##################################### Heading 2
%% #################### Heading 3
%% Heading 4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ################################################################# %% FUNCTION BODY
    % Figure settings
    fig1 = figure('visible','on','Position',[50,60,1400,900]); 
    set(gcf,'InvertHardCopy','off'); % this stops white lines being plotted black      
    set(gcf,'color','w'); % makes the background colour white
    colormap(jet(256)); % to make sure the colormap is not the horrible default one

%% #################### Spikes and trajectory
    xnow = 150;
    ynow = 500;
    excell = 54;
    ax1 = axes('Units','pixels','Position',[xnow ynow 150 150]);    
        ah = add_panel_title('a',sprintf('Ripleys K demonstration'),'yoffset',10,'xoffset',0,'width',400);  

        pindx = sdata.pos_index(excell);
        pos_now = all_pos{pindx,1};
        pox = pos_now(:,1);
        poy = pos_now(:,2);
        pot = pos_now(:,3);
        spk = sdata.spk{excell};
        spx = spk(:,1);
        spy = spk(:,2);
        spt = spk(:,3);
        pox = pox(pot<16*60);
        poy = poy(pot<16*60);
        spx = spx(spt<16*60);
        spy = spy(spt<16*60);

        p1 = plot(pox,poy,'Color',[.5 .5 .5 .5]); hold on;
        s1 = scatter(spx,spy,12,'r','filled','MarkerEdgeColor','none','MarkerFaceAlpha',0.75);
        plot(epoly(:,1),epoly(:,2),'k','LineWidth',1.5);
        
        daspect([1 1 1])
        axis xy off  
        ax1.XLim = [min(epoly(:,1)) max(epoly(:,1))];
        ax1.YLim = [min(epoly(:,2)) max(epoly(:,2))];  
        text_offset = 1.07;
        text(0,text_offset,'16 minute session','Units','normalized','FontSize',8,'HorizontalAl','left');
    
        centr = [-165 -178];
        rads = 20;
        ds = NaN(rads,3);
        spx2 = pox(randperm(length(pox),numel(spx)));
        spy2 = poy(randperm(length(poy),numel(spy)));
        rng(999); % for reproducibility
        for i = 1:rads
            r = i*60;
            rectangle('Position',[centr(1)-r centr(2)-r r*2 r*2],'Curvature',[1 1],'Clipping','on','EdgeColor','b');
            d(i,1) = r;    
            ds = pdist2(centr,double([spx(:) spy(:)]));
            d(i,2) = sum( ds<r );
            ds = pdist2(centr,double([spx2(:) spy2(:)]));
            d(i,3) = sum( ds<r );            
        end
        re = plot([-10000 -11000],[0 0],'bo');

        [~,leg] = legendflex([p1 s1 re]...
            ,{'positions','spikes','test radii'}...
            ,'anchor',{'s','s'},'ncol',1,'box','off','buffer',[-30,-60],'xscale',1,'fontsize',9);         
        
    axd = axes('Units','pixels','Position',[ax1.Position(1)+ax1.Position(3)+70 ax1.Position(2) ax1.Position(3) ax1.Position(4)]);    
        p1 = plot(d(:,1),d(:,2),'b','Marker','o'); hold on;
        p2 = plot(d(:,1),d(:,3),'k'); hold on;
        
        axd.YLim = [0 1200];
        axd.XLim = [0 1200];
        xlabel('Radius (mm)')
        ylabel('Total spikes')
        box off
        
        [~,leg] = legendflex([p1 p2]...
            ,{'Cell','Uniform random'}...
            ,'anchor',{'se','se'},'ncol',1,'box','off','buffer',[0,-80],'xscale',1,'fontsize',9);        
        
    axk = axes('Units','pixels','Position',[axd.Position(1)+axd.Position(3)+70 axd.Position(2) axd.Position(3) axd.Position(4)]);    
        dnew = d(:,2)-d(:,3);
        plot(d(:,1),dnew,'b'); hold on;      
        x = find(dnew==max(dnew));
        
        axk.YLim = axd.YLim;        
        axk.XLim = axd.XLim;        
        xlabel('Radius (mm)')
        ylabel('Cell - Random')
        box off
        li = line([d(x,1) d(x,1)],axk.YLim,'Color',[.5 .5 .5]);
        text(0.29,0.95,sprintf('Estimated radius: %.1f mm\n(actual = %.1f)',d(x,1),fs2(3)),'Units','norm','Hori','left','FontSize',8)
             
%% #################### Ripley K plot
    xnow = 150;
    ynow = ynow-290;
    excell = 54;
    ax1 = axes('Units','pixels','Position',[xnow ynow 220 150]);     
        ah = add_panel_title('b',sprintf('Estimating field size'),'yoffset',0,'xoffset',0,'width',400);  
    
        fname = [scratch_space '\' ename '_histogram_field_sizes.mat'];
        if exist(fname,'file') && 1
            load(fname,'dat')   
        else
            field_sizes_now = [config.analysis_pairs{:,1}];    
            dat = NaN(size(sdata,1),2,3);
            for ff = 1:length(field_sizes_now) 
                mean_ssnow = field_sizes_now(ff);
                sdir = [scratch_space '\' ename '_sigma' num2str(mean_ssnow) '_duration64_sdata.mat'];            
                load(sdir,'sdata');  
                dat(:,1,ff) = sdata.actual_field_size(:,1);            
                dat(:,2,ff) = sdata.estimated_field_size(:,1);
            end  
            save(fname,'dat','-v7.3')           
        end        
        
        s1 = scatter(dat(:,1,1),dat(:,2,1),10,cols(round(fs2(1)),:),'filled','Marker','o'); hold on;
        s2 = scatter(dat(:,1,2),dat(:,2,2),10,cols(round(fs2(2)),:),'filled','Marker','o');     
        s3 = scatter(dat(:,1,3),dat(:,2,3),10,cols(round(fs2(3)),:),'filled','Marker','o');     
    
        ax1.XLim = [0 320];
        ax1.YLim = [0 320];
        xlabel('Actual field radius (mm)')
        ylabel('Estimated field radius (mm)')
        line([0 ax1.XLim],[0 ax1.YLim],'Color',[.5 .5 .5],'LineStyle','--')
    
        x = reshape(dat(:,1,:),[],1);
        y = reshape(dat(:,2,:),[],1);
        p = polyfit(x, y, 1);
        r = refline(p(1),p(2));
        r.Color = 'k';
        r.LineStyle = '-';

        text(1.04,0.91,sprintf('r = %.1f',corr(x,y,'rows','pairwise','type','Pearson')),'Units','normalized','FontSize',8)
        text(1.04,1.05,sprintf('x = y'),'Units','normalized','FontSize',8)
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ################################################################# %% Save the figure
    % Save the figure    
    disp(sprintf('\tSaving figure...'))    
    figname1 = [fig_dir2 '\ripley_k.png'];
    print_qual = '-r250';
    export_fig(figname1,'-png',print_qual)  % [top,right,bottom,left]
    close(gcf);
        
    
    
    
    
    
  
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    