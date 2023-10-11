%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> DESCRIPTION
% DIS_rev2_fig_1  
% Script to analyse data and produce Fig. 1 for:
% Grieves, Shinder, Rosow, Kenna and Taube (2022) 
% The neural correlates of spatial disorientation in head direction cells. eNeuro
%
% USAGE:
%       MAP_generate_maps process with default settings
%
% See also: graphtest MAP_get_random_walks rate_mapper

% HISTORY:
% version 1.0.0, Release 22/11/21 Initial release
%
% Author: Roddy Grieves
% Dartmouth College, Moore Hall
% eMail: roddy.m.grieves@dartmouth.edu
% Copyright 2021 Roddy Grieves

%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Heading 3
%% >>>>>>>>>>>>>>>>>>>> Heading 2
%% >>>>>>>>>> Heading 1
%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> INPUT ARGUMENTS CHECK
    if 0
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
            wlength = 64;
            sdir = [scratch_space '\' config.environs{1} '_sigma8000_duration64_sdata.mat'];            
            disp(sprintf('\t\t...loading %s',sdir));            
            load(sdir,'sdata'); 
        end
    end
    ename = 'arena120cm';
    rmap_colormap = 'turbo';
    
%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> FUNCTION BODY
    % create a figure
    figure('Position',[100 50 1400 800]);
    set(gcf,'InvertHardCopy','off'); % this stops white lines being plotted black      
    set(gcf,'color','w'); % makes the background colour white

%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Protocol Exp 1
    xnow = 150;
    ynow = 550;

    excell = 31;

    wdir = [scratch_space '\' ename '_64_walk1.mat'];
    load(wdir,'epoly');  

    % environment
    axe = axes('Units','pixels','Position',[xnow ynow 180 180]); 
        ah = add_panel_title('a',sprintf('Simulated environment'),'yoffset',5,'xoffset',0,'width',400);  

        plot(epoly(:,1),epoly(:,2),'k','LineWidth',1.5);
        axis off xy
        daspect([1 1 1])

        annotation('arrow',[0.096 0.096],[0.91 0.688]) % y-arrow
        annotation('arrow',[0.108 0.238],[0.67 0.67]) % x-arrow        

        text(0.5,-0.11,'1200 mm','Units','normalized','FontSize',8,'HorizontalAl','center')
        text(-0.12,0.5,'1200 mm','Units','normalized','FontSize',8,'HorizontalAl','center','Rotation',90)

    % rat scale inset
    h = 200;
    scale = (axe.Position(4)/1200)*h;
    axp = axes('Units','pixels','Position',[axe.Position(1)+10 axe.Position(2)+10 scale*2 scale]); 
        ir = imread('C:\Users\F004KS7\OneDrive - Dartmouth College\Projects in prep\2019 Mapping project\associated media\rat_silh.png');
        image(ir);
        axis off ij
        daspect([1 1 1])

    % true PDF for the cell at 1mm resolution
    xnow = xnow+220;
    axp = axes('Units','pixels','Position',[xnow ynow 180 180]);  
        ah = add_panel_title('b',sprintf('Place cell firing probability map'),'yoffset',5,'xoffset',0,'width',400);  

        pmap = sdata.pmap_1mm{excell,1};
        im = imagesc(pmap); hold on;
        daspect([1 1 1])
        axis xy off tight
        colormap(gca,rmap_colormap)  

        cents = sdata.centroids{excell};
        sigma =  sdata.sigs{excell};            
        text(0.1,-0.1,sprintf('Centroid = %d, %d\n%c_{x} = %.f, %c_{y} = %.f, %c_{xy} = %.f',cents(1),cents(2),963,sqrt(sigma(1)),963,sqrt(sigma(2)),963,sigma(3)),'Units','normalized','FontSize',8,'HorizontalAl','left')




    ax = axes('Units','pixels','Position',[axp.Position(1)+200 axp.Position(2)+40 10 120]);
        x = 1:100;
        imagesc(x');
        colormap(gca,rmap_colormap)  
        axis xy off
        text(0.5,-0.1,'0','Units','normalized','FontSize',9,'HorizontalAl','center');
        text(0.5,1.1,'Max','Units','normalized','FontSize',9,'HorizontalAl','center');

    % spikes and position data
    xnow = xnow+250;
    axs = axes('Units','pixels','Position',[xnow ynow 180 180]);   
        ah = add_panel_title('c',sprintf('Random walk & simulated spikes'),'yoffset',5,'xoffset',0,'width',400);  

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
        daspect([1 1 1])
        axis xy off  
        plot(epoly(:,1),epoly(:,2),'k','LineWidth',1.5);
        axis off xy
        daspect([1 1 1])
        text(0,1.05,'16 minutes','Units','normalized','FontSize',8,'HorizontalAl','left');

        [~,leg] = legendflex([p1 s1]...
            ,{'trajectory','spikes'}...
            ,'anchor',{'s','s'},'ncol',1,'box','off','buffer',[-30,-40],'xscale',1,'fontsize',9); 

    % 3 field sizes
    xnow = 150;
    ynow = ynow-260;
    hbuff = 30;
    wid = 180;
    posplots = [xnow ynow wid wid];
    excells = [3 1 21];
    sigmas = [4000 8000 16000];
    
    fs2 = NaN(size(sigmas));    
    for i = 1:length(sigmas)
        sg = [sigmas(i) 0;0 sigmas(i)];                
        e = 2.*sqrt( eig(sg) );
        fs2(1,i) = double( geomean(e(:)) ); % field size in mm
    end  
            
    for i = 1:3
        axn = axes('Units','pixels','Position',[xnow+((i-1)*hbuff)+((i-1)*wid) ynow wid wid]);   
            if i==1
                ah = add_panel_title('d',sprintf('Simulated field sizes'),'yoffset',-5,'xoffset',0,'width',400);                  
            end
            sdir = [scratch_space '\' config.environs{1} '_sigma' num2str(sigmas(i)) '_duration64_sdata.mat'];            
            load(sdir,'sdata');         

            pmap = sdata.pmap_1mm{excells(i),1};
            im = imagesc(pmap); hold on;
            daspect([1 1 1])
            axis xy off tight
            colormap(gca,rmap_colormap)  

            cents = sdata.centroids{excells(i)};
            sigma =  sdata.sigs{excells(i)};            
            text(0.1,-0.1,sprintf('Centroid = %d, %d\n%c_{x} = %.f, %c_{y} = %.f, %c_{xy} = %.f',cents(1),cents(2),963,sqrt(sigma(1)),963,sqrt(sigma(2)),963,sigma(3)),'Units','normalized','FontSize',8,'HorizontalAl','left')
            text(0,1.05,sprintf('Radius %d',i),'Units','normalized','FontSize',8,'HorizontalAl','left');            
            text(0.1,0.95,sprintf('Approx. 2%c radius: %d mm',963,round(fs2(i))),'Units','normalized','FontSize',8,'HorizontalAl','left','Color','w')

    end

    % 3 durations
    xnow = 150;
    ynow = ynow-270;
    hbuff = 30;
    wid = 180;
    posplots = [xnow ynow wid wid];
    durs = [4 16 64];
    exwalk = 2;
    for i = 1:3
        axn = axes('Units','pixels','Position',[xnow+((i-1)*hbuff)+((i-1)*wid) ynow wid wid]);   
            if i==1
                ah = add_panel_title('e',sprintf('Simulated trial durations'),'yoffset',0,'xoffset',0,'width',400);                  
            end        
            pos_now = all_pos{exwalk,1};
            pox = pos_now(:,1);
            poy = pos_now(:,2);
            pot = pos_now(:,3);
            pox = pox(pot<durs(i)*60);
            poy = poy(pot<durs(i)*60);

            p1 = plot(pox,poy,'Color',[.5 .5 .5 .5]); hold on;
            daspect([1 1 1])
            axis xy off  
            plot(epoly(:,1),epoly(:,2),'k','LineWidth',1.5);
            text(0,1.05,sprintf('Duration %d: %d minutes',i,durs(i)),'Units','normalized','FontSize',8,'HorizontalAl','left');

    end

%% >>>>>>>>>> Save figure so far (the final figure is too large to save in one piece)
    if print_now
        disp(sprintf('\tSaving figure...'))    
        figname = [fig_dir2 '\Fig 1 methods.png'];
        [~,~,~] = mkdir(fig_dir2);    
        exportgraphics(gcf,figname,'BackgroundColor','w','ContentType','image','Resolution',250);  
        close(gcf)        
    end         









    
    