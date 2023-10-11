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

%% ##################################### Heading 2
%% #################### Heading 3
%% Heading 4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ################################################################# %% FUNCTION BODY

%% #################### Example cells
for ll = 1:size(config.analysis_pairs,1)
    wlength = 64;
    sdir = [scratch_space '\' config.environs{1} '_sigma' num2str(config.analysis_pairs{ll,1}) '_duration64_sdata.mat'];            
    disp(sprintf('\t\t...loading %s',sdir));            
    load(sdir,'sdata'); 
    
    ncol = 15;
    nrow = 10;
    ncells = nrow*ncol;    
    
    fig1 = figure('visible','off','Position',[50,60,1300,900]); 
    set(gcf,'InvertHardCopy','off'); % this stops white lines being plotted black      
    set(gcf,'color','w'); % makes the background colour white
    colormap(jet(256)); % to make sure the colormap is not the horrible default one
    annotation('textbox',[0.01 1 1 0],'String',sprintf('nFields: %d - Field size: %d',ncells,config.analysis_pairs{ll,1}),'EdgeColor','none')

    rindx = randperm(ncells,ncells);
    
    xnow = 50;
    ynow = 700;
    xwidth = 70;
    ywidth = 70;
    xbuff = 10;
    ybuff = 10;
    
    [x,y] = ndgrid(xnow:(xwidth+xbuff):ncol*(xwidth+xbuff),ynow:-(ywidth+ybuff):ynow-(nrow-1)*(ywidth+ybuff));
    
    for rr = 1:ncells
        excell = rindx(rr);

        pindx = sdata.pos_index(excell);
        pos_now = all_pos{pindx,1};
        pox = pos_now(:,1);
        poy = pos_now(:,2);
        pot = pos_now(:,3); 
                        
        spk = sdata.spk{excell};
        spx = spk(:,1);
        spy = spk(:,2);
        spt = spk(:,3);

        ax = axes('Units','pixels','Position',[x(rr) y(rr) xwidth ywidth]); 
            plot(pox,poy,'Color',[.5 .5 .5 .5]); hold on;
            plot(spx,spy,'r.','MarkerSize',4);
            daspect([1 1 1])
            axis xy off   
            text(0,1,sprintf('Cell %d',excell),'Units','normalized');
    end
    
    % Save the figure    
    figname = [fig_dir '\Cells_fieldsize_' num2str(config.analysis_pairs{ll,1}) '.png'];
    exportgraphics(gcf,figname,'BackgroundColor','w','ContentType','image','Resolution',200);
end
    
    
    
    
    

    
    
    
    
    
    
    
    
    
    