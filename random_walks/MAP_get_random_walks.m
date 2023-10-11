%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ################################################################# %% DESCRIPTION
% Script to analyse data for:
% Grieves et al. (2020) 
% Disorientation in HD cells
%
% This script 
% NOTE: This script is dependent on disoritest to load and prepare data and indices
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
        pp = 1;
        wlength = 64;       
%         [epoly,emat] = contENV(config.environs{pp},0);
%         epoly = epoly .* 10;            
%         epoly = epoly-mean( [ max(epoly,[],1);min(epoly,[],1) ] );    
        arena_size_mm = 1200; % we want a 1.2m or 1200mm square arena
        epoly = [0 0; 0 arena_size_mm/2; 0 arena_size_mm; arena_size_mm arena_size_mm; arena_size_mm arena_size_mm/2; arena_size_mm 0; 0 0]; % maze polygon in mm
        epoly = epoly-mean( [ max(epoly,[],1);min(epoly,[],1) ] );            
        emat = padarray(ones(arena_size_mm./10-2,'logical'),[1 1],false,'both'); % maze matrix map, each pixel is 1cm
            
        nwalks = config.nwalks;
        for ww = 1:nwalks
            wdir = [scratch_space '\' config.environs{pp} '_' num2str(wlength) '_walk' num2str(ww) '.mat'];
            
            if ~exist(wdir,'file') || overwrite_walks
                fname = [fig_dir '\' config.environs{pp} '_' num2str(wlength) '_plot_params' num2str(ww) '.png'];

                % prepare environment matrix for walk function
                emat2 = imresize(emat,2);       
                b = bwdist(~emat2);
                emat2(b<10) = false;      
                [pox,poy,pot,wmap] = getWALK2('env_map',emat2,'walk_length',wlength,'step_dist',[64 128],'central_drag',1024,'wall_drag',1024,'jit',100,'smoo',128,'plot_params_dir',fname);
                pox = (pox ./ 2); % position data are in mm
                poy = (poy ./ 2);  

                disp(sprintf('\t\t...saving %s',wdir));
                save(wdir,'pox','poy','pot','wmap','epoly','emat','wlength','-v7.3'); % save walk data     

                if 0
                    figure
                    plot(epoly(:,1),epoly(:,2),'r'); hold on;
                    plot(pox(:),poy(:),'k');
                    daspect([1 1 1])

                    figure
                    subplot(1,4,1)
                    tindx = pot<8*60;
                    plot(pox(tindx),poy(tindx),'k');
                    daspect([1 1 1])

                    subplot(1,4,2)
                    tindx = pot<16*60;
                    plot(pox(tindx),poy(tindx),'k');
                    daspect([1 1 1])            

                    subplot(1,4,3)
                    tindx = pot<32*60;
                    plot(pox(tindx),poy(tindx),'k');
                    daspect([1 1 1])   

                    subplot(1,4,4)
                    tindx = pot<64*60;
                    plot(pox(tindx),poy(tindx),'k');
                    daspect([1 1 1])               
                    keyboard

                    figure
                    tindx = pot>32*60 & pot<64*60;
                    plot(pox(tindx),poy(tindx),'k');
                    daspect([1 1 1]) 
                    keyboard
                end
            else
                disp(sprintf('\t\t...found %s',wdir));
            end
        end   

%     wdir = [scratch_space '\' config.environs{pp} '_' num2str(wlength) '_walk2.mat'];
%     save(wdir,'pox','poy','pot','wmap','epoly','emat','wlength','-v7.3'); % save walk data     

    
    
    
    
    
    