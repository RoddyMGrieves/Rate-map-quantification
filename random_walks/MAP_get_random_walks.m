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
            wdir = [scratch_space '\' config.environs{pp} '_' num2str(wlength) '_walk' num2str(ww) config.append{1} '.mat'];
            
            if ~exist(wdir,'file') || overwrite_walks
                fname = [fig_dir '\' config.environs{pp} '_' num2str(wlength) '_plot_params' num2str(ww) config.append{1} '.png'];

                % prepare environment matrix for walk function
                emat2 = imresize(emat,2);       
                b = bwdist(~emat2);
                emat2(b<10) = false;   
                if config.biased_walk==1 % biased to a point
                    push_point = zeros(size(emat2),'logical');
                    push_point(120,60) = true;
                    pull_point = zeros(size(emat2),'logical');
                    pull_point(120,180) = true;                    
                    [pox,poy,pot,wmap] = getWALK2('env_map',emat2,'walk_length',wlength,'step_dist',[64 128],'central_drag',512,'wall_drag',512,'jit',100,'smoo',128,'plot_params_dir',fname,'pull_point',pull_point,'push_point',push_point);
                elseif config.biased_walk==2 % thigmotaxis
                    [pox,poy,pot,wmap] = getWALK2('env_map',emat2,'walk_length',wlength,'step_dist',[64 128],'central_drag',0,'wall_drag',100,'jit',100,'smoo',128,'plot_params_dir',fname);
                else
                    [pox,poy,pot,wmap] = getWALK2('env_map',emat2,'walk_length',wlength,'step_dist',[64 128],'central_drag',512,'wall_drag',512,'jit',100,'smoo',128,'plot_params_dir',fname);
                end
                pox = (pox ./ 2); % position data are in mm
                poy = (poy ./ 2);  

                disp(sprintf('\t\t...saving %s',wdir));
                save(wdir,'pox','poy','pot','wmap','epoly','emat','wlength','-v7.3'); % save walk data     

                if 0
                    figure
                    subplot(3,3,1)
                    plot(epoly(:,1),epoly(:,2),'r'); hold on;
                    plot(pox(:),poy(:),'k');
                    daspect([1 1 1])

                    subplot(3,3,4)
                    tindx = pot<8*60;
                    plot(epoly(:,1),epoly(:,2),'r'); hold on;                    
                    plot(pox(tindx),poy(tindx),'k');
                    daspect([1 1 1])

                    subplot(3,3,5)
                    tindx = pot<16*60;
                    plot(epoly(:,1),epoly(:,2),'r'); hold on;                    
                    plot(pox(tindx),poy(tindx),'k');
                    daspect([1 1 1])             

                    subplot(3,3,6)
                    tindx = pot<64*60;
                    plot(epoly(:,1),epoly(:,2),'r'); hold on;                    
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

    
    
    
    
    
    