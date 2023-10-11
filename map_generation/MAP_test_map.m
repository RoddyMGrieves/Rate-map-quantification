




% close all;
% np = 1000; % number of position points to simulate
% pos = rand(np,2)*1000; % simulate position data
% spk = normrnd(500,100,ceil(np/10),2); % simulate spike data
% [ratemap,dwellmap,spikemap,rmset] = rate_mapper(pos,spk); % generate map
% figure
% subplot(2,2,1)
% imagesc(ratemap,'alphadata',~isnan(ratemap)); hold on; % plot ratemap
% daspect([1 1 1]); axis xy; title('Ratemap');
% subplot(2,2,2)
% imagesc(dwellmap,'alphadata',~isnan(dwellmap)); hold on; % plot dwellmap
% daspect([1 1 1]); axis xy; title('Dwellmap');
% subplot(2,2,3)
% imagesc(spikemap,'alphadata',~isnan(spikemap)); hold on; % plot spikemap
% daspect([1 1 1]); axis xy; title('Spikemap');
% subplot(2,2,4)
% imagesc(ratemap,'alphadata',~isnan(ratemap)); hold on; % plot ratemap
% mapXY = rmset.points_to_map(pos); % convert position data to map coordinates
% plot(mapXY(:,1),mapXY(:,2),'w'); % plot map-scaled position data
% daspect([1 1 1]); axis xy; title('Ratemap + positions');
% 
% 
% return
% 

        bvnow = 30;
        svnow = 2;  
        tcut = 32;
        excell = 3;

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

        tcut = tcut*60;
        ppox = pox(pot<tcut);
        ppoy = poy(pot<tcut);
        pspx = spx(spt<tcut);
        pspy = spy(spt<tcut);
        pos = [ppox ppoy]; % positions in mm
        spk = [pspx pspy]; % spikes in mm
                    
        
        rmset = struct;
        rmset.method = 'kyadaptive';
        rmset.binsize = 10;
        rmset.ssigma = 1;
        rmset.maplims = [-max(abs(epoly(:,1))) -max(abs(epoly(:,2))) max(abs(epoly(:,1))) max(abs(epoly(:,2)))];      
        tic
        [ratemap0,dwellmap0,spikemap0,rmset0,speedlift0] = rate_mapper(pos,spk,rmset);
        toc
        
        rmset = struct;
        rmset.method = 'kyadaptive';
        rmset.binsize = 10;
        rmset.ssigma = 1;
        rmset.maplims = [-max(abs(epoly(:,1))) -max(abs(epoly(:,2))) max(abs(epoly(:,1))) max(abs(epoly(:,2)))];      
        tic
        [ratemap,dwellmap,spikemap,rmset,speedlift] = rate_mapper(pos,spk,rmset);
        ratemap = imgaussfilt(ratemap,2.5);
        toc

        figure
        subplot(2,2,1)
        imagesc(ratemap0);
        daspect([1 1 1])
        axis xy on tight
        colormap(gca,flipud(viridis))         

        subplot(2,2,2)        
        imagesc(ratemap);
        daspect([1 1 1])
        axis xy on tight
        colormap(gca,flipud(viridis))  
        
        % subplot(2,2,3)        
        % emat = (1-rmset.coeffv);
        % imagesc(emat);
        % daspect([1 1 1])
        % axis xy on tight
        % colorbar
        % 
        % subplot(2,2,4)        
        % imagesc(ratemap,'alphadata',emat);
        % daspect([1 1 1])
        % axis xy on tight
        % colormap(gca,flipud(viridis))         
        
        return
        
        
        
        
% [r,d1,d2] = ripleyk(pos,spk,10);
% 
% figure
% plot(r.radius_mm,r.h_derivative)
% sdata.actual_field_size(excell)        
%              
% return        
        
        
        
        
        
        
        
        
        rmset = struct;
        rmset.method = 'ksde';
        rmset.binsize = 20;
        rmset.ssigma = 100;
        rmset.ash = 4;        
        rmset.maxdist = 640;
        rmset.maplims = [-max(abs(epoly(:,1))) -max(abs(epoly(:,2))) max(abs(epoly(:,1))) max(abs(epoly(:,2)))];      
        tic
        [ratemap1,dwellmap1,spikemap1,rmset1,speedlift1] = rate_mapper(pos,spk,rmset);
        %[sum(spikemap1(:)) numel(pspx)]
        toc

        figure
        subplot(1,3,1)
        imagesc(sdata.pmap_1mm{excell});
        daspect([1 1 1])
        axis xy on tight
        colormap(gca,flipud(viridis))          
        
        subplot(1,3,2)
        im = imagesc(ratemap1); hold on;
        daspect([1 1 1])
        axis xy on tight
        colormap(gca,flipud(viridis))  

        rmset.method = 'ksde';  
        rmset.binsize = 20;
        rmset.ssigma = 200;        
        tic
        [ratemap2,dwellmap2,spikemap2,rmset2,speedlift2] = rate_mapper(pos,spk,rmset);
        %[sum(spikemap1(:)) numel(pspx)]        
        toc          

        subplot(1,3,3)        
        im = imagesc(ratemap2); hold on;
        daspect([1 1 1])
        axis xy on tight
        colormap(gca,flipud(viridis))  





















