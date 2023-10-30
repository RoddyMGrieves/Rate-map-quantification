% function MAP_generate_maps
%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> DESCRIPTION
% MAP_generate_maps  script to analyse data for:
% Grieves et al. (2020) Quantification of firing rate map procedures
% This script is dependent on graphtest for inputs and data
%
% USAGE:
%       MAP_generate_maps process with default settings
%
% See also: graphtest MAP_get_random_walks rate_mapper

% HISTORY:
% version 1.0.0, Release 12/07/20 Initial release/comments added (script generated before this date)
%
% Author: Roddy Grieves
% Dartmouth College, Moore Hall
% eMail: roddy.m.grieves@dartmouth.edu
% Copyright 2021 Roddy Grieves

%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Heading 3
%% >>>>>>>>>>>>>>>>>>>> Heading 2
%% >>>>>>>>>> Heading 1
%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> INPUT ARGUMENTS CHECK
%% Parse inputs
    multiWaitbar( 'CloseAll' );
    multiWaitbar( 'Environments', 0 );
    multiWaitbar( 'Method', 0 );
    multiWaitbar( 'Field sizes', 0 );   
    multiWaitbar( 'Session length', 0 );    
    multiWaitbar( 'Cell', 0 );
    multiWaitbar( 'Bin Size', 0 );
    multiWaitbar( 'Smoothing', 0 );

%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> FUNCTION BODY
    for ee = 1:length(config.environs) % for every environment        

        % get position data (random walks)
        nwalks = config.nwalks;
        all_pos = cell(nwalks,1);
        for ww = 1:nwalks
            wlength = 64;
            wdir = [scratch_space '\' config.environs{ee} '_' num2str(wlength) '_walk' num2str(ww) config.append{1} '.mat'];
            load(wdir,'pox','poy','pot','epoly','emat','wlength'); % load walk data, positions are in mm      

            pox = fillmissing(pox(:),'linear');
            poy = fillmissing(poy(:),'linear');

            duration = numel(pox)*(1/50);
            all_pos(ww,1) = {[pox(:) poy(:) pot(:)]};
        end        
        
%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> for every mapping method requested     
        for mm = 1:length(mapset.methods)
            disp(sprintf('Method: %s...',mapset.methods{mm}));            
            analysis_pairs = config.analysis_pairs;
            
            % map settings for this loop
            method_now = mapset.methods{mm};
            bvect = mapset.([method_now '_binsize']);
            svect = mapset.([method_now '_smoothing']); 
 
            field_sizes_now = [config.analysis_pairs{:,1}];
            disp(sprintf('\tAnalysing field sizes: %s',mat2str(field_sizes_now)));      
            
%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> for every field size requested           
            for ff = 1:length(field_sizes_now)             
                mean_ssnow = field_sizes_now(ff);
                disp(sprintf('\tField size: %d...',mean_ssnow));      

                trial_lengths_now = [config.analysis_pairs{ff,2}];
                disp(sprintf('\t\tAnalysing trial lengths: %s',mat2str(trial_lengths_now)));      
                
                % quickly test to see if we should load the field data for analysis
                for ww = 1:length(trial_lengths_now)              
                    data_dir = [scratch_space '\' config.environs{ee} '_' method_now '_sigma' num2str(mean_ssnow) '_duration' num2str(trial_lengths_now(ww)) config.append{1} '_datamats.mat']; 
                    if ~exist(data_dir,'file') || overwrite_maps || append_maps
                        % load cell and field data
                        fdir = [scratch_space '\' config.environs{ee} '_sigma' num2str(mean_ssnow) '_duration' num2str(wlength) config.append{1} '_fields.mat'];
                        disp(sprintf('\t\t\t...loading fields: %s',fdir));            
                        load(fdir,'field_dat');  
                        sdir = [scratch_space '\' config.environs{ee} '_sigma' num2str(mean_ssnow) '_duration' num2str(wlength) config.append{1} '_sdata.mat'];            
                        disp(sprintf('\t\t\t...loading cells: %s',sdir));            
                        load(sdir,'sdata');   
                        break
                    end
                end

%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> for every trial length              
                for ww = 1:length(trial_lengths_now)              
                    disp(sprintf('\t\tTrial length: %d...',trial_lengths_now(ww)));                                

                    % if the map data has already been processed and we don't want to overwrite it skip to the next loop
                    data_dir = [scratch_space '\' config.environs{ee} '_' method_now '_sigma' num2str(mean_ssnow) '_duration' num2str(trial_lengths_now(ww)) config.append{1} '_datamats.mat']; 
                    if exist(data_dir,'file') && ~overwrite_maps && ~append_maps
                        disp(sprintf('\t\t\t...completed already'));                              
                        ss = 1;
                        continue
                    end                      

                    % if we want to append new values to an existing error matrix, rather than completely overwrite it
                    % load dwellmap data if there is any and we don't want to overwrite it
                    slift_name = [config.environs{ee} '_' method_now '_sigma' num2str(mean_ssnow) '_duration' num2str(trial_lengths_now(ww)) config.append{1} '_dwellmaps.mat'];
                    ddir = [];
                    if exist([scratch_space '\' slift_name],'file') % one drive
                        ddir = [scratch_space '\' slift_name];
                    elseif exist(['C:\Users\F004KS7\Desktop\Mapping dwellmaps\' slift_name],'file') % bonsai or desktop computer
                        ddir = ['C:\Users\F004KS7\Desktop\Mapping dwellmaps\' slift_name];
                    elseif exist(['C:\Users\admin\Desktop\Mapping dwellmaps\' slift_name],'file') % cheetah computer
                        ddir = ['C:\Users\admin\Desktop\Mapping dwellmaps\' slift_name];                     
                    end
                    if isempty(ddir) % if the speedlift file is not found, default to onedrive
                        ddir = [scratch_space '\' slift_name];
                    end

                    if append_maps && exist(data_dir,'file')
                        load(data_dir,'data_matrix_cputime','data_matrix_error','data_matrix_fields','data_matrix_missing'); % save field data
                        load(ddir,'all_speedlifts'); % save field data                        
                        bvect_old = mapset.([method_now '_binsize_old']);
                        svect_old = mapset.([method_now '_smoothing_old']);  
                        
                        % preallocate some arrays
                        d1 = NaN(length(svect),length(bvect),size(sdata,1),'single'); % will hold time to process maps
                        d2 = NaN(length(svect),length(bvect),size(sdata,1),4,'single'); % will hold error values
                        d3 = NaN(length(svect),length(bvect),size(sdata,1),'single'); % fields in map minus fields in PDF
                        d4 = NaN(length(svect),length(bvect),size(sdata,1),3,'single'); % missing pixels    
                        dd2 = cell(length(svect),length(bvect),max(sdata.pos_index(:))); % speedlifts                       

                        for bb = 1:length(bvect) 
                            bb_idx = find(bvect_old==bvect(bb),1);     
                            if isempty(bb_idx)
                                continue
                            end
                            for ss = 1:length(svect)  
                                ss_idx = find(svect_old==svect(ss),1);
                                if isempty(ss_idx)
                                    continue
                                end                                
                                d1(ss,bb,:) = data_matrix_cputime(ss_idx,bb_idx,:);
                                d2(ss,bb,:,:) = data_matrix_error(ss_idx,bb_idx,:,:);
                                d3(ss,bb,:,:) = data_matrix_fields(ss_idx,bb_idx,:,:);
                                d4(ss,bb,:) = data_matrix_missing(ss_idx,bb_idx,:);                                
                                dd2(ss,bb,:) = all_speedlifts(ss_idx,bb_idx,:);                                    
                            end
                        end
                        
                        data_matrix_cputime = d1;
                        data_matrix_error = d2;
                        data_matrix_fields = d3;
                        data_matrix_missing = d4;
                        all_speedlifts = dd2;

                    else
                        % just preallocate some empty arrays
                        data_matrix_params = NaN(length(svect),length(bvect),2,'single'); % will hold time to process maps                        
                        data_matrix_cputime = NaN(length(svect),length(bvect),size(sdata,1),'single'); % will hold time to process maps
                        data_matrix_error = NaN(length(svect),length(bvect),size(sdata,1),3,'single'); % will hold error values
                        data_matrix_fields = NaN(length(svect),length(bvect),size(sdata,1),'single'); % fields in map minus fields in PDF
                        data_matrix_missing = NaN(length(svect),length(bvect),size(sdata,1),2,'single'); % missing data
                        
                        if exist(ddir,'file') && ~overwrite_maps % if we have some dwellmap information saved
                            load(ddir,'all_speedlifts'); % save field data
                        else
                            all_speedlifts = cell(length(svect),length(bvect),max(sdata.pos_index(:)));                         
                        end                             
                    end                     

%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> for every cell              
                    for uu = 1:size(sdata,1)                       
                        savd = 0; % whether to save dwellmaps or not

                        % get this cell's settings and data
                        mus = double( sdata.centroids{uu} );
                        sds = double( sdata.sigs{uu} );                       
                        spk = sdata.spk{uu};

                        % cut the position and spike data to the trial length
                        % we want
                        pindx = sdata.pos_index(uu);
                        pos_now = all_pos{pindx,1};
                        pox = pos_now(:,1); % in mm
                        poy = pos_now(:,2); % in mm
                        pot = pos_now(:,3); % in seconds                 

                        max_t = trial_lengths_now(ww);
                        
                        tindx = pot<(max_t.*60);
                        ppox = pox(tindx);
                        ppoy = poy(tindx);
                        ppot = pot(tindx);
                        tindx = spk(:,3)<(max_t.*60);
                        pspx = spk(tindx,1);
                        pspy = spk(tindx,2);
                        pspt = spk(tindx,3);
                        
                        %% Get the PDF values
                        switch mapset.methods{mm}
                            case {'ash'}
                                % get grid of bin centres for probability maps
                                xvec = 0 : (min(bvect)/max(svect)) : ( max(abs(epoly(:,1)),[],'omitnan') + mapset.padding ); 
                                yvec = 0 : (min(bvect)/max(svect)) : ( max(abs(epoly(:,2)),[],'omitnan') + mapset.padding );                     
                                xvec = unique(sort([-xvec xvec],'ascend')); % mirror these vectors and then sort them in ascending order
                                yvec = unique(sort([-yvec yvec],'ascend'));     
                                [X1,Y1] = meshgrid(movmean(xvec,2,'EndPoints','discard'),movmean(yvec,2,'EndPoints','discard')); 
                                X = [X1(:) Y1(:)];                                
                                
                                % get the underlying PDF for this field
                                mus = sdata.centroids{uu};
                                sds = sdata.sigs{uu};                       
                                ptot = zeros([size(X1),length(mus(:,1))]);
                                for munow = 1:length(mus(:,1))       
                                    p = mvnpdf(X, mus(munow,:), [sds(munow,1) sds(munow,3);sds(munow,3) sds(munow,2)]);                     
                                    p = reshape(p,size(X1));
                                    ptot(:,:,munow) = p;                    
                                end            
                                ptot = max(ptot,[],3,'omitnan'); % this way gaussians will not overlap and accumulate                                                               
                                
                            otherwise
                                ptot = sdata.pmap_1mm{uu};
                                % pdf coordinate values
                                xvec = 0 : 1 : ( max(abs(epoly(:,1)),[],'omitnan') + mapset.padding ); 
                                yvec = 0 : 1 : ( max(abs(epoly(:,2)),[],'omitnan') + mapset.padding );                     
                                xvec = unique(sort([-xvec xvec],'ascend')); % mirror these vectors and then sort them in ascending order
                                yvec = unique(sort([-yvec yvec],'ascend'));     
                                [X1,Y1] = meshgrid(movmean(xvec,2,'EndPoints','discard'),movmean(yvec,2,'EndPoints','discard')); 
                                Xpdf = [X1(:) Y1(:)];
                        end
                        p = sdata.pmap_1mm{uu};
                        m1 = p ./ sum(p(:),'omitnan') ./ 1;     

%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> for every binsize                           
                        for bb = 1:length(bvect)                        
                            % sort out where our real bin edges will be
                            bnow = bvect(bb);
                            pmap_now = NaN;

%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> for every smoothing value
                            if strcmp(method_now,'fyhn') && uu>10 % with the fyhn method it is faster to produce all smoothed variants at once
                                % if we are appending values rather than overwriting, skip this value if it already exists
                                if append_maps
                                    if ~isnan(data_matrix_cputime(1,bb,uu))
                                        continue
                                    end
                                end

                                %% generate the firing rate map
                                maplims = NaN(1,4);
                                maplims([1 3]) = [-max(abs(epoly(:,1))) max(abs(epoly(:,1)))];
                                maplims([2 4]) = [-max(abs(epoly(:,2))) max(abs(epoly(:,2)))];                                   

                                pos = [ppox ppoy]; % positions in mm
                                spk = [pspx pspy]; % spikes in mm
                                
                                rmset = mapset; % rate mapper settings structure - passing mapset directly to ratemapper causes a memory leak
                                rmset.method = mapset.methods{mm};
                                rmset.binsize = bnow;
                                rmset.ssigma = svect;
                                rmset.maplims = maplims;  
                                rmset.pot = ppot;
                                rmset.spt = pspt;
                                rmset.maxdist = 320;
                                rmset.mindist = 50;
                                
                                %% >>>>>>>>>>>>>>>>> Ratemapper code
                                [rmaps,~,~,rmset,~] = rate_mapper(pos,spk,rmset,[]);
                            end

                            for ss = 1:length(svect)  
                                if strcmp(method_now,'fyhn') && uu>10 
                                    ratemap = rmaps(:,:,ss);
                                else
                                    % if we are appending values rather than overwriting, skip this value if it already exists
                                    if append_maps
                                        if ~isnan(data_matrix_cputime(ss,bb,uu))
                                            continue
                                        end
                                    end

                                    %% generate the firing rate map and dwell map
                                    maplims = NaN(1,4);
                                    maplims([1 3]) = [-max(abs(epoly(:,1))) max(abs(epoly(:,1)))];
                                    maplims([2 4]) = [-max(abs(epoly(:,2))) max(abs(epoly(:,2)))];                                   
                                    if ~isempty(all_speedlifts{ss,bb,sdata.pos_index(uu)})
                                        speedlift = all_speedlifts{ss,bb,sdata.pos_index(uu)};
                                    else
                                        speedlift = [];
                                    end   

                                    pos = [ppox ppoy]; % positions in mm
                                    spk = [pspx pspy]; % spikes in mm

                                    rmset = mapset; % rate mapper settings structure - passing mapset directly to ratemapper causes a memory leak
                                    rmset.method = mapset.methods{mm};
                                    rmset.binsize = bnow;
                                    rmset.ssigma = svect(ss);
                                    rmset.ash = svect(ss);
                                    rmset.maplims = maplims;  
                                    rmset.pot = ppot;
                                    rmset.spt = pspt;
                                    if strcmp(rmset.method,'fyhn')
                                        rmset.maxdist = 320;
                                        rmset.mindist = 50;
                                    end
                                
                                    %%----%% For inspection only
                                        if 0
                                            rmset.binsize = 20;
                                            rmset.ssigma = 1000;
                                            rmset.ash = 10;                                        
                                            rmset.method = 'ash';
                                            dwellmap = [];
                                            speedlift = [];
                                            disp(sprintf('Method %s, Binsize %d, smoothing %.1f, filter size %.2f',rmset.method,rmset.binsize,rmset.ssigma,2*ceil(2*(rmset.ssigma./rmset.binsize))+1))
                                            %keyboard
                                        end
                                    %%----%%

                                    %% >>>>>>>>>>>>>>>>> Ratemapper code
                                    tic;
                                    [ratemap,dwellmap,spikemap,rmset,speedlift] = rate_mapper(pos,spk,rmset,speedlift);
                                    t1 = toc;
                                end
                                
                                if strcmp(rmset.method,'ash')
                                    delta = rmset.binsize/rmset.ash;
                                    d = delta.^2;                                                                

                                    [X2,Y2] = meshgrid(movmean(rmset.xgrid_ash,2,'EndPoints','discard'),movmean(rmset.ygrid_ash,2,'EndPoints','discard')); 
                                    
                                    pmap_now = interp2(X1,Y1,ptot,X2,Y2,'nearest',NaN);
                                    p2 = pmap_now ./ sum(pmap_now(:),'omitnan') ./ d; % normalize to unit integral (pdf)                                     
                                else
                                    delta = rmset.binsize;
                                    d = delta.^2;     
                                    
                                    [X2,Y2] = meshgrid(movmean(rmset.xgrid,2,'EndPoints','discard'),movmean(rmset.ygrid,2,'EndPoints','discard')); 

                                    if isnan(pmap_now)
                                        pmap_now = interp2(X1,Y1,ptot,X2,Y2,'nearest',NaN);
                                        p2 = pmap_now ./ sum(pmap_now(:),'omitnan') ./ d; % normalize to unit integral (pdf) 
                                    end                                    
                                end
                                if isempty(all_speedlifts{ss,bb,sdata.pos_index(uu)})
                                    all_speedlifts{ss,bb,sdata.pos_index(uu)} = single(speedlift);
                                    savd = 1;
                                end          

                                % normalize ratemap
                                r2 = ratemap ./ sum(ratemap(:) , 'omitnan') ./ d; % normalize to unit integral (pdf)
                                
                                % calculate error estimates                                
                                % data_matrix_error(ss,bb,uu,1) = mean( (p2(:)-r2(:)).^2 , 'omitnan' ) .* d; % mean integrated squared error                                
                                % data_matrix_error(ss,bb,uu,3) = sum( sum( (p2 .* log(p2./r2)) , 'omitnan' ) , 'omitnan' ) .* d; % KL divergence      

                                m2 = imresize(ratemap,size(m1),'nearest'); % resize to match 1mm probability map
                                m2 = m2 ./ sum(m2(:),'omitnan') ./ 1; % normalize    
                                data_matrix_error(ss,bb,uu,4) = mean( (m1(:)-m2(:)).^2 , 'omitnan' ) .* 1; % calculate MISE
                                data_matrix_error(ss,bb,uu,3) = corr(m1(:),m2(:),'rows','pairwise','type','Pearson'); % correlation
                                % data_matrix_error(ss,bb,uu,1) = mi(double(m1),double(m2)); % mutual information

% m1 = m1 ./ max(m1(:),[],'omitnan');
% m2 = m2 ./ max(m2(:),[],'omitnan');
% rmap8 = im2uint8(m1);
% dmap8 = im2uint8(m2);
% p = imhist(rmap8(:));
% q = imhist(dmap8(:));
% p = p ./ numel(rmap8);
% q = q ./ numel(dmap8);
% p = p + eps; % add smallest possible amount to avoid log of 0
% q = q + eps;
% data_matrix_error(ss,bb,uu,1) = MutualInformation(rmap8(:),dmap8(:)); % mutual information
% data_matrix_error(ss,bb,uu,1) = sum(p .* (log2(p) - log2(q))); % KLD

                                % euclidean distance
                                m1b = m1(:);
                                m2b = m2(:);
                                idx = ~isnan(m1b) & ~isnan(m2b);
                                data_matrix_error(ss,bb,uu,2) = pdist([m1b(idx)';m2b(idx)'],'euclidean');     

                                %%----%% For inspection only                                
                                    if 0
                                        figure
                                        subplot(2,2,1)
                                        imagesc(p2);
                                        daspect([1 1 1])
                                        title('pdf')
                                        axis xy

                                        subplot(2,2,2)
                                        imagesc(r2,'alphadata',~isnan(r2))
                                        daspect([1 1 1])
                                        title('rmap')
                                        axis xy
                                        e = mean( (p2(:)-r2(:)).^2 , 'omitnan' ) .* d;
                                        title(sprintf('error = %.10f',e))

                                        subplot(2,2,3)
                                        plot(pos(:,1),pos(:,2),'k'); hold on;
                                        plot(spk(:,1),spk(:,2),'r.','MarkerSize',25);
                                        daspect([1 1 1]);
                                        axis xy
                                        return   
                                    end
                                %%----%%
                                
                                %% accumulate data
                                % time taken to generate the map
                                data_matrix_params(ss,bb,1) = bnow;
                                data_matrix_params(ss,bb,2) = svect(ss);                               
                                data_matrix_cputime(ss,bb,uu) = t1;                                    
                                data_matrix_missing(ss,bb,uu,1) = sum( isnan(r2(:)) ) ./ numel(r2(:)); % total missing pixels                                                                                                      
                                data_matrix_missing(ss,bb,uu,2) = sum( r2(:)==0 ) ./ numel(r2(:)); % total missing pixels                                                                                                      
                                data_matrix_missing(ss,bb,uu,3) = numel(r2); % total pixels                                                                                                      

                                % place field count
                                %fnum1 = graphPEAK(sdata.pmap_1mm{uu},'method','threshold','lowpass',0,'threshold',0.2,'binsize',0.1,'min_area',36);
                                fnum1 = sdata.detected_fields(uu);
                                fnum2 = graphPEAK(ratemap,'method','threshold','lowpass',0,'threshold',0.2,'binsize',delta./10,'min_area',36);
                                data_matrix_fields(ss,bb,uu) = fnum1-fnum2;  
                                multiWaitbar( 'Smoothing', ss/length(svect) );                                 
                            end
                            multiWaitbar( 'Bin Size', bb/length(bvect) );                        
                        end
                        
                        %%----%% For inspection only                                
                            if 1
                                figure
                                % subplot(2,3,1)
                                x = data_matrix_error(:,:,1,1);
                                im = imagesc(x);
                                set(im,'alphadata',~isnan(x))
                                daspect([1 1 1])
                                axis ij
                                title('MSE')
                                caxis('auto')
                                colormap(flipud(inferno(64)))
                                set(gca,'ColorScale','log') 
                                colorbar
                                keyboard
                            end
                        %%----%%

                        if savd % if we added some information to our dwellmap arrays
                            save(ddir,'all_speedlifts','-v7.3'); % save field data
                        end
                        multiWaitbar( 'Cell', uu/size(sdata,1) );
                    end
                    % save the map output data
                    save(data_dir,'data_matrix_params','data_matrix_cputime','data_matrix_error','data_matrix_fields','data_matrix_missing','-v7.3'); % save field data
                    multiWaitbar( 'Session length', ww/length(trial_lengths_now) ); 
                    disp(sprintf('\t\t...done'));                                
                end
                multiWaitbar( 'Field sizes', ff/length(field_sizes_now) );                
            end
            multiWaitbar( 'Method', mm/length(mapset.methods) );
        end
        multiWaitbar( 'Environments', ee/length(config.environs) );              
    end


































