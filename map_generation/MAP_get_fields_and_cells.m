%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ################################################################# %% DESCRIPTION
%FUNCTION  short desc.
%
% USAGE:
%     out = template(in) process with default settings
% 
%     out = template(in,optional1) process using optional argument 1
% 
%     out = template(___,Name,Value,...) process with Name-Value pairs used to control aspects 
%     of the process
% 
%     Parameters include:
% 
%     'param1'          -   (default = X) Scalar value, parameter to do something
% 
%     'param2'          -   (default = X) Scalar value, parameter to do something
% 
% INPUT:
%     in    - input as a vector
%
% OUTPUT:
%     out   - output as a vector
%
% EXAMPLES:
%
%     % run function using default values
%     out = template(in,varargin)
%
% See also: FUNCTION2 FUNCTION3

% HISTORY:
% version 1.0.0, Release 00/00/00 Initial release
%
% Author: Roddy Grieves
% UCL, 26 Bedford Way
% eMail: r.grieves@ucl.ac.uk
% Copyright 2020 Roddy Grieves

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ################################################################# %% INPUT ARGUMENTS CHECK
%% Prepare default settings

%% ##################################### Heading 2
%% #################### Heading 3
%% Heading 4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ################################################################# %% FUNCTION BODY
    % run through different environments
    for ee = 1:length(config.environs)   
        
        % get position data (random walks)
        nwalks = config.nwalks;
        all_pos = cell(nwalks,1);
        for ww = 1:nwalks
            wlength = 64;
            wdir = [scratch_space '\' config.environs{ee} '_' num2str(wlength) '_walk' num2str(ww) config.append{1} '.mat'];
            load(wdir,'pox','poy','pot','epoly','emat','wlength'); % load walk data       

            pox = fillmissing(pox(:),'linear');
            poy = fillmissing(poy(:),'linear');

            duration = numel(pox)*(1/50);
            all_pos(ww,1) = {[pox(:) poy(:) pot(:)]};
        end        
        
        % run through different field sizes
        for ss = 1:size(config.analysis_pairs,1) % for every field size requested          
            mean_ssnow = config.analysis_pairs{ss,1};

            % generate fields and spikes
            fdir = [scratch_space '\' config.environs{ee} '_sigma' num2str(mean_ssnow) '_duration' num2str(wlength) config.append{1} '_fields.mat'];            

            if ~exist(fdir,'file') || overwrite_fields
                loopout = looper(config.nfields);
                field_dat = table;
                for ff = 1:config.nfields
                    % generate probability map
                    % get working sigmas (covariance)
                    v1 = normal_rand(config.covariance(ss,1),1000,1000,inf,1);
                    v2 = normal_rand(config.covariance(ss,1),1000,1000,inf,1);   

                    while 1
                        v3 = mean([v1 v2]).*normal_rand(0,0.25,-1,1,1);
                        sigma = round([v1 v3;v3 v2]);
                        [~,p] = cholcov(sigma);
                        if ~p
                            break
                        end
                    end

                    % get working mean (centres of mass)
                    mu = [randi(round([min(epoly(:,1)) max(epoly(:,1))]),1); randi(round([min(epoly(:,2)) max(epoly(:,2))]),1)]';                  

                    % accumulate data
                    field_datp = table;     
                    field_datp.centroid(1,:) = single(mu);
                    field_datp.sd(1,:) = single( [v1 v2 v3] );
                    field_dat = [field_dat; field_datp];

                    loopout = looper(loopout);                
                end         
                disp(sprintf('\t\t...saving %s',fdir));            
                save(fdir,'field_dat','-v7.3'); % save field data      
            else
                disp(sprintf('\t\t...loading %s',fdir));            
                load(fdir,'field_dat'); % save field data              
            end

%% #################### Make cells from fields
            sdir = [scratch_space '\' config.environs{ee} '_sigma' num2str(mean_ssnow) '_duration' num2str(wlength) config.append{1} '_sdata.mat'];            

            if ~exist(sdir,'file') || overwrite_cells
                loopout = looper(config.npcells);
                sdata = table;
                
                % get grid of bin centres for probability maps
                xvec = 0 : 1 : ( max(abs(epoly(:,1)),[],'omitnan') + mapset.padding*10 ); 
                yvec = 0 : 1 : ( max(abs(epoly(:,2)),[],'omitnan') + mapset.padding*10 );                     
                xvec = unique(sort([-xvec xvec],'ascend')); % mirror these vectors and then sort them in ascending order
                yvec = unique(sort([-yvec yvec],'ascend'));     
                [X1,Y1] = meshgrid(histcents(xvec),histcents(yvec));                
                X = [X1(:) Y1(:)];
                    
                for pp = 1:config.npcells
                    % get number of fields
                    R = round(gamrnd(config.field_dist(1),config.field_dist(2),100,1));
                    R(R<1) = [];
                    nf = R(1);

                    % choose fields randomly
                    findx = randperm(config.nfields,nf);
                   
                    % get the underlying PDF for this field
                    mus = field_dat.centroid(findx,:);
                    sds = field_dat.sd(findx,:);                       
                    ptot = zeros([size(X1),length(mus(:,1))]);
                    for ff = 1:length(mus(:,1))       
                        p = mvnpdf(X, mus(ff,:), [sds(ff,1) sds(ff,3);sds(ff,3) sds(ff,2)]);                     
                        p = reshape(p,size(X1));
                        %p = imgaussfilt(p,50);
                        ptot(:,:,ff) = p;                    
                    end            
                    pmap = max(ptot,[],3,'omitnan'); % this way gaussians will not overlap and accumulate

                    % get the position data
                    pindx = round(1 + (nwalks-1) .* rand(1,1)); % random index for random walk
                    pos_now = all_pos{pindx,1};
                    pox = pos_now(:,1);
                    poy = pos_now(:,2);
                    pot = pos_now(:,3);
                    
                    % get the firing probability at every position data point
                    coeff1 = 0.5; % weighting for pmap probability, bigger = more spikes
                    denom1 = 256; % denominator for spike lag, bigger = smaller lag, zero disables
                    Fr = ( pmap./max(pmap(:)) ) .* coeff1;                    
                    Vp = interp2(X1,Y1,Fr,pox,poy,'nearest',0); 

                    if config.overdisperse(1) % if we want to include overdispersion
                        % create a 1Hz sinusoid
                        f = 1/(config.overdisperse(2).*2);
                        t = pot;
                        time_modulation = sin(2*pi*f*t);
                        time_modulation(time_modulation>0) = 1;
                        time_modulation(time_modulation<0) = -1;
                        time_modulation = (config.overdisperse(3).*time_modulation)+1;

                        % time_modulation = (config.overdisperse(3).*time_modulation)+1;
                        % Therefore, the values for the amplitude and the average duration of the input modulation 
                        % that produce histograms like what was experimentally observed are around 10% and 1 s.
                        % https://www.sciencedirect.com/science/article/pii/S0306452201005863

                        % Vp = Vp .* time_modulation; % modulate the firing by the 1Hz sinusoid

                        % figure
                        % cline(pox,poy,Vp)
                        % daspect([1 1 1])
                        % keyboard
                        % figure
                        % plot(pot,time_modulation)
                        % keyboard
                    else
                        time_modulation = ones(size(pot));
                    end


                    nout = normal_rand(128,256,0,inf,length(pox)) ./ denom1;
                    idx = 1:length(Vp);
                    idx2 = idx(:) + round(nout);
                    idx2(idx2>max(idx(:),[],'omitnan')) = max(idx(:),[],'omitnan');
                    Vp2 = Vp(idx2);

                    % convert this to spikes
                    ind = 1:length(pox(:));
                    noiselam = 0.001; % lambda for poisson noise, higher = more noise    
                    spk_prob = poissrnd(Vp2) + poissrnd(noiselam,size(Vp2));
                    spk_prob = round(spk_prob .* time_modulation);
                    spk_prob(spk_prob<0) = 0;
                    ind2 = repelem(ind,spk_prob);    
                    spx = pox(ind2,:);                 
                    spy = poy(ind2,:);                 
                    spt = pot(ind2,:);        

% figure
% ax1 = subplot(1,2,1);
% t = 5*60;
% cline(pox(pot<t),poy(pot<t),Vp2(pot<t))
% daspect([1 1 1])
% 
% ax2 = subplot(1,2,2);
% plot(pox(pot<t),poy(pot<t),'k'); hold on;
% plot(spx(spt<t),spy(spt<t),'r.','MarkerSize',20)
% daspect([1 1 1])
% 
% linkaxes([ax1 ax2],'xy')
% 
% figure
% plot(pot,Vp2);
% hold on
% plot(spt,zeros(size(spt)),'ro')
% keyboard

                    % downsample to get desired firing rate
                    num = round( normal_rand(config.frate(1),config.frate(2),config.frate(3),config.frate(4),1) * duration );
                    if numel(spx)>num
                        rindx = randperm(numel(spx),num);
                        spx = spx(rindx);
                        spy = spy(rindx);
                        spt = spt(rindx);
                    end                

%                     figure
%                     subplot(1,3,1)
%                     imagesc(Fr); hold on;
%                     daspect([1 1 1])
%                     axis xy
%                                         
%                     subplot(1,3,2)
%                     plot(pox,poy,'k'); hold on;
%                     plot(spx,spy,'r.','MarkerSize',25);
%                     daspect([1 1 1]);
%                     axis xy
%                     
%                     subplot(1,3,3)
%                     config2 = struct;
%                     config2.dwellmap = [];
%                     config2.speedlift = [];                            
%                     config2.method = 'histogram';
%                     config2.binsize = 1;           
%                     config2.ssigma = 5;
%                     config2.smethod = 1;
%                     config2.ash = svect(ss); 
%                     config2.kern = 'biweight';                                
%                     config2.padding = mapset.padding;
%                     config2.mindwell = mapset.mindwell;
%                     config2.mindist = mapset.mindist;
%                     config2.maxdist = mapset.maxdist;
%                     config2.steps = mapset.steps;                                
%                     config2.srate = mapset.srate;
%                     config2.maplims = NaN(1,4);
%                     config2.maplims([1 3]) = [-max(abs(epoly(:,1))) max(abs(epoly(:,1)))]./10;
%                     config2.maplims([2 4]) = [-max(abs(epoly(:,2))) max(abs(epoly(:,2)))]./10;
%                     delta = config2.binsize;
%                     [rmap,dmap,smap,cfig] = graphDATA([pox poy]./10,[spx spy]./10,[],config2);  
%                     imagesc(rmap);
%                     axis xy
%                     daspect([1 1 1])
%                                         
%                     keyboard                    
                    
                    %% accumulate data
                    sdatap = table;            
                    sdatap.uci = { [datestr(now,'dd.mm.yyyy') '.' num2str(pp)] };
                    sdatap.nspikes = numel(spx);
                    sdatap.frate = numel(spx)/duration;    
                    sdatap.nfields = numel(findx(:));
                    sdatap.centroids = {field_dat.centroid(findx,:)};
                    sdatap.sigs = {field_dat.sd(findx,:)};        
                    sdatap.spk = { single([spx spy spt]) };
                    sdatap.pmap_1mm = { pmap };                      
                    sdatap.pos_index = pindx;                      
                    fnum1 = graphPEAK(pmap,'method','threshold','lowpass',0,'threshold',0.2,'binsize',0.1,'min_area',36);
                    sdatap.detected_fields = fnum1;

                    % get actual field size
                    sigs = sdatap.sigs;
                    es = NaN(size(sigs,1),1);
                    for nn = 1:size(sigs,1)
                        sg = sigs{nn,:};
                        sg = [sg(1,1) sg(1,3);sg(1,3) sg(1,2)];
                        e = 2.*sqrt( eig(sg) );
                        es(nn,1) = double( geomean(e(:)) );          
                    end
                    sdatap.actual_field_size = nanmean(es);   
                    
                    % accumulate
                    sdata = [sdata; sdatap];

                    loopout = looper(loopout);    
                end

                disp(sprintf('\t\t...saving %s',sdir));                        
                save(sdir,'sdata','-v7.3'); % save field data         

            else
                disp(sprintf('\t\t...found %s',sdir));            

            end
        end
    end


































