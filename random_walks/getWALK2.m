function [pox,poy,pot,wmap] = getWALK2(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ################################################################# %% DESCRIPTION
%FUNCTION  short descr.
% long descr.
%
% USAGE:
%         [out] = template(in,in2)
%
% INPUT:
%         in - input 1
%         in2 - input 2
%
% OUTPUT:
%    p - output
%
% EXAMPLES:
%
% See also: FUNCTION2 FUNCTION3

% HISTORY:
% version 1.0.0, Release 00/00/00 Initial release
%
% Author: Roddy Grieves
% UCL, 26 Bedford Way
% eMail: r.grieves@ucl.ac.uk
% Copyright 2019 Roddy Grieves

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ################################################################# %% INPUT ARGUMENTS CHECK
%% Prepare default settings
    def_config              = struct;
    def_env_size            = 128;    
    def_env_map             = padarray(ones(def_env_size-2,'logical'),[1 1],false,'both'); 
    def_p_map               = padarray(zeros(def_env_size-2,'logical'),[1 1],false,'both');     
    def_walk_length         = 20;  
    def_steps_per_min       = 60;  
    def_step_dist           = [24 12]; % [mean SD] (cm) distribution of preferred distances relative to previous position
    def_step_s_fudge        = [0 0]; % [min max] (cm) a random value between these will be added to step_dist SD at each step
    def_step_fudge          = [-8 8]; % [min max] (cm) random distance added to each step will fall between these values
    def_angle_dist          = [0 1]; % [mean kappa] von mises distribution of preferred angles relative to previous angle, lower kappa = flatter distribution
    def_angle_k_fudge       = [0 0]; % [min max] (cm) a random value between these will be added to angle_dist kappa at each step
    def_angle_fudge         = [-2 2]; % [min max] (rads) random angle added at each step will fall between these two values
    def_smoo                = 100; % smoothing factor for smoothn
    def_jit                 = 20; % (pixels) amount by which to jitter position data before smoothing - for a more realistic look
    def_central_drag        = 100; % (SD) strength of gaussian drag to centre of arena, smaller means stronger
    def_wall_drag           = 100; % (SD) strength of gaussian drag away from walls, smaller means stronger
    def_push                = 50; % (SD) strength of gaussian drag away from bias points, smaller means stronger
    def_pull                = 25; % (SD) strength of gaussian drag towards bias points, smaller means stronger
    
%     def_corner_drag         = 100; % (SD) strength of gaussian drag away from walls, smaller means stronger    
    def_plot_params         = 1; % set to 1 to plot parameter curves
    def_plot_params_dir     = pwd; % directory to save parameter figure in
    def_plot_walk           = 0; % set to 1 to plot walk as it is calculated
    def_vid_walk            = 0; % set to 1 to plot walk and save as a video 
    def_fig_vis             = 'off'; % set to 1 to show figures
    def_cname               = 'egosim'; % set to 1 to show figures
    def_date                = datestr(now,30);
    def_pname                = 'a';

%% Parse inputs
    p = inputParser;
    addParameter(p,'config',def_config,@(x) isstruct(x));   
    addParameter(p,'env_size',def_env_size,@(x) isnumeric(x) && isscalar(x));   
    addParameter(p,'env_map',def_env_map,@(x) islogical(x));  
    addParameter(p,'push_point',def_p_map,@(x) islogical(x));  
    addParameter(p,'pull_point',def_p_map,@(x) islogical(x));  
    addParameter(p,'push',def_push,@(x) isnumeric(x) && isscalar(x));
    addParameter(p,'pull',def_pull,@(x) isnumeric(x) && isscalar(x));
    addParameter(p,'walk_length',def_walk_length,@(x) isnumeric(x));
    addParameter(p,'steps_per_min',def_steps_per_min,@(x) isnumeric(x));    
    addParameter(p,'step_dist',def_step_dist,@(x) isnumeric(x) && length(x)==2);
    addParameter(p,'step_s_fudge',def_step_s_fudge,@(x) isnumeric(x) && length(x)==2);
    addParameter(p,'step_fudge',def_step_fudge,@(x) isnumeric(x) && length(x)==2);    
    addParameter(p,'angle_dist',def_angle_dist,@(x) isnumeric(x) && length(x)==2);
    addParameter(p,'angle_k_fudge',def_angle_k_fudge,@(x) isnumeric(x) && length(x)==2);
    addParameter(p,'angle_fudge',def_angle_fudge,@(x) isnumeric(x) && length(x)==2);
    addParameter(p,'smoo',def_smoo,@(x) isnumeric(x) && isscalar(x));
    addParameter(p,'jit',def_jit,@(x) isnumeric(x) && isscalar(x));
    addParameter(p,'central_drag',def_central_drag,@(x) isnumeric(x) && isscalar(x));
    addParameter(p,'wall_drag',def_wall_drag,@(x) isnumeric(x) && isscalar(x));    
    addParameter(p,'plot_params',def_plot_params,@(x) isnumeric(x) && isscalar(x));
    addParameter(p,'plot_params_dir',def_plot_params_dir,@(x) isstring(x) | ischar(x));
    addParameter(p,'plot_walk',def_plot_walk,@(x) isnumeric(x) && isscalar(x));
    addParameter(p,'vid_walk',def_vid_walk,@(x) isnumeric(x) && isscalar(x));    
    addParameter(p,'fig_vis',def_fig_vis,@(x) isstring(x) | ischar(x)); 
    addParameter(p,'cname',def_cname,@(x) ~isempty(x) && ~all(isnan(x(:))) && ischar(x)); 
    addParameter(p,'date',def_date,@(x) isstring(x) || ischar(x)); 
    addParameter(p,'pname',def_pname,@(x) isstring(x) || ischar(x));        
    parse(p,varargin{:});

%% Retrieve parameters 
    % if a config structure (like the one output by this function or graphDATA) was passed in as the 2nd argument
    % use those settings. Otherwise use the default settings. This is useful if you want to manually
    % specify the settings elsewhere in a structure or if you want to pass a structure to the function 
    % to ensure it uses exactly the same settings
    if ~isempty(fieldnames(p.Results.config))
        config = p.Results.config;
        if ~isfield(config,'env_map') || isempty(config.env_map) || all(isnan(config.env_map(:)))
            config.env_map = p.Results.env_map;
        end
    else
        config = p.Results;
    end
    steps = config.walk_length*config.steps_per_min;
    config.steps = steps;
    
%% Display info
    disp(sprintf('\t\t\t\tGenerating random walk...'))
    disp(sprintf('\t\t\t\t...%.1f mins, %d steps',config.walk_length,config.steps))
    
%% ##################################### Heading 2
%% #################### Heading 3
%% Heading 4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ################################################################# %% FUNCTION BODY
%% Plot the input parameters
%     if config.plot_params
%         figure('visible','off')
%         subplot(2,2,1)
%         x = 0:1:40;
%         y = normpdf(x,config.step_dist(1),config.step_dist(2));
%         plot(x,y); hold on
%         y2 = normpdf(x,config.step_dist(1),config.step_dist(2)+config.step_s_fudge(1));
%         y3 = normpdf(x,config.step_dist(1),config.step_dist(2)+config.step_s_fudge(2));
%         plot(x,y2,'r');
%         plot(x,y3,'g');
%         title('Step distance')
% 
%         subplot(2,2,2)
%         x = 0:0.1:2*pi;
%         y = circ_vmpdf(x,config.angle_dist(1),config.angle_dist(2));
%         plot(x,y); hold on
%         y2 = circ_vmpdf(x,config.angle_dist(1),config.angle_dist(2)+config.angle_k_fudge(1));
%         y3 = circ_vmpdf(x,config.angle_dist(1),config.angle_dist(2)+config.angle_k_fudge(2));
%         plot(x,y2,'r');
%         plot(x,y3,'g');
%         title('Angle distance')
% 
%         print('random_walk_parameter_summary','-dpng');
%         close(gcf)
%     end

%% ##################################### Calculate walk    
    % make the map of the environment
    wmap = config.env_map;

    % ready the start point and angle
    lpoint = round(size(wmap)./2);
    langle = -pi + (2*pi*rand(1,1));

    % ready central drag matrix (to gently pull the agent towards the centre)
    if config.central_drag==0
        centre_drag = zeros(size(wmap));
    else
        wmap2 = zeros(size(wmap),'logical');
        wmap2(lpoint(1),lpoint(2)) = true;
        ds = bwdistgeodesic(wmap,wmap2,'quasi-euclidean');
        centre_drag = normpdf(ds,0,config.central_drag); 
        centre_drag = centre_drag./nanmax(centre_drag(:));
    end

    % ready wall push matrix (to gently push the agent off the walls)
    if config.wall_drag==0
        wall_drag = zeros(size(wmap));
    else
        ds = bwdist(~wmap,'quasi-euclidean');
        wall_drag = normpdf(ds,0,config.wall_drag); 
        wall_drag = wall_drag./nanmax(wall_drag(:)); 
    end

    % prepare push/pull matrices (for biasing movements)
    if any(config.push_point(:))
        ds = bwdistgeodesic(wmap,config.push_point,'quasi-euclidean');
        push_mat = normpdf(ds,0,config.push); 
        push_mat = 1-push_mat;
    else
        push_mat = zeros(size(wmap));
    end
    if any(config.pull_point(:))
        ds = bwdistgeodesic(wmap,config.pull_point,'quasi-euclidean');
        pull_mat = normpdf(ds,0,config.pull); 
        pull_mat = pull_mat./nanmax(pull_mat(:));   
    else
        pull_mat = zeros(size(wmap));
    end

    if config.plot_walk || config.vid_walk
        figure('visible',config.fig_vis,'position',[100 100 900 600]);
        set(gcf,'InvertHardCopy','off');   
        set(gcf,'color','w'); % makes the background colour white
        
        if config.vid_walk
            vname = [pwd '\klustest\' config.cname '\' config.cname '_pathvid_' datestr(now,30) '.mat'];
            writerObj = VideoWriter(vname,'MPEG-4');
            writerObj.FrameRate = 2;
            writerObj.Quality = 100;
            open(writerObj);
        end
    end

%% #################### run through every step
    pos = NaN(steps,3);
    loopout = looper(steps);
    sample_rate = 50;
    occ_history_steps = config.steps_per_min * 8;
    for ss = 0:steps
        cpoint = zeros(size(wmap),'logical');
        cpoint(lpoint(1),lpoint(2)) = true;

        % get distance from point to all pixels
        ds1 = bwdistgeodesic(wmap,cpoint,'quasi-euclidean');

        % convert distances to a probability centered on the mean step length
        fudge_factor = config.step_fudge(1)+(config.step_fudge(2)-config.step_fudge(1)).*rand(1,1);
        fudge_factor_s = config.step_s_fudge(1)+(config.step_s_fudge(2)-config.step_s_fudge(1)).*rand(1,1);    
        ds = normpdf(ds1,config.step_dist(1)+fudge_factor,config.step_dist(2)+fudge_factor_s); 
        ds = ds./nanmax(ds(:));    
        idx = ds<.1;
        ds(idx) = NaN;

        % calculate the angle of each pixel in the mask to the current point
        [cc,rr] = meshgrid(1:size(wmap,1),1:size(wmap,2));
        as = atan2(rr-lpoint(1),cc-lpoint(2));
        as(~wmap) = NaN;
        as(idx) = NaN;

        % convert angles to show the difference to the last angle (we want the rat to continue on a similar path if possible)
        as2 = angdiff(as(~idx),ones(sum(~idx(:)),1).*langle);
        fudge_factor = config.angle_fudge(1)+(config.angle_fudge(2)-config.angle_fudge(1)).*rand(1,1);
        fudge_factor_k = config.angle_k_fudge(1)+(config.angle_k_fudge(2)-config.angle_k_fudge(1)).*rand(1,1);    
        as2 = circ_vmpdf(as2, config.angle_dist(1)+fudge_factor, config.angle_dist(2)+fudge_factor_k);
        as3 = NaN(size(wmap));
        as3(~idx) = as2;
        as3 = reshape(as3,size(wmap));
        as3 = as3./nanmax(as3(:));

        % calculate occupancy so far
        idx = pos(:,1)>(ss-occ_history_steps) & pos(:,1)<ss;
        occ = hist3(pos(idx,2:3),{1:size(wmap,1) 1:size(wmap,2)});
        occ = imboxfilt(occ,[17 17]);
        occ = 1-(occ./nanmax(occ(:)));    
        occ = occ+eps;

        % make a decision matrix (by taking the product of the angle and distance matrices and finding the maximum)
        as3 = as3./nanmax(as3(:));
        ds = ds./nanmax(ds(:));
        occ = occ./nanmax(occ(:));
        % dm = as3 .* ds .* occ .* centre_drag .* push_mat .* pull_mat;
        dm = (1.*as3) + (1.*ds) + (1.*occ) + (1.*centre_drag) + (1.*wall_drag) + (1.*push_mat) + (0.9.*pull_mat);

if 0
figure
subplot(3,3,1)
imagesc(as3)
daspect([1 1 1])
title('angle')

subplot(3,3,2)
imagesc(ds)
daspect([1 1 1])
title('distance')

subplot(3,3,3)
imagesc(occ)
daspect([1 1 1])
title('occupancy')

subplot(3,3,4)
imagesc(centre_drag)
daspect([1 1 1])
title('centre drag')

subplot(3,3,5)
imagesc(wall_drag)
daspect([1 1 1])
title('wall drag')

subplot(3,3,6)
imagesc(push_mat)
daspect([1 1 1])
title('push')

subplot(3,3,7)
imagesc(pull_mat)
daspect([1 1 1])
title('pull')


subplot(3,3,8)
imagesc(dm)
daspect([1 1 1])
title('result')

keyboard
end
        lpoint0 = lpoint;
        if all(isnan(dm(:)))
            rn = lpoint(1);
            cn = lpoint(2);
        else
            [rn,cn] = find(dm==max(dm(:)),1,'first');
            lpointn = [rn,cn]; 
            langlen = as(rn,cn);
            if ~isempty(lpointn)
                lpoint = lpointn;
            end
            if ~isempty(langlen)
                langle = langlen;
            end
        end

        if ~ss
            continue
        end
        
        % find the best path between the current point and the next point (avoiding obstacles)
        cpoint = zeros(size(wmap),'logical');
        cpoint(rn,cn) = true;
        ds2 = bwdistgeodesic(wmap,cpoint,'quasi-euclidean');        
        ds3 = round(ds1+ds2);
        ds3(isnan(ds3)) = inf;
        paths = imregionalmin(ds3);
        paths_thinned_many = bwmorph(paths,'thin',inf);        

        ds4 = ds1;
        ds4(~paths_thinned_many) = NaN;
        [I,J,V] = find(ds4);
        nindx = isnan(V);
        dat = [I(~nindx) J(~nindx) V(~nindx)];
        [~,sidx] = sort(dat(:,3));
        dat = dat(sidx,:);

        if isempty(dat) || size(dat,1)==1
            ipath = repmat([rn cn],50,1);
        else
            ipath = [interp1(1:size(dat,1),dat(:,1),linspace(1,size(dat,1),sample_rate+1))' interp1(1:size(dat,1),dat(:,2),linspace(1,size(dat,1),sample_rate+1))'];
            ipath = ipath(2:end,:);
        end

        % collect data
        pos((ss*sample_rate)+1:((ss+1)*sample_rate),:) = [ones(sample_rate,1).*ss ipath(:,1) ipath(:,2)];
        loopout = looper(loopout);        
        
% figure
% subplot(2,2,1)
% imagesc(ds1)
%         
% subplot(2,2,2)
% imagesc(ds2)   
% 
% subplot(2,2,3)
% imagesc(ds3)  
%         
% subplot(2,2,4)
% imshow(~wmap | paths_thinned_many);
 
        if config.plot_walk || config.vid_walk
            clf
            subplot(2,3,1); 
                imshow(wmap); 
                daspect([1 1 1]); 
                hold on;
                plot(pos(:,3),pos(:,2),'k'); 
                axis off
                title('path so far');
                
            subplot(2,3,2); 
                imshow(wmap); 
                daspect([1 1 1]); 
                hold on;
                plot(lpoint0(2),lpoint0(1),'rx');
                plot(lpoint(2),lpoint(1),'ko');                
                plot(ipath(:,2),ipath(:,1),'r');
                axis off
                title('next step');

            subplot(2,3,3); 
                im = imshow(dm); 
                set(im,'alphadata',wmap);
                colormap(gca,jet);                
                daspect([1 1 1]); 
                title('decision');   
                axis off
                
            subplot(2,3,4); 
                im = imagesc(ds); 
                set(im,'alphadata',wmap);                
                daspect([1 1 1]); 
                title('distance');
                axis off
            
            subplot(2,3,5); 
                im = imagesc(as3); 
                set(im,'alphadata',wmap);   
                colormap(gca,jet);
                daspect([1 1 1]); 
                title('angle');
                axis off

            subplot(2,3,6); 
                im = imagesc(occ); 
                set(im,'alphadata',wmap);  
                colormap(gca,jet);                
                daspect([1 1 1]); 
                title('occupancy');
                axis off
                
                drawnow;                
                pause(0.001);
%             keyboard
            if config.vid_walk
                drawnow;
                cframe = getframe(gcf); % take shot
                writeVideo(writerObj,cframe); % add it to video
            end
        end
    end
    if config.vid_walk
        close(writerObj);
    end    

%% Smooth the walk   
    config.posraw = pos;
%     bmap = bwdist(~wmap);
%     posidx = [round(pos(:,2)),round(pos(:,3))];
%     nindx = any(isnan(posidx),2);
%     posidx(nindx,:) = 1;
%     dvals = bmap( sub2ind(size(bmap),posidx(:,1),posidx(:,2)) );
%     dvals(nindx,:) = NaN;
%     jitter_tol = 1;

    pox0 = ( pos(:,3) - size(wmap,2)/2 ).*10./3.5;
    poy0 = ( pos(:,2) - size(wmap,1)/2 ).*10./3.5;

    pox(:) = pox0(:) + (rand(size(pox0(:)))-0.5) .* config.jit;
    poy(:) = poy0(:) + (rand(size(poy0(:)))-0.5) .* config.jit; 
    
    if config.smoo
        d = smoothn({pox poy},config.smoo,'robust');
        pox = d{1};
        poy = d{2};  
    end
    pox = pox(:).*3.5;
    poy = poy(:).*3.5;    
    
    step_time_interval = 60/config.steps_per_min/sample_rate;
    pot = (0:length(pox0)-1) .* step_time_interval;
    pot = pot(:);

    xi = 0 : (1/50) : (config.walk_length.*60);
    pox = interp1(pot,pox,xi,'linear',NaN);
    poy = interp1(pot,poy,xi,'linear',NaN);
    
    pox = pox(:);
    poy = poy(:);
    pot = xi(:);
    
%% Print the walk  
    if config.plot_params
        fig_walk_parameters
    end
    
    


































