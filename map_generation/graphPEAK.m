function [fields,fstats,config] = graphPEAK(ratemap,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ################################################################# %% DESCRIPTION
%graphPEAK Measure properties of image regions.
%   fields = graphPEAK(ratemap) counts the number of fields in the firing rate map
%   ratemap using defaults settings (threshold method at 20% maximum)
%
%   [fields,fstats] = graphPEAK(ratemap) additionally returns an fstats table containing
%   information about each detected place field:
%         'Field'                   field number
%         'Area'                    area (total pixels)
%         'FilledArea'              area (total pixels within perimeter)
%         'Area_cm2'                Area but in cm2
%         'FilledArea_cm2'          FilledArea but in cm2
%         'Centroid'                Center of mass
%         'WeightedCentroid'        Center of mass taking into account firing rate
%         'MaxIntensity'            Maximum rate value in field
%         'MajorAxisLength'         Length of longest axis of fitted ellipsoid
%         'MinorAxisLength'         Length of smallest axis of fitted ellipsoid
%         'Orientation'             Orientation (angle with x axis) of fitted ellipsoid
%         'Eccentricity'            Ratio of the distance between the foci of the ellipse and its major axis length
%         'SNR'                     Signal to noise ratio of field (MaxIntensity / mean map value)
%
%   [fields,fstats,config] = graphPEAK(ratemap) additionally returns a config structure containing
%   settings used to detect place fields, can also be passed to graphPEAK to ensure same settings
%
%   [fields,fstats,config] = graphPEAK(ratemap,____) addtionally use option-value pair arguments:
%          'method'                 [default: 'threshold'] Method to use, choose from 'threshold','denclue','watershed'
%          'lowpass'                [default: 1] Ignore all ratemap bins with a firing rate less than this value
%          'threshold'              [default: 0.2] Firing rate threshold for fields, proportion of peak value in map (i.e. 0.2 = 20% of peak firing)
%          'min_area'               [default: 8] Minimum continguous bins a region must have to be considered a field
%          'max_area'               [default: inf] Maximum continguous bins a region must have to be considered a field
%          'binsize'                [default: 2.5] Binsize in cm used to make the ratemap, used to calculate values in cm and cm2
%
% USAGE:
%   [fields,fstats,config] = graphPEAK(ratemap);
%
% INPUT:
%   ratemap - firing rate map, like the one produced by graphDATA
%   name value pairs - see above
%
% OUTPUT:
%   fields - scalar, number of detected fields
%   fstats - table containing place field info, see above
%   config - structure containing configuration settings
%
% EXAMPLES:
%
%     [fields,fstat,config] = graphPEAK(peaks(256),'method','threshold','lowpass',2.5,'threshold',0.2);
%     figure
%     subplot(1,2,1)
%     imagesc(config.binmap); 
%     daspect([1 1 1]);
%     title('Thresholded map')
%     subplot(1,2,2)
%     imagesc(peaks(256)); 
%     hold on; 
%     daspect([1 1 1]);
%     plot(fstat{:,3}(:,1),fstat{:,3}(:,2),'ko')
%     title('Detected fields')
%
% See also: graphDATA regionprops

% HISTORY
% version 1.0.0, Release 31/10/19 Initial release
%
% Author: Roddy Grieves
% UCL, 26 Bedford Way
% eMail: r.grieves@ucl.ac.uk
% Copyright 2019 Roddy Grieves

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ################################################################# %% INPUT ARGUMENTS CHECK
%% Prepare default settings
    def_config       = struct;
    def_method       = 'threshold';    
    expectedmethods  = {'threshold','denclue','watershed'};
    def_lowpass      = 1;             
    def_threshold    = 0.2;    
    def_min_area     = 36;
    def_max_area     = inf;    
    def_binsize      = 2.5;

%% Parse inputs
    p = inputParser;
    addRequired(p,'ratemap',@(x) ~isempty(x) && isnumeric(x));   
    addOptional(p,'config',def_config,@(x) isstruct(x));  
    addParameter(p,'method',def_method,@(x) any(validatestring(x,expectedmethods)));   
    addParameter(p,'lowpass',def_lowpass,@(x) isnumeric(x) && isscalar(x));
    addParameter(p,'threshold',def_threshold,@(x) isnumeric(x) && isscalar(x));
    addParameter(p,'min_area',def_min_area,@(x) isnumeric(x) && isscalar(x));
    addParameter(p,'max_area',def_max_area,@(x) isnumeric(x) && isscalar(x));    
    addParameter(p,'binsize',def_binsize,@(x) isnumeric(x) && isscalar(x));
    parse(p,ratemap,varargin{:});

%% Retrieve parameters 
    % if a config structure (like the one output by this function or graphDATA) was passed in as the 2nd argument
    % use those settings. Otherwise use the default settings. This is useful if you want to manually
    % specify the settings elsewhere in a structure or if you want to pass a structure to the function 
    % to ensure it uses exactly the same settings
    if ~isempty(fieldnames(p.Results.config))
        config = p.Results.config;
        if ~isfield(config,'ratemap') || isempty(config.ratemap) || all(isnan(config.ratemap(:)))
            config.ratemap = p.Results.ratemap;
        end
    else
        config = p.Results;
    end

%% ##################################### Heading 2
%% #################### Heading 3
%% Heading 4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ################################################################# %% Prepare intial data
%% Prepare spike and position data
    map_mean = mean(ratemap(:),'omitnan');
    fields = 0;
    fstats = table;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ################################################################# %% Detect place fields
    switch config.method
%% ##################################### Threshold method
        case {'threshold'}
            %% REF

            %% DESCRIPTION IN REF      

            %% SIMPLE DESCRIPTION
            
            % some necessary values
            thresh = config.threshold.*max(ratemap(:),[],'omitnan');

            % threshold ratemap to generate a binary map
            binmap = zeros(size(ratemap));
            binmap(ratemap>=thresh) = 1;
            binmap(ratemap<config.lowpass) = 0;
            config.binmap = binmap;
            if ~any(binmap(:))
                return
            end
            binmap = bwlabel(binmap,4); % 4 connectivity means pixels have to share a side, they can't be diagonally connected
            
            % detect regions of interest using regionprops
            fstats = regionprops('table',binmap,ratemap,'Area','FilledArea','Centroid','WeightedCentroid','MaxIntensity','MajorAxisLength','MinorAxisLength','Orientation','Eccentricity');
            
            % delete fields that are too small or too big
            fstats.Area_cm2 = fstats.Area .* (config.binsize.^2);
            fstats(fstats.Area_cm2<=config.min_area,:) = [];
            fstats(fstats.Area_cm2>=config.max_area,:) = [];
            fields = size(fstats,1);
            if nargout==1 % if we only want to count fields, end here
                return
            end
                
            % arrange output and add some additional info
            fstats.Field = reshape(1:fields,[],1);
            fstats = fstats(:,[end 1:end-1]);
            fstats.SNR = fstats.MaxIntensity ./ map_mean;
            fstats.FilledArea_cm2 = fstats.FilledArea .* (config.binsize.^2);

%% ##################################### Watershed method
        case {'watershed'}
            %% REF

            %% DESCRIPTION IN REF      

            %% SIMPLE DESCRIPTION






    end






























































