function [r,d1,d2] = ripleyk(pos,spk,binsize,padding)
%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> DESCRIPTION
% FUNCTION  short desc.
% long description
%
% USAGE:
%       [out] = template(in) process with default settings
% 
%       [out] = template(in,optional1) process using optional argument 1
% 
%       [out] = template(___,Name,Value,...) process with Name-Value pairs used to control aspects 
%       of the process
% 
%       Parameters include:
% 
%       'param1'          -   (default = X) Scalar value, parameter to do something
% 
%       'param2'          -   (default = X) Scalar value, parameter to do something
% 
% INPUT:
%       in    - input as a vector
%
% OUTPUT:
%       out   - output as a vector
%
% EXAMPLES:
%       % run function using default values
%       out = template(in,varargin)
%
% See also: FUNCTION2 FUNCTION3

% HISTORY:
% version 1.0.0, Release 00/00/00 Initial release
%
% Author: Roddy Grieves
% Dartmouth College, Moore Hall
% eMail: roddy.m.grieves@dartmouth.edu
% Copyright 2021 Roddy Grieves

%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Heading 3
%% >>>>>>>>>>>>>>>>>>>> Heading 2
%% >>>>>>>>>> Heading 1
%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> INPUT ARGUMENTS CHECK
    if ~exist('binsize','var') || isempty(binsize) || isnan(binsize) || ~isnumeric(binsize) 
        binsize = 10; % bin size in mm
    end
    if ~exist('padding','var') || isempty(padding) || isnan(padding) || ~isnumeric(padding) 
        padding = 5; % padding in bins
    end

%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> FUNCTION BODY
    % prepare point data
    pox = double( pos(:,1) );
    poy = double( pos(:,2) );    
    spx = double( spk(:,1) );
    spy = double( spk(:,2) );

    %% create rough bivariate histograms of point data
    xedg = min(pox) : binsize : max(pox);
    yedg = min(poy) : binsize : max(poy);
    pos_binned = histcounts2(poy,pox,yedg,xedg);
    spk_binned = histcounts2(spy,spx,yedg,xedg);

    % pad the maps
    pos_binned = padarray(pos_binned,[padding padding],0,'both');
    spk_binned = padarray(spk_binned,[padding padding],0,'both');

    % find valid positions 
    % valid positions are regions where it is possible for spikes to have been
    % emitted (i.e. bins containing position data)
    bounds = double( pos_binned>0 );

    % prepare important values
    % these will be used in the radius loop and it saves some computation to
    % prepare them here instead
    n = numel(spx); % total (n)umber of spikes
    A = sum(bounds(:)) .* (binsize.^2); % total (A)rea of field
    E = A ./ (n.*(n-1)); % (E)xpected value/density of points
    radius_bins = 3 : 2 : max([range(pox) range(poy)]./binsize); % radii we want to test in bins
    radius_bins = radius_bins(:);

%% >>>>>>>>>> run through different radii and calculate Ripley's K   
    % for each radius we will use a convolution window with that size, by doing
    % the computation through convolution instead of geometrically we can make
    % a significant speed increase
    radius_mm = radius_bins(:) .* binsize; % each radius in mm
    area_mm2 = pi.*(radius_mm(:).^2); % each radius area in mm2  
    
    ripley_k = NaN(length(radius_bins),1);
    for rr = 1:length(radius_bins) % for every test radius
        % create a circular convolution kernel
        kern_now = strel('disk',radius_bins(rr));
        kern_now = double(kern_now.Neighborhood);    

        % convolve the valid position map with this kernel
        % for every pixel this will count the number of pixels within
        % the test radius that are valid positions, then we convert it to
        % area in mm2
        area_map = imfilter(bounds,kern_now,0,'same','conv') .* (binsize.^2);

        % convolve the spike map with this kernel (ignoring the centre)   
        % for every pixel (p) this will count the number of spikes within the
        % test radius (excluding spikes in p)
        kern_now( ceil(size(kern_now,1)/2) , ceil(size(kern_now,1)/2) ) = 0;
        spk_in_r = imfilter(spk_binned,kern_now,0,'same','conv');

        % calculate Ripley's K for this radius 
        ripley_k(rr,1) = E .* sum(sum( area_mm2(rr,1) ./ area_map .* spk_binned .* (spk_binned - 1 + spk_in_r),'omitnan'),'omitnan'); % Ripley's k  
    end

    % this is the extended method reported by Amgad et al. (2015) that
    % is adapted to work on grayscale images and includes Besag's edge
    % correction. We also normalise the K value to get the H-function      
    besag_norm = sqrt( ripley_k/pi ); % Besag's normalization
    h_function = besag_norm - radius_mm; % H-function     
    
%% >>>>>>>>>> calculate the derivative of H
    % Kiskowski et al. (2009) report using this as suggested by Weigand and Moloney (2004): 
    % "although the K-function may be used to detect the spatial range of repulsive 
    % and attractive effects, functions based on the derivative of the K-function 
    % should be used to describe the extent of aggregation at a particular distance.
    % This is the case because the rate of change of a function does not depend on the
    % function's magnitude, and thus is not affected by accumulative effects (21). 
    % We therefore tested the use of the derivative of H(r) to subtract
    % accumulative effects."
    h_derivative = smooth( gradient(h_function) ./ gradient(radius_mm(:)), 3, 'moving');

    % combine everything into one output table
    r = table(radius_bins,radius_mm,area_mm2,ripley_k,besag_norm,h_function,h_derivative);

%% >>>>>>>>>> estimate domain size
    % Kiskowski et al. (2009):
    % "In previous studies, the value yielding the maximum of H(r) was used to provide a measure 
    % of domain radius"
    % but they also note that:
    % "the domain radius R∗ predicted by [HMAX] monotonically increases from the domain radius 
    % R to the domain diameter 2R as the separation is increased from 4R to arbitrarily large 
    % values. We evaluated the radius of maximal aggregation predicted as S → ∞ by applying 
    % Ripley's K to a single domain without periodic boundary conditions, and found that 
    % [HMAX] approaches 2R."
    [~,midx] = max(h_function,[],'omitnan'); % HMAX
    d1 = radius_mm(midx) / 2; % divide by 2 because place fields are likely seperated a lot  

    % Kiskowski et al. (2009):    
    % "Thus, the minimum value of r that yields H′(r) = −1, divided by 2, accurately yields 
    % the domain radius independent of the domain separation for the case of idealized domains."
    % However, I have found that the derivative of H does not always fall below -1 so this method
    % does not seem widely applicable. Instead I am just using the location of the minimum.
    idx1 = find(h_derivative<0,1,'first');
    [~,idx2] = min(h_derivative,[],'omitnan'); % 
    %d2 = radius_mm(idx2) / 2; % divide by 2 because place fields are likely seperated a lot  
    idx3 = mean([idx1 idx2]);
    d2 = interp1(1:length(radius_mm),radius_mm,idx3,'linear') / 2;

































