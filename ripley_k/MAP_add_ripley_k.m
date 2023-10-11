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
    wlength = 64;
    
%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> FUNCTION BODY
    for ee = 1:length(config.environs) % for every environment        
        field_sizes_now = [config.analysis_pairs{:,1}];
        disp(sprintf('\tEnvironment: %s...',config.environs{ee}));      

%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> for every field size requested           
        for ff = 1:length(field_sizes_now) 
            mean_ssnow = field_sizes_now(ff);
            disp(sprintf('\t\tField size: %d...',mean_ssnow));      

            sdir = [scratch_space '\' config.environs{ee} '_sigma' num2str(mean_ssnow) '_duration' num2str(wlength) '_sdata.mat'];            
            disp(sprintf('\t\t\t...loading cells: %s',sdir));            
            load(sdir,'sdata');  
            sdata.estimated_field_size = NaN(size(sdata,1),2);

%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> for every cell     
            loopout = looper(size(sdata,1));           
            for uu = 1:size(sdata,1)
                % get the position and spike data
                pindx = sdata.pos_index(uu);
                pos_now = all_pos{pindx,1};
                pos = pos_now(:,1:2);
                spk = sdata.spk{uu};

                [r,d1,d2] = ripleyk(pos,spk,10);
                sdata.estimated_field_size(uu,:) = [d1 d2];
                loopout = looper(uu);
            end
            save(sdir,'sdata','-v7.3');              
        end
    end







