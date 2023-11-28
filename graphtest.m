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
% version 1.0.0, Release 12/07/20 Initial release
%
% Author: Roddy Grieves
% UCL, 26 Bedford Way
% eMail: r.grieves@ucl.ac.uk
% Copyright 2019 Roddy Grieves

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ################################################################# %% INPUT ARGUMENTS CHECK
% deal with initial variables
    close all;
    scratch_space = 'C:\Users\F004KS7\OneDrive - Dartmouth College\Projects in prep\2019 Mapping project\associated data\outputs'; % desk computer + bonsai computer
    [~,~,~] = mkdir(scratch_space);
    fig_dir = 'C:\Users\F004KS7\OneDrive - Dartmouth College\Projects in prep\2019 Mapping project\associated media\outputs'; % desk computer + bonsai computer
    fig_dir2 = 'C:\Users\F004KS7\OneDrive - Dartmouth College\Projects in prep\2019 Mapping project\main figures'; % desk computer + bonsai computer
    
    % scratch_space = 'C:\Users\roddy\OneDrive - Dartmouth College\Projects in prep\2019 Mapping project\associated data\outputs'; % desk computer + bonsai computer
    % [~,~,~] = mkdir(scratch_space);
    % fig_dir = 'C:\Users\roddy\OneDrive - Dartmouth College\Projects in prep\2019 Mapping project\associated media\outputs'; % desk computer + bonsai computer
    % fig_dir2 = 'C:\Users\roddy\OneDrive - Dartmouth College\Projects in prep\2019 Mapping project\main figures'; % desk computer + bonsai computer

    try
        cd(scratch_space)
        [~,~,~] = mkdir(fig_dir);        
    catch
        scratch_space = 'C:\Users\admin\Documents\OneDrive - Dartmouth College\Projects in prep\2019 Mapping project\associated data\outputs'; % cheetah computer       
        fig_dir = 'C:\Users\admin\Documents\OneDrive - Dartmouth College\Projects in prep\2019 Mapping project\associated media\outputs'; % cheetah computer
        [~,~,~] = mkdir(fig_dir);        
    end
    
    init = 'C:\Users\F004KS7\OneDrive - Dartmouth College\Projects in prep\2019 Mapping project'; % desktop pc
    % init = 'C:\Users\admin\Documents\OneDrive - Dartmouth College\Projects in prep\2019 Mapping project'; % cheetah pc
    % init = 'C:\Users\F004KS7\OneDrive - Dartmouth College\Projects in prep\2019 Mapping project'; % bonsai PC

    % %% overdispersion settings
    % % close all;
    % scratch_space = [init '\associated data\outputs\overdispersion'];
    % [~,~,~] = mkdir(scratch_space);
    % fig_dir = [init '\associated media\outputs\overdispersion'];
    % fig_dir2 = [init '\main figures\overdispersion'];
    % cd(scratch_space)

    % %% biased settings
    % % close all;
    % scratch_space = [init '\associated data\outputs\biased'];
    % [~,~,~] = mkdir(scratch_space);
    % fig_dir = [init '\associated media\outputs\biased'];
    % fig_dir2 = [init '\main figures\biased'];
    % cd(scratch_space)

    %% error method settings
    % close all;
    scratch_space = [init '\associated data\outputs\errors'];
    [~,~,~] = mkdir(scratch_space);
    fig_dir = [init '\associated media\outputs\errors'];
    fig_dir2 = [init '\main figures\errors'];
    cd(scratch_space)

%% overwrite settings
    overwrite_walks                                     = 0;
    overwrite_fields                                    = 0;    
    overwrite_cells                                     = 0;
    overwrite_maps                                      = 0;    
    append_maps                                         = 0;   
    print_now                                           = 1;
    print_qual                                          = '-r250';
    fig_vis                                             = 'on';
    
%% Functions to run
    funconfig = struct;
    funconfig.MAP_get_random_walks                      = 0; % create random walks (position data)
    funconfig.MAP_get_fields_and_cells                  = 1; % create fields and cells (spike data)
    funconfig.MAP_generate_maps                         = 1; % create firing rate maps and compare to original distribution 
    funconfig.MAP_add_ripley_k                          = 0; % add Ripley's k values (estimate field sizes)
    funconfig.MAP_test_map                              = 0; % small function to run tests on rate_mapper
    funconfig.MAP_fix_datamats                          = 0; % add map parameter matrix to datamat files
    funconfig.MAP_eg_cells                              = 0; % plot example cells from population

%% Figures to make
    figfuns = struct;

    % intro figures
    figfuns.general.MAP_fig_1                           = 0; % Example spikes, trajectories, general method    
    figfuns.general.MAP_fig_ripley_k                    = 0; % Ripley's K method and result 
    
    % summary figures
    figfuns.general.MAP_fig_errors                      = 0; % show MISE, all methods
    figfuns.general.MAP_fig_time_etc                    = 0; %    
    figfuns.general.MAP_fig_summary                     = 0; % compare different methods
    figfuns.general.MAP_fig_maps                        = 0; % show error maps, all methods
    figfuns.general.MAP_fig_v2_pareto                   = 0; % show pareto and regression results, all methods
    figfuns.general.MAP_fig_histogram_lit               = 0; % histogram literature review
    figfuns.general.MAP_fig_overdispersion              = 0; % overdispersion analysis
    figfuns.general.MAP_fig_sampling                    = 0; % biased sampling analysis
    figfuns.general.MAP_fig_error_measures              = 0; % compare different error measures, must be in biased directory
    figfuns.general.MAP_fig_ex_pcells                   = 0; % example real place cells

    % histogram figures
    figfuns.histogram.MAP_fig_1_multi                   = {0,'histogram'}; % example maps, MISE, additional factors, Pareto fronts, regression outcomes
    figfuns.histogram.MAP_fig_ex_multi                  = {0,'histogram'}; % example maps made using pareto-optimal and minimum error settings
    figfuns.histogram.MAP_fig_histogram_meth            = 0; % histogram method for transforming spikes and positions into ratemap
    figfuns.histogram.MAP_fig_histogram_smoo            = 0; % smoothing before vs after division
    
    % ASH figures
    figfuns.ash.MAP_fig_1_multi                         = {0,'ash'};    
    figfuns.ash.MAP_fig_ex_multi                        = {0,'ash'};  
    figfuns.ash.MAP_fig_ash_demo                        = 0;  

    % Skaggs adaptive figures
    figfuns.kadaptive.MAP_fig_1_multi                   = {0,'kadaptive'}; % example maps, MISE, additional factors, Pareto fronts, regression outcomes
    figfuns.kadaptive.MAP_fig_ex_multi                  = {0,'kadaptive'}; % example maps made using pareto-optimal and minimum error settings    
    figfuns.kadaptive.MAP_fig_adaptive_meth             = 0; % adaptive method for transforming spikes and positions into ratemap 
    figfuns.kadaptive.MAP_fig_kadaptive_meth            = 0; % kadaptive method for transforming spikes and positions into ratemap     
    figfuns.kadaptive.MAP_fig_adapt_comparison          = 0; % compare pixelwise adaptive, convolution adaptive methods
    
    % Yartsev adaptive figures
    figfuns.kyadaptive.MAP_fig_1_multi                  = {0,'kyadaptive'}; % example maps, error distribution, additional mapping factors         
    figfuns.kyadaptive.MAP_fig_ex_multi                 = {0,'kyadaptive'}; % example maps made using pareto-optimal and minimum error settings      
        
    % ksde figures
    figfuns.ksde.MAP_fig_1_multi                        = {0,'ksde'}; % example maps, MISE, additional factors, Pareto fronts, regression outcomes
    figfuns.ksde.MAP_fig_ex_multi                       = {0,'ksde'}; % example maps made using pareto-optimal and minimum error settings     
    figfuns.ksde.MAP_fig_ksde_comparison                = 0; % compare leutgeb, pixelwise, matlab KSDE methods
    figfuns.ksde.MAP_fig_ksde_meth                      = 0; % adaptive method for transforming spikes and positions into ratemap 
    
    % fyhn figures
    figfuns.fyhn.MAP_fig_1_multi                        = {0,'fyhn'}; % example maps, MISE, additional factors, Pareto fronts, regression outcomes
    figfuns.fyhn.MAP_fig_ex_multi                       = {0,'fyhn'}; % example maps made using pareto-optimal and minimum error settings      
    figfuns.fyhn.MAP_fig_fyhn_time                      = 0; % smoothing time window, short vs long
        
%% general settings
    config                                              = struct;
    config.environs                                     = {'arena120cm'};
%     config.trial_lengths                                = [4,16,64];      
    config.append                                       = {'errors'};    
    config.nfields                                      = 512; % 512
    config.npcells                                      = 64; % 256
    config.nwalks                                       = 8;
    config.biased_walk                                  = 0; % 0 = off, 1 = biased walk, 2 = thigmotaxis     
    config.plot_fields                                  = 100;
    config.plot_cells                                   = 100;
   
%% cell settings    
    config.frate                                        = [1 1 0.5 10]; % [mean SD min max] normal distribution from which to pull mean firing rate values for cells
    config.field_dist                                   = [5.73,0.26]; % [alpha beta] gamma distribution from which to pull number of field values for cells    
    config.nfactor                                      = [0.3 0.07]; % [mean SD] value added to spike probability, increasing the mean increases the overall probability of spikes, increasing the SD increases the background noise spikes
    config.tfactor                                      = [0 0.02]; % [mean SD] value added to spike times to simulate prospective/retrospective firing
    config.covariance                                   = [4000; 8000; 16000]; % [mean] distribution from which to pull covariance values
    config.overdisperse                                 = [0 1.5 0.8]; % include overdispersion when simulating spikes, time interval of modulation (i.e. time spent in 'on' periods), amplitude of modulation
    
%% generation settings
% trial length x field size parameter pairs
    % config.analysis_pairs = {4000,  [4 16 64];...
    %                          8000,  [4 16 64];...
    %                          16000, [4 16 64]}; % {field size, [durations]}
%     config.analysis_pairs = {4000,  64}; % {field size, [durations]}                         
% % trial length x field size parameter pairs
    % config.analysis_pairs = {16000, 16}; % {field size, [durations]}
    config.analysis_pairs = {16000, 8}; % {field size, [durations]} - overdispersion setup    
%     config.analysis_pairs = {16000,  [4 16 64]}; % {field size, [durations]}                         
        
%% mapping settings
    mapset                                              = struct;
    mapset.padding                                      = 0; % in mm
    mapset.mindwell                                     = 0; % (default 0.1) total number of seconds that rat has to be in a bin for it to be filled, for adaptive method this controls how big the bins will be expanded
    mapset.mindist                                      = 40; % (mm, default 1) used by adaptive and gaussian method, only bins within this distance of some position data will be filled, bigger means more filled (blue) bins
    mapset.maxdist                                      = 640; % (mm, default 50) used by kadaptive, adaptive and gaussian method, only bins within this distance of some position data will be filled, bigger means more filled (blue) bins
    mapset.srate                                        = 50; % (default 50) sampling rate of data in Hz, used to calculate time
    mapset.steps                                        = 32; % the number of convolution size steps to use for kadaptive    
    mapset.kern                                         = 'biweight'; % kernel
    mapset.smethod                                      = 1; % smoothing method, 1 = before division, 2 = after, 3 = no smoothing
    mapset.methods                                      = {'histogram','ash','ksde','kyadaptive','kadaptive','fyhn'}; % histogram ash ksde kyadaptive kadaptive fyhn
    % mapset.methods                                      = {'histogram'};    
    mapset.twindow                                      = 0.25; % time window (s) over which to estimate instantaneous firing for temporal methods
    
    % histogram settings
    mapset.histogram_binsize                            = [1 2.5 5 10:10:100 160 240 320 640]; % binsize in mm
    mapset.histogram_smoothing                          = [0 1 2.5 5 10 15 20 25 30 40 50 60 80 100 160 240 320 640]; % smoothing sigma in mm
    
    % ash settings    
    mapset.ash_binsize                                  = [20 30 40 50 60 70 80 90 100 160 240 320 640]; % binsize in mm
    mapset.ash_smoothing                                = [1 2 3 4 5 6 7 8 9 10 16 24 32 64]; % subdivision count
        
    % kadaptive settings    
    % kadaptive and adaptive methods are interchangeable but kadaptive is
    % the fastest
    mapset.kadaptive_binsize                            = [1 2.5 5 10:10:100 160 240 320 640]; % binsize in mm
    mapset.kadaptive_smoothing                          = [0.1 0.2 0.3 0.4 0.5 0.6 0.8 1 1.5 2 3 4 5 6 8 16 24 32].*1000; % smoothing alpha parameter
    
    mapset.kyadaptive_binsize                           = [1 2.5 5 10:10:100 160 240 320 640]; % binsize in mm
    mapset.kyadaptive_smoothing                         = [0.5:0.5:10]; % smoothing alpha parameter

    % KSDE settings  
    % KSDE, leutgeb method (mvksdensity implementation) and leutgeb method
    % (pixelwise implementation) are interchangeable but KSDE is the
    % fastest
    %[4 and 16 duration sessions]
    mapset.ksde_binsize                                 = [2.5 5 10:10:100 132 160 240 320 640]; % binsize in mm  [2.5 5 10:10:100 132 160 240 320 640] [20:10:100 132 160 240 320 640]
    mapset.ksde_smoothing                               = [10 15 20 25 30 40 50 60 70 80 90 100 120 160 240 320]; % smoothing sigma in mm [10:10:100 120 160 240 320] [10:10:100 120 160 240]
    
    mapset.fyhn_binsize                                 = [5 10:10:100 132 160 240 320 640]; % binsize in mm
    mapset.fyhn_smoothing                               = [10 15 20 25 30 35 40 50 60 70 80 90 100 132 160 240 320 640]; % binsize in mm
        
    
%% ##################################### Heading 2
%% #################### Heading 3
%% Heading 4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ################################################################# %% FUNCTION BODY
%% run the requested analyses
    fnames = fieldnames(funconfig);
    for ff = 1:length(fnames)
        val = funconfig.(fnames{ff});
        if val
            disp(sprintf('\trunning: %s',fnames{ff}))
            eval([fnames{ff} ';']);
        end
    end

%% create the requested figures
    hnames = fieldnames(figfuns);
    for hh = 1:length(hnames)
        %disp(sprintf('\t\trunning %s functions',hnames{hh}))
        
        fnames = fieldnames(figfuns.(hnames{hh}));
        for ff = 1:length(fnames)
            %disp(sprintf('\t\t\t...%s',fnames{ff}))
            
            if iscell( figfuns.(hnames{hh}).(fnames{ff}) )
                val = figfuns.(hnames{hh}).(fnames{ff}){1};
                if val
                    disp(sprintf('\t\t\trunning: %s(%s)',fnames{ff},figfuns.(hnames{hh}).(fnames{ff}){2}))
                    eval(['mname = ''' figfuns.(hnames{hh}).(fnames{ff}){2} ''';']);
                    eval([fnames{ff} ';']);          
                end            
            else
                val = figfuns.(hnames{hh}).(fnames{ff});
                if val
                    disp(sprintf('\trunning: %s',fnames{ff}))
                    eval([fnames{ff} ';']);          
                end            
            end
        end
    end














































