%% Main CHAPTER ************ 4 ********************
% Purpose:  This script will call other functions to create the figures
%           generated for Chapter 4 of Tobias Hasse's dissertation.  This
%           script was written to streamline the analysis process but is
%           customized for the simulation data and figure publication sizes
%           for Tobias Hasse's dissertation. When using this on other
%           simulation output data sets, you will need to step through the
%           analysis and make adjustments.
%
% Dependencies: This code requires some output files generated in the
%               process of simulation and analysis for Chapter 2
% Author:   Tobias Hasse tobiack@udel.edu
% Date:     October 2021

%% Variability plots of storage time: Calendar Plot & Violins
% Read in the storage time variables created in extract_triangles from 
% Chapter 2 then create the Calendar plot of storage times and the violin 
% plot of storage times and the autocorrelation of median storage time.
% These plots are slightly different than the dissertation
% plots due to an error in the dissertation plots, explained in the
% Ch4_storage_variability.m code file and Appendix A.1.5

% estimated run time ::::::::::::::::::: 30 seconds

% load in triangle
cd C:\Users\User\Documents\MATLAB\MyLib\Pub\Ch2\other
% fname='TRH 205k upstream storage triangles.mat';
% fname='TRH 205k downstream storage triangles.mat';
fname='TRH 205k up & downstream storage triangles.mat';
load(fname, 'pt_bar_dist', 'vert_dist' )
load('params_storage.mat', 'time_step_years' )

% set the output folder for Chapter 4
cd C:\Users\User\Documents\MATLAB\MyLib\Pub\Ch4\Figures
Ch4_storage_variability( pt_bar_dist, vert_dist, time_step_years )

%% Channel length, deposition and erosion volumes, and cross correlation 
% between channel length and deposition and erosion volumes

% estimated run time ::::::::::::::::::: 2 minutes

clear
cd C:\Users\User\Documents\MATLAB\MyLib\Pub\Ch2\other
load('Hasse_211ka_30yr_A3_Cfo24_2Eo_offset.mat','riv2')
load('Hasse_211ka_30yr_A3_Cfo24_2Eo.mat','riv','B')
% load('code cleanup_30yr_A16_Cfo0.010_2Eo.mat')
% load('code cleanup_30yr_A16_Cfo0.010_2Eo_offset.mat')

load params_storage.mat pix_per_chan x_start x_end time_step_years
load params_reaches.mat

% set the output folder for Chapter 4
cd C:\Users\User\Documents\MATLAB\MyLib\Pub\Ch4\Figures

beep; pause(1); beep

% get the max and min coordinates (in meters) from the floodplain area
limits = domain_limits(riv,riv2,pix_per_chan,B);

% Obtain channel lengths in each reach and create figure Ch4fig25
reach_length = channel_reach_length(riv, B, limits(1), pix_per_chan, ...
    x_start, x_end );

clearvars riv riv2

% Create the cross correlation and deposit volume distribution plots
% Ch4fig26 Ch4fig27, not just for the combined reaches but for each reach
files = {
    'TRH 205k upstream storage triangles.mat'...
    'TRH 205k downstream storage triangles.mat'...
    'TRH 205k up & downstream storage triangles.mat'
    };

for i=1:numel(files)
    
    if isequal(i,3)
        rl = sum(reach_length);
    else
        rl = reach_length(i,:);
    end
    % input folder
    cd C:\Users\User\Documents\MATLAB\MyLib\Pub\Ch2\other
    load(files{i});
    % output folder
    cd C:\Users\User\Documents\MATLAB\MyLib\Pub\Ch4\Figures

length_v_volume( files{i}, rl, time_step_years, ...
    pt_bar_dist, vert_dist, pb_age_dist, vd_age_dist )
disp('One file completed, on to the next')
end 
disp('finished')

%% Create triangular plots of age, storage time, and life expectancy
% load in the data and create the distributions

% estimated run time ::::::::::::::::::: 1 minute

% load in age and storage distributions
cd C:\Users\User\Documents\MATLAB\MyLib\Pub\Ch2\other
fname = 'TRH 205k upstream storage triangles.mat';
load(fname)
load('params_storage.mat', 'time_step_years' )

% set the output folder for Chapter 4
cd C:\Users\User\Documents\MATLAB\MyLib\Pub\Ch4\Figures

% remove terraces ( possibly redundant ) the original storage distribution 
% includes 'terrace' from prior to model start. removed from updated files
pt_bar_dist = fliplr(tril(fliplr(pt_bar_dist),-2)); %-2
pb_age_dist = fliplr(tril(fliplr(pb_age_dist),-2)); %-2 

% convert storage time distributions to CCDF survivor functions
vert_dist   = bsxfun(@rdivide,cumsum(vert_dist,2),sum(vert_dist,2));
pt_bar_dist = bsxfun(@rdivide,cumsum(pt_bar_dist,2),sum(pt_bar_dist,2));
% Forward transit time distributions
ftt=(rot90(spdiags(fliplr(vd_age_dist))));  % rotate data
% make individual ccdfs as well   
vd_ftt = bsxfun(@rdivide,ftt,ftt(:,1));
% for point bar as well
ftt=(rot90(spdiags(fliplr(pb_age_dist))));
pb_ftt = bsxfun(@rdivide,ftt,ftt(:,1));

%% create figures 28-39 of Chapter 4
% this function will save the images in multiple formats and automatically
% close the matlab figures because that shortens run time and avoids memory
% availability problems.  You can edit this behavior within plot_triangles

% estimated run time ::::::::::::::::::: 7 minutes 

plot_triangles(fname, time_step_years, pb_age_dist, vd_age_dist, ...
    pt_bar_dist, vert_dist, pb_ftt, vd_ftt)

% and here my story ends, if you have traveled this far it is certainly
% time to venture in new directions and go beyond the scope of this work!
