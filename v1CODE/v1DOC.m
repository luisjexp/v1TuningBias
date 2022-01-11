%% INITIALIZE PATH AND CLASS
clear;
clc;
close all;
addpath /Users/luis/Box/prjV1TB/v1DATA/   
addpath /Users/luis/Box/prjV1TB/v1CODE/
cd /Users/luis/Box/prjV1TB/v1CODE

% Instantiate the class 
v1  = v1anz();

%% - - - -  - - - -  - - - -  - - - -  - - - -  - - - -  - - - -  - - - - 
%% - - - -  TABLE OF V1 UNIT PROPERTIES- - - - 
%% LOAD A TABLE OF V1 UNITS
clf
v1.get_roi_table

%% TAKE A LOOK AT THE TABLE
clc
fprintf('\nData of the V1 imaging fields that was collected\n')
disp(v1.roi_stack_uif_list()')


fprintf('\n\nThe table of V1 units and their properties for each imaging field\n\n')
head(v1.roi_stack)
fprintf('\n\n\t...\n\n')
tail(v1.roi_stack)



clc
close all
figure('WindowStyle', 'docked')
figure(gcf)
v1.PLT_single_units(v1.roi_stack)

%% - - - -  - - - -  - - - -  - - - -  - - - -  - - - -  - - - -  - - - - 
%% - - - - COMPUTE SIMILARITY AND GENERATE TABLE
clc
v1.get_distance_table();
%% TAKE A LOOK
clc
fprintf('\n\nThe Similarity table, containing properties of *[%d]* pairs V1 units\n\n',...
    height(v1.dist_stack))

head(v1.dist_stack)
tail(v1.dist_stack)

%% - - - -  - - - -  - - - -  - - - -  - - - -  - - - -  - - - -  - - - - 
%% - - -  ASSESS RELATIONSHIP BETWEEN TUNING SIMILARITY AND RF OVERLAP 

% TUNING SIMILARITY ~ SUBREGION OVERLAP
clc
close all
figure('WindowStyle', 'docked')

% We want to show the relationship between tuning similarity and
% ON-subregion overlap and tuning similarity and OFF-subregion overlap
all_uifs    = v1.roi_stack_uif_list;
dY_name     = 'dTunKern_corr';

% Plot the relationship for Tuning Similiraty ~ ON subregion overlap
tbl_on_only = v1.get_distance_table(v1.roi_stack_uif_list,'ON_ONLY' );
v1.fit_plot_lm(tbl_on_only, dY_name, 'dON_corr',...
    Marker = '.', MarkerSize = 5, MarkerEdgeColor = 'r');
ax = gca;

%
% Plot the relationship for Tuning Similiraty ~ OFF subregion overlap
tbl_off_only = v1.get_distance_table(v1.roi_stack_uif_list,'OFF_ONLY' );
v1.fit_plot_lm(tbl_off_only, dY_name, 'dOFF_corr',...
     Marker = '.', MarkerSize = 5, MarkerEdgeColor = 'b');

% - - - -  - - - -  - - - -  - - - -  - - - -  - - - -  - - - -  - - - - 
% TUNING SIMILARITY ~ SUBREGION OVERLAP, ONLY CELLS W BOTH SUBREGIONS
figure(gcf);

tbl_off_and_on = v1.get_distance_table(v1.roi_stack_uif_list,'OFF_AND_ON' );
v1.fit_plot_lm(tbl_off_and_on, dY_name, 'dON_corr',...
     Marker = 'o', MarkerSize = 12, MarkerEdgeColor = 'r');    
ax = gca;

v1.fit_plot_lm(tbl_off_and_on, dY_name, 'dOFF_corr',...
     Marker = 'o', MarkerSize = 12, MarkerEdgeColor = 'b');    
ax = gca;


% Label the figures
set(ax.XLabel, 'String', 'Overlap Between Subregions')
set(ax.YLabel, 'String', 'Tuning Similarity')
legend(ax, {'ON-Overlap', 'OFF-Overlap',...
    'Units w/ both ON+OFF Subregions)','Units w/ both ON+OFF Subregions)'},...
    'Location', 'best')
title('Units with overlapping subregions show similar tuning profiles')





