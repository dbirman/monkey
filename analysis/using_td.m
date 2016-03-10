%% Temporal Difference Learning Model
% This script will load the human and monkey data that it has available and
% attempt to fit a TD learning model (with some params). It will output a
% few specific figures:
% The max likelihood fit based on the best params??
% The learning matrix over time for the model fit vs. actual data

%% List of TD Learning models
% Model 1: Pure Rescorla-Wagner
% The simplest possible model is that each pair of stimuli gives you a unit
% of information, and you use that to update the transition probabilities
% in the matrix
% Model 2: Bleed Learning
% When you learn a given pair that information 'bleeds' to other pairs in
% the matrix
% Model 3: 

%% Convert all data to long-form

% Header:

% Format will be LONGFORM:
% #Hoomans x Trial x SDir x TDir x Resp x Corr x Known x DMS=1
[mdata,mheader] = monkey_tolong();
[hdata,header] = human_tolong();

%% Compute Statistics
stats(mdata);
stats(hdata);

%% DMS or DMC
data = fil(data,1,'>',7); % for DMC
%data = fil(data,1,'<=',7); % for DMS
