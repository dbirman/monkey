%% Monkey Analysis
%
% The goal here is to set up a set of models that can be fit on arbitrary
% data (from humans or monkeys) for tasks like the DMS/DMC task. We'll fit
% these models to both the human and monkey data to try to understand what
% exactly the humans and monkeys are doing.
%
% DMS Task: Match to Sample
% DMC Task: Match to Category
%
% Models:
%
% 1 Parameter
% Uniform: Choose L/R on a proportion of trials
%       Trivial model
%
% Degrees: "Match" when the second trial is within a certain range
% of the original data, otherwise non-match. Fit with Gaussian (or von
% mises since it's polar angles)
%       Prediction implies that there is sensory noise, predicts that
%       errors will be largest near a particular angle idfference.
%
% Grid Learning: Learn all the pairs of data that are shown.
%       Predcition: this would suggest that any pairs that haven't been
%       seen, or that had errors before, will be less well learned than
%       other pairs. Also it implies that showing a specific pattern of
%       pairs will result in a specific kind of learning.
%
% Model List:
%   Uniform (k changes over time)
%   Degrees (deg changes over time)
%   Grid Learning (g*g choices, change over time, basic RL algorithm)
%   Grid-Angle Learning (g*g choices, n groups, grid learning with the
%   ability to share information within groups defined by angle->angle 
%   ranges)
%
% Fitting
%
%   Each of the models will be fit to the full set of training data in
%   order for each animal individually. Each model will be fit sequentially
%   so Train(1:t-1) Predict(t), the % correct predictions will be used to
%   evaluate model success. Note that these are binary outputs so the
%   models will all represent choice probabilities (m vs. nm) rather than
%   continuous outputs.