function [data_, header]  = monkey_tolong()

header = {'SessNum','Trial','SessTrial','StimLoc','SamDir','TesDir','Lever','Corr','Rewarded','Blink','Valid'};

%% Load Data
load('~/proj/monkey/learning_experiment_behavior_matrix.mat');


data1 = x.beh{1};
data1(:,end+1) = 1;
data2 = x.beh{2};
data2(:,end+1) = 2;

data = [data1;data2];
%% Clean
% Remove invalid trials
data = sel(data,11,1);
data = data(:,[1:10 12]); % remove valid trial col

%% Transform into human space

% Format will be LONGFORM:
% #Hoomans x Trial x SDir x TDir x Resp x Corr x Known x DMS=1

data_ = zeros(size(data,1),8);

for i = 1:size(data,1)
    if data(i,1)<=7
        t = 1;
    else
        t= 2;
    end
    data_(i,:) = [data(i,11) data(i,2) data(i,5) data(i,6) data(i,7) data(i,8) 0 t];
end