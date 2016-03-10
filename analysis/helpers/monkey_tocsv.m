function monkey_tocsv()
load('~/proj/monkey/learning_experiment_behavior_matrix.mat');
header = {'SessNum','Trial','SessTrial','StimLoc','SamDir','TesDir','Lever','Corr','Rewarded','Blink','Valid'};

%% Re-write as csv
fname1 = '~/proj/monkey/wahwah.csv';
fname2 = '~/proj/monkey/quincy.csv';

data1 = x.beh{1};
data2 = x.beh{2};

csvwriteh(fname1,data1,header);
csvwriteh(fname2,data2,header);