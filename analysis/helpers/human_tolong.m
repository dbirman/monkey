function [ data, header ] = human_tolong()
%% Load hooman DMC data
files = dir('~/proj/monkey/human/dmc/*.mat');

% Format will be LONGFORM:
% #Hoomans x Trial x SDir x TDir x Resp x Corr x Known
data = zeros(length(files)*360,7);

pos = 1;
for fi = 1:length(files)
    fname = fullfile('~/proj/monkey/human/dmc',files(fi).name);
    load(fname);
    
    if length(jglData.responses)<360
        jglData.responses(jglData.responses==-1) = 0;
        jglData.correct(jglData.correct==-1) = 0;

        for t = 1:length(jglData.rot1)
            dat = [fi, t, jglData.rot1(t)/pi*4+1, jglData.rot2(t)/pi*4+1, jglData.responses(t), jglData.correct(t), jglData.known(t)];
            data(pos,:) = dat;
            pos = pos+1;
        end
    end
end
data = data(1:pos-1,:);

data(:,end+1) = 2;

header = {'Human','Trial','Rot1','Rot2','Response','Correct','Known','Task'};
