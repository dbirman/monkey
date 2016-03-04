%% Monkey Analysis Script

load('~/proj/monkey/learning_experiment_behavior_matrix.mat');

% Header:
%    1        2         3          4           5         6        7
% Session - Trial - SessTrial - StimLoc - SampleDir - TestDir - Lever - 
%    8        9         10         11
% Correct - Rewarded - Blink - Valid Trial

header = {'SessNum','Trial','SessTrial','StimLoc','SamDir','TesDir','Lever','Corr','Rewarded','Blink','Valid'};
%% Re-write as csv
fname1 = '~/proj/monkey/wahwah.csv';
fname2 = '~/proj/monkey/quincy.csv';

data1 = x.beh{1};
data2 = x.beh{2};

csvwriteh(fname1,data1,header);
csvwriteh(fname2,data2,header);

%% Load Data

data1 = x.beh{1};
data1(:,end+1) = 1;
data1 = fil(data1,1,'>',7); % DMC
data1(:,2) = data1(:,2)-min(data1(:,2))+1; % revert trial to 1
data2 = x.beh{2};
data2(:,end+1) = 2;
data2 = fil(data2,1,'>',7); % DMC
data2(:,2) = data2(:,2)-min(data2(:,2))+1;

data = [data1;data2];
%% Clean
% Remove invalid trials
data = sel(data,11,1);
data = data(:,1:10); % remove valid trial col
% x = fil(data,1,'<=',7); % for DMS
data = fil(data,1,'>',7); % for DMC

%% Compute, for trial window, average across monkeys, DMC performance.
%%
max_t = max(data(:,2));

% Now the idea is to build a big matrix that is trial x rot1 x rot2 where
% each value is the across-subject proportion of pressing for that
% rotation. We don't have enough data on each trial, so we bin across 9
% trials
win = 72;

out = zeros(max_t,6,6);
s = win+1;
e = max_t-win;

for t = s:e
    t_ = t-win;
    t__ = t+win;
    
    % select out relevant trials
    dat = fil(data,2,'>=',t_);
    dat = fil(dat,2,'<=',t__);
    % select out relevant directions
    for r1 = 1:6
        for r2 = 1:6
            dat_ = sel(dat,5,r1);
            dat_ = sel(dat_,6,r2);
            press_prop = nanmean(dat_(:,7));
            out(t,r1,r2) = press_prop;
        end
    end
end
% % %% Smooth with 1-d gaussian 
% % gp = [1,0,1,0];
% % pos = -3:3;
% % g = [];
% % for pi = 1:length(pos)
% %     g(pi) = gauss(gp,pos(pi));
% % end
% % g = g/sum(g);
% % % requires 3 pad
% % for i = 1:4
% %     cout = conv(hmat(:,i),g);
% %     hmat(:,i) = cout(4:end-3);
% % end

%% Plot moonkey DMC data!
figure
colormap winter
for i = s:50:e
    imagesc(squeeze(nanmean(out(i:i+50,:,:))));
    colorbar
    caxis([0 1])
    pause(.0000001);
end

%% Write DMC gif
e = 12500;
out_ = round(out*100);
c = 0;
for i = s:50:e-50
    fname = sprintf('~/proj/monkey/moonkey/img%03d.png',c);
    c = c + 1;
    imwrite(imresize(squeeze(nanmean(out_(i:i+50,:,:))),50,'nearest'),winter(100),fname,'png');
    disp(c);
end
lf = imresize(squeeze(nanmean(out_(i:i+50,:,:))),50,'nearest');
for j = 0:100
    fname = sprintf('~/proj/monkey/moonkey/img%03d.png',c+j);
    imwrite(lf,winter(100),fname,'png');
        disp(c+j);

end

%% Plot
plot(data(:,1),data(:,8),'*');

%%
gs = linspace(.9,.999999,10);
loss = zeros(1,length(gs));
parfor gi = 1:length(gs)
    gamma = gs(gi);
    [~,~,~,loss(gi)] = dms_td(x,gamma,0);
end
[m,i] = min(loss);
disp(sprintf('Best Loss %4.2f gamma = %0.4f',m,gs(i)));

%%
best = gs(i);

%% Compute with best

[wmat, emat, dif, loss] = dms_td(x,.5,0);

%% Plot EMAT?
    figure
i = 1
for learned=round(linspace(1,size(emat,3),10))
    subplot(11,1,i);
    i = i+1;
    imagesc(emat(:,:,learned))
    colorbar
    caxis([0 1])
end

%% plot specifc emat value
figure
imagesc(emat(:,:,end-1000));
colorbar
caxis([0 1]);

%% Monkey Video
figure
colormap winter
for i=1:size(emat,3)
    imagesc(emat(:,:,i));
    caxis([0 1])
    pause(.000000001);
end

%% Print monkey video
out_ = round(emat*100);
s = 201; e = 600;
for i = s:e
    fname = sprintf('~/proj/monkey/moonkey/img%03d.png',i-s);
    imwrite(imresize(squeeze(out_(:,:,i)),50,'nearest'),winter(100),fname,'png');
    disp(i);
end
for j = 1:50
    fname = sprintf('~/proj/monkey/moonkey/img%03d.png',i-s+j);
    imwrite(imresize(squeeze(out_(:,:,i)),50,'nearest'),winter(100),fname,'png');
end

%% monkey mask

mask = [1 0 0 1 0 0
    0 1 0 0 1 0
    0 0 1 0 0 1
    1 0 0 1 0 0
    0 1 0 0 1 0
    0 0 1 0 0 1];

mask = 1-mask;
imagesc(mask)

colormap gray

%% Load hooman DMC data
files = dir('~/proj/monkey/human/dmc/*.mat');

% Format will be LONGFORM:
% #Hoomans x Trial x SDir x TDir x Resp x Corr x Known
hooman_data = zeros(length(files)*360,7);

pos = 1;
for fi = 1:length(files)
    fname = fullfile('~/proj/monkey/human/dmc',files(fi).name);
    load(fname);
    
    if length(jglData.responses)<360
        jglData.responses(jglData.responses==-1) = 0;
        jglData.correct(jglData.correct==-1) = 0;

        for t = 1:length(jglData.rot1)
            dat = [fi, t, jglData.rot1(t)/pi*4+1, jglData.rot2(t)/pi*4+1, jglData.responses(t), jglData.correct(t), jglData.known(t)];
            hooman_data(pos,:) = dat;
            pos = pos+1;
        end
    end
end
hooman_data = hooman_data(1:pos-1,:);

%%
max_t = max(hooman_data(:,2));

% Now the idea is to build a big matrix that is trial x rot1 x rot2 where
% each value is the across-subject proportion of pressing for that
% rotation. We don't have enough data on each trial, so we bin across 9
% trials
n = length(hooman_data);
win = 14;
% #Hoomans x Trial x SDir x TDir x Resp x Corr x Known

% the last two dimensions are for data, # of subjects, and # of examples
out = zeros(max_t,8,8);
s = win+1;
e = max_t-win;

use_data = fil(hooman_data,1,65);

for t = s:e
    t_ = t-win;
    t__ = t+win;
    
    % select out relevant trials
    dat = fil(hooman_data,2,'>=',t_);
    dat = fil(dat,2,'<=',t__);
    % select out relevant directions
    for r1 = 1:8
        for r2 = 1:8
            dat_ = sel(dat,3,r1);
            dat_ = sel(dat_,4,r2);
            press_prop = nanmean(dat_(:,5));
            out(t,r1,r2) = press_prop;
        end
    end
end

%% Plot hooman DMC data!
figure
colormap winter
for i = s:e
    imagesc(squeeze(out(i,:,:)));
    colorbar
    caxis([0 1])
    pause(.001);
end

%% Write DMC gif
out_ = round(out*100);
for i = s:e
    fname = sprintf('~/proj/monkey/hooman/img%03d.png',i-s);
    imwrite(imresize(squeeze(out_(i,:,:)),50,'nearest'),winter(100),fname,'png');
    disp(i);
end
for j = 1:100
    fname = sprintf('~/proj/monkey/hooman/img%03d.png',i-s+j);
    imwrite(imresize(squeeze(out_(i,:,:)),50,'nearest'),winter(100),fname,'png');
end

%% Plot ideal

ideal = [1,1,1,1,0,0,0,0
    1,1,1,1,0,0,0,0
    1,1,1,1,0,0,0,0
    1,1,1,1,0,0,0,0
    0,0,0,0,1,1,1,1
    0,0,0,0,1,1,1,1
    0,0,0,0,1,1,1,1
    0,0,0,0,1,1,1,1];

imwrite(imresize(ideal,50,'nearest')*100,winter(100),'~/proj/repo/Grad School/ideal_hooman.png','png');

imagesc(ideal)
colormap(winter)
colorbar

%% Load Hooman data
files = dir('~/proj/monkey/human/dms/*.mat');

% Format will be:
% #Hoomans x Trial x SDir x TDir x Dif x Resp x Corr
hooman_data = zeros(length(files),100,4);

for fi = 1:length(files)
    fname = fullfile('~/proj/monkey/human/dms',files(fi).name);
    load(fname);
    
    hooman_data(fi,:,1) = jglData.rot1/pi*4+1;
    hooman_data(fi,:,2) = jglData.rot2/pi*4+1;
    hooman_data(fi,:,3) = round(abs(jglData.rot1-jglData.rot2)/pi*8+1);
    hooman_data(fi,:,4) = jglData.responses;
    hooman_data(fi,:,5) = jglData.correct;
end

% split by training vs. known
hd_train = hooman_data(:,1:75,:);
hd_known = hooman_data(:,76:end,:);

%% Build Hooman Video

hmat = zeros(75,4);
for t = 1:75
    dat = squeeze(hd_train(:,t,:));
    
    for dif = 1:4
        cdat = sel(dat,3,dif);
        
        if ~isempty(cdat)
            press_prop = sum(cdat(:,4)==1)/size(cdat,1);
            hmat(t,dif) = press_prop;
        else
            if t > 1
                hmat(t,dif) = hmat(t-1,dif);
            else
                hmat(t,dif) = 0.5;
            end
        end
    end
end

%% Smooth with 1-d gaussian 
gp = [1,0,1,0];
pos = -3:3;
g = [];
for pi = 1:length(pos)
    g(pi) = gauss(gp,pos(pi));
end
g = g/sum(g);
% requires 3 pad
for i = 1:4
    cout = conv(hmat(:,i),g);
    hmat(:,i) = cout(4:end-3);
end
%%
figure

%%
for i = 1:125
imagesc(hmat(i,:));
colorbar
caxis([0 1])
pause(.01);
end

%% Simulate human data using RL

% that is, given the sequence of events each human saw, try to reconstruct
% (based on a btwn-subj learning rate) the overall pattern of results.
% Compare this to a 'jump' in learning.