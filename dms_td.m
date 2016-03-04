function [ wmat,emat,dif,loss ] = dms_td( data,gamma,figs )
%% TD-Learning Algorithm

% We will initialize a weight matrix 6x6x2 which is the sample dir x test
% direction x choice (press, nopress). On each action we will update the
% choices by using gamma*old + (1-gamma)*correct, based on what the
% data saw and the choice made by the monkey. This weight matrix therevore
% is a 6x6xt matrix of choice probabilities over time, with one parameter
% gamma, which ist he learning rate.

wmat = ones(6,6,size(data,1)+1)/2; % initialize everything to 0.5

for r = 1:size(data,1) % for each row
    dat = data(r,:);
    
    sdir = dat(5);
    tdir = dat(6);
    press = dat(7);
    corr = dat(8);
    
    cfl = [-1 1];
    corr = cfl(corr+1);
    
    for sd = 1:6
        for td = 1:6
            if sdir==-1 || tdir==-1 || sd~=sdir || td~=tdir
                wmat(sd,td,r+1) = wmat(sd,td,r);
            else
                if press
                    wmat(sdir,tdir,r+1) = gamma*wmat(sdir,tdir,r) + (1-gamma)*corr;
                else
                    wmat(sdir,tdir,r+1) = gamma*wmat(sdir,tdir,r) + (1-gamma)*(-1*corr);
                end
            end
        end
    end
end
wmat = wmat(:,:,2:end);

if figs
    %% Figure showing learning rate
    for learned=round(linspace(1,size(data,1),10))
        figure
        imagesc(wmat(:,:,learned))
        colorbar
        caxis([0 1])
    end
end
%% Build empirical TD data

% The idea here is we will build a running average of the probability of
% pressing the bar for every choice combination, we will look in a 360
% trial window and compute at each window the probabilities of all button
% press combinations leading to a press

window = 200;

emat = zeros(6,6,size(data,1));

disppercent(-inf,'Calculating');

evals = (1+window):(size(data,1)-window);
for r = evals
    dat = data(r:r+window,:);
    sdirs = dat(:,5);
    tdirs = dat(:,6);
    press = dat(:,7);
    corr = dat(:,8);
    
    for sd = 1:6
        f1 = sdirs==sd;
        for td = 1:6
            f2 = tdirs==td;
            filter = logical(f1.*f2);
            % filter out only values at sd,rd
            cp = press(filter);
            cc = corr(filter);
            if isempty(cp)
                emat(sd,td,r) = -1;
            else
                emat(sd,td,r) = nanmean(cp); % value = probability of pressing
            end
        end
    end
    disppercent(r/length((1+window):(size(data,1)-window)));
end
disppercent(inf);

%% Figure showing empircal probabilities
if figs
    for learned=round(linspace(1,size(emat,3),10))
        figure
        imagesc(emat(:,:,learned))
        colorbar
        caxis([0 1])
    end
end

%% Diff
dif = wmat(evals)-emat(evals);

%%
% if figs
%     for learned=round(linspace(1,size(dif,3),10))
%         figure
%         imagesc(dif(:,:,learned))
%         colorbar
%         caxis([0 1])
%     end
% end
%%

loss = nansum(dif(:).^2);