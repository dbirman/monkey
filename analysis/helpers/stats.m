function stats( data )
%MONKEY_STATS Compute statistics for the monkey dataset

%%

% Header:

% Format will be LONGFORM:
% #Hoomans x Trial x SDir x TDir x Resp x Corr x Known
tasks = {'DMS','DMC'};
% split by task
for t = 1:2
    disp(sprintf('Data for Task: %s',tasks{t}));
    % split by monkey
    ids = unique(data(:,1));
    for ii = 1:length(ids)
        cid = ids(ii);
        disp(sprintf('Stats for %i',cid));

        cdata = sel(data,1,cid);

        disp(sprintf('Saw %i trials. Valid: %i',max(cdata(:,2)),size(cdata(:,2),1)));

        disp('Performance by third');
        perf = [];
        spc = round(linspace(1,size(cdata,1),4));
        for t = 2:length(spc)
            dat_ = cdata(spc(t-1):spc(t),:);
            perf = [perf mean(dat_(:,6))];
        end
        disp(sprintf('Performance: 1/3 %0.2f 2/3 %0.2f 3/3 %0.2f',perf(1),perf(2),perf(3)));
        
    end
end