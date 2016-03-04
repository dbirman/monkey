function [ wmat,emat,dif,loss ] = dms_ju( data,time,alpha,figs )
%% Jump algorithm

% the assumption behind this is that at some point the monkey's behavior
% "jumps", i.e. it goes from not understanding the task to fully
% understanding the task. This is what it feels like when humans do the
% task. This has two parameters, an error rate (alpha) and a jump time when
% the monkey "learns" the task.

% Reps is an internal parameter that determines how many simulations of the
% alpha rate / time to run

[T,~] = size(data);

reps = 100;

% Let's build a random dataset of choice probabilities based on our alpha
% and time

% Stay / switch probabilities are just .5 at the beginning, and then 1 once
% the algorithm learns"

for r = 1:reps
    % build probability matrix
    probs = zeros(2,size(data,1));
    
    probs(:,1:time) = .5 + alpha * 2 * rand(2,time) - alpha / 2;
    probs(1,time+1:end) = 1 - alpha*rand(1,T-time);
    probs(2,time+1:end) = 0 + alpha*rand(1,T-time);
    
    % simulate this probability matrix
end