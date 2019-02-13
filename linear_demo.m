% DESCRIPTION: Demo code for running seqNMF on simulated and real data,
% including how to test significance of each factor on held-out data, and
% how to select lambda
% 
% ------------------------------------------------------------------------
% Andrew Bahle and Emily Mackevicius 1.26.2018
%
% See paper: 
% https://www.biorxiv.org/content/early/2018/03/02/273128
%% Vars
COOKED_DATA_PATH = 'Z:\experiments\projects\bambi\linear_track_2\analysis\1_linear_track\c51m4\day_1\cooked';
NUMBER_OF_FRAMES = 2000; % Choose a small number of frames for demo to work fast
VIDEOfs = 10;

%% load data from linear track
load([COOKED_DATA_PATH, '\finalTracesMat']);
neural = abs(allTracesMat(1:NUMBER_OF_FRAMES, :))';
behavior = csvread([COOKED_DATA_PATH, '\behavior.csv'], 1, 1);
bins = behavior(:, 3);
bins = bins(1:NUMBER_OF_FRAMES)';
display('loaded data')
%% break data into training set and test set
splitN = floor(size(neural,2)*.75); 
train_neural = neural(:,1:splitN); 
train_bins = bins(:,1:splitN); 
test_neural = neural(:,(splitN+1):end); 
test_bins = bins(:,(splitN+1):end); 
%% plot one example factorization
rng(235); % fixed rng seed for reproduceability
X = train_neural;
K = 16;
L = 1/10; % units of seconds
Lneural = ceil(L*VIDEOfs);
shg
display('Running seqNMF on real neural data')
[W, H, ~,loadings,power]= seqNMF(X,'K',K,'L',Lneural,...
            'lambdaL1W', .1, 'lambda', .015, 'maxiter', 100, 'showPlot', 0,...
            'lambdaOrthoW', 0); 
p = .05; % desired p value for factors

display('Testing significance of factors on held-out data')
[pvals,is_significant] = test_significance(test_neural,W,p);

W = W(:,is_significant,:); 
H = H(is_significant,:); 

%% plot, sorting neurons by latency within each factor
[max_factor, L_sort, max_sort, hybrid] = helper.ClusterByFactor(W(:,:,:),1);
indSort = hybrid(:,3);
tstart = 1; % plot data starting at this timebin
figure; WHPlot(W(indSort,:,:),H(:,tstart:end), X(indSort,tstart:end), ...
    1,train_bins)
title('Significant seqNMF factors, with raw data')
figure; WHPlot(W(indSort,:,:),H(:,tstart:end), ...
    helper.reconstruct(W(indSort,:,:),H(:,tstart:end)),...
    1,train_bins)
title('SeqNMF reconstruction')

%% Procedure for choosing lambda
nLambdas = 20; % increase if you're patient
X = train_neural;
lambdas = sort([logspace(-1,-5,nLambdas)], 'ascend'); 
loadings = [];
regularization = []; 
cost = []; 
for li = 1:length(lambdas)
    [N,T] = size(X);
    [W, H, ~,loadings(li,:),power]= seqNMF(X,'K',K,'L',Lneural,...
        'lambdaL1W', .1, 'lambda', lambdas(li), 'maxiter', 100, 'showPlot', 0); 
    [cost(li),regularization(li),~] = helper.get_seqNMF_cost(X,W,H);
    display(['Testing lambda ' num2str(li) '/' num2str(length(lambdas))])
end
%% plot costs as a function of lambda
windowSize = 3; 
b = (1/windowSize)*ones(1,windowSize);
a = 1;
Rs = filtfilt(b,a,regularization); 
minRs = prctile(regularization,10); maxRs= prctile(regularization,90);
Rs = (Rs-minRs)/(maxRs-minRs); 
R = (regularization-minRs)/(maxRs-minRs); 
Cs = filtfilt(b,a,cost); 
minCs =  prctile(cost,10); maxCs =  prctile(cost,90); 
Cs = (Cs -minCs)/(maxCs-minCs); 
C = (cost -minCs)/(maxCs-minCs); 

clf; hold on
plot(lambdas,Rs, 'b')
plot(lambdas,Cs,'r')
scatter(lambdas, R, 'b', 'markerfacecolor', 'flat');
scatter(lambdas, C, 'r', 'markerfacecolor', 'flat');
xlabel('Lambda'); ylabel('Cost (au)')
set(legend('Correlation cost', 'Reconstruction cost'), 'Box', 'on')
set(gca, 'xscale', 'log', 'ytick', [], 'color', 'none')
set(gca,'color','none','tickdir','out','ticklength', [0.025, 0.025])

%% chose lambda; run multiple times, see number of sig factors
loadings = [];
pvals = []; 
is_significant = []; 
X = train_neural;
nIter = 20; % increase if patient
lambda = 0.015;
display(['Running seqNMF multiple times for lambda=', num2str(lambda)])

for iteri = 1:nIter
    [W, H, ~,loadings(iteri,:),power]= seqNMF(X,'K',K,'L',Lneural,...
            'lambdaL1W', .1, 'lambda', lambda, 'maxiter', 100, 'showPlot', 0); 
    p = .05;
    [pvals(iteri,:),is_significant(iteri,:)] = test_significance(test_neural,W,p);
    W = W(:,is_significant(iteri,:)==1,:); 
    H = H(is_significant(iteri,:)==1,:); 
%     [max_factor, L_sort, max_sort, hybrid] = helper.ClusterByFactor(W(:,:,:),1);
%     indSort = hybrid(:,3);
%     tstart = 300; 
%     figure; WHPlot(W(indSort,:,:),H(:,tstart:end), X(indSort,tstart:end), 1, train_bins)
%     display(['seqNMF run ' num2str(iteri) '/' num2str(nIter)])
end
%%
figure; hold on
h = hist(sum(is_significant,2), 'edgecolor', 'w', 'facecolor', .7*[1 1 1]); 
h.BinCounts = h.BinCounts/sum(h.BinCounts)*100; 
xlim([0 10]); 
xlabel('# significant factors')
ylabel('% seqNMF runs')

%% Plot factor-triggered song examples and rastors
addpath(genpath('misc_elm')); 
figure; HTriggeredSpec(H,train_bins,VIDEOfs,SONGfs,Lsong); 

figure; HTriggeredRaster(H,train_neural(indSort,:),Lneural);

