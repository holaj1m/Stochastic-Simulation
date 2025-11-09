%% weighted histogram analysis method (WHAM)
% The idea of this code is to simulate data for a particle
% under the anharmonic potential \phi(x) = -0.5*(x-1)^3 at
% temperature T. 
% We use the canonical model and Métropolis method to sampling
% the position of particle. However the nature of potential
% energy makes very difficult to sample the space properly. 
% So to fix this issue we divided the domain in different 
% intervals. Given that each interval has its own distribution
% we have to plug the intervals together using different factor.
% The method described above corresponds to WHAM.

% Create anonymous function with potential energy
potentialEnergy = @(x) 0.5*(x-1).^3;

% Prob. dist. in canonical model
canonicalProb = @(x,beta) exp(-beta.*potentialEnergy(x));

% Number of points considered in the sample
N = 1E7;

% define global limits of the function
limGlobalInf = -5;
limGlobalSup = 5;

% define pieces of domain to sample

lims = [
    3, 5;        % interval 1
    0, 4;        % interval 2
    -1.0, 1.0;   % interval 3
    -1.5, -0.5;  % interval 4
    -2.0, -1.0;  % interval 5
    -2.5, -1.5;  % interval 6
    -3.2, -2.0;  % interval 7
    -3.5, -2.8;  % interval 8
    -4.0, -3.3;  % interval 9
    -4.5, -3.8;  % interval 10
    -4.8, -4.3;  % interval 11
    -5.0, -4.5   % interval 12
    ];

% to glue the intervals together we need to define
% global bins to sample

% define a width for counting data
boxwidth = 0.05;

% create vector for bins
binsX = limGlobalInf : boxwidth : limGlobalSup;

% get the number of bins
numBins = length(binsX) - 1;

% each interval will have a vector with sampled data and it will
% be saved into a matrix:

% get the toal number of intervals
numIntervals = size(lims, 1);

% Rows: selected interval, Cols: counts in the selected bin
allHistograms = zeros(numIntervals, numBins);

% initialize variables for Monte Carlo

% Create variable to store new position
xNew = NaN;

% define a fixed temperature
T = 5;
beta = 1/T;

% define the Monte Carlo step
deltaX = 35;


tic
for i=1:numIntervals

    % initialize data for interval i
    limInf = lims(i, 1);
    limSup = lims(i, 2);

    % define initial seed near of the upper limit
    x0 = limSup - 0.1;
    x = x0;

    % save data of the interval
    dataInterval = zeros(N, 1);
    rejectionsCount = 0;

    for j = 1:N

        % generate a new position
        while true
            xNew = x + (2*rand - 1)*deltaX;
            % break the while if the proposal is valid
            if xNew <= limSup && xNew >= limInf
                break
            end
        end

        % compute the probs. associated to old and new values
        pOld = canonicalProb(x,beta);
        pNew = canonicalProb(xNew,beta);

        if (pNew > pOld) || (rand < pNew/pOld)
            % accept the change
            x = xNew;
        else
            rejectionsCount = rejectionsCount + 1;
        end

        % save the sampled data
        dataInterval(j) = x;
    end

    % do the hist and save

    % Count data using the global bins defined
    [Hk, ~] = histcounts(dataInterval, binsX);

    % save the histogram in the global matrix
    allHistograms(i, :) = Hk;

    disp(['Interval ' num2str(i) ': ' num2str(rejectionsCount*100/N) '% Rejections'])
end
toc

%----------------------------------------------------------
% **at this point we can go to the plot section to observe
% the raw histograms per interval**
%----------------------------------------------------------


% After observing the difference in superposed histograms
% we have to fix it by adding a multiplicative factor on hists.

% numIntervals give us the total number of intervals
% we have to create the pairs of intervals to analyze

% build pairs of histograms that overlap

% define index of first histogram
idxHist_1 = (1 : numIntervals - 1)';

% index of the next histogram to compare
idxHist_2 = (2 : numIntervals)';

% unified index
histsIdx = [idxHist_1, idxHist_2];

% create a vector to store the ratio factor between histograms
ratioFactor = zeros(length(histsIdx),1);

% create a vector to store the final histogram
combinedHist = zeros(length(binsX)-1,1);

for i=1:length(idxHist_1)

    % select intervals to study
    interval1_idx = histsIdx(i,1);
    interval2_idx = histsIdx(i,2);

    % get the lims where there is overlap
    overlapMin = max(lims(interval1_idx, 1), lims(interval2_idx, 1));
    overlapMax = min(lims(interval1_idx, 2), lims(interval2_idx, 2));

    % get idx of cells overlaping
    idx_startOverlap = find(binsX <= overlapMin,1,'last');
    idx_endOverlap = find(binsX >= overlapMax, 1, 'first') - 1;

    ratio = [];

    for k = idx_startOverlap:idx_endOverlap
        % El índice del histograma es k-1 ya que histcounts devuelve un vector de N-1
        H1 = allHistograms(interval1_idx, k);
        H2 = allHistograms(interval2_idx, k);

        % Evitar división por cero
        ratio = [ratio; H1 / H2];

    end
    
    % save the mean ratio factor
    ratioFactor(i) = mean(ratio);

    % update the data using the computed ratio
    allHistograms(interval2_idx,:) = ratioFactor(i)*allHistograms(interval2_idx,:);
    
    % Save the data without overlaping in a single vector
    % notice that we do this with pairs and considering
    % the overlaping data of the secon histogram
    % this lead to rewriting of data twice per interval. One
    % at each limit of the interval. For the studied case is not
    % relevant because the data is continuous.

    idxFinH1 = find(allHistograms(interval1_idx,:),1,"last");

    idxStartH2 = find(allHistograms(interval2_idx,:),1,"first");
    idxFinH2 = find(allHistograms(interval2_idx,:),1,"last");

    combinedHist(idxStartH2:idxFinH2) = allHistograms(interval2_idx,idxStartH2:idxFinH2)';
    combinedHist(idxFinH2:idxFinH1) = allHistograms(interval1_idx,idxFinH2:idxFinH1)';
end

% at this point combinedHist stores the generated sample following
% the given potential. So we can compute for example the average
% potential by evaluationg the function on each bin and weighting
% each bin by frequency

% get the center of each bin
binCenter = (binsX(1:end-1) + binsX(2:end)) / 2;

% compute the value of potential energy per bin
potentialBin = potentialEnergy(binCenter);

% get the average potential energy for the selected temperature
meanPotential = sum(potentialBin.*combinedHist')/sum(combinedHist);

% and the average energy using equipartition theorem is given by
meanEnergy = meanPotential + 0.5*T;


%% PLOT SECTION

% we can use the following plot to observe the raw histograms
% per interval and the glue unormalized histograms 

% Plot the histograms to watch the superpositions
figure;
hold on;

% get the center of each bin
binCenter = (binsX(1:end-1) + binsX(2:end)) / 2;

for i = 1:numIntervals
    % plot the counts of each histogram
    bar(binCenter, allHistograms(i, :), 'LineWidth', 1.5, 'DisplayName', ['Interval ' num2str(i)]);
end

legend('Location', 'Best');
xlim([limGlobalInf limGlobalSup]);
hold off;

% logarithmic scale in the Y axes is reccomended for the glued
% histograms
set(gca,'YScale', 'log');

% get points from the analytical PDF

% create points to evaluate
xTeo = linspace(limGlobalInf, limGlobalSup, 10000);

% get the unormalized version
pNonNorm = canonicalProb(xTeo, beta);

% compute the normalization factor using numerical integration
Z = integral(@(x) canonicalProb(x,beta),limGlobalInf,limGlobalSup);

% get the analytical normalized PDF
pNormalized = pNonNorm / Z;

% To compare the sampled data after the treatment we have to 
% normalize per area under the curve

% get the normalization factor of data
normFactor = sum(combinedHist)*boxwidth;

% create plot with the sampled PDF
figure
bar(binCenter, combinedHist/normFactor, 'LineWidth', 1.5, 'DisplayName', ['Interval ' num2str(i)]);
hold on; 
% put the analytical pdf together
plot(xTeo, pNormalized, 'r-', 'LineWidth', 2);
set(gca,'YScale', 'log');