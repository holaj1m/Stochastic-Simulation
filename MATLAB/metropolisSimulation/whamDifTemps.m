%% weighted histogram analysis method (WHAM)
% The idea of this code is to simulate data for a particle
% under the anharmonic potential \phi(x) = -0.5*(x-1)^3 at
% different temperatures T between 3 and 10. 
% We use the canonical model and Métropolis method to sampling
% the position of particle. However the nature of potential
% energy makes very difficult to sample the space properly. 
% So to fix this issue we divided the domain in different 
% intervals. Given that each interval has its own distribution
% we have to plug the intervals together using different factor.
% The method described above corresponds to WHAM.

% define the potential energy
potentialEnergy = @(x) 0.5*(x-1).^3;

% define the probability of a position given beta
canonicalProb = @(x,beta) exp(-beta.*potentialEnergy(x));

% Number of points considered in the sample
N = 1E7;

% 1. define limits of the function
limGlobalInf = -5; 
limGlobalSup = 5; 

% define pieces of domain to sample

lims = [
    3, 5;       % interval 0
    0, 4;       % interval 1 
   -1.0, 1.0;   % interval 2
   -1.5, -0.5;  % interval 3
   -2.0, -1.0;  % interval 4
   -2.5, -1.5;  % interval 5
   -3.2, -2.0;  % interval 6
   -3.5, -2.8;  % interval 7
   -4.0, -3.3;  % interval 8 
   -4.5, -3.8;  % interval 9
   -4.8, -4.3;  % interval 10
   -5.0, -4.5   % interval 11
];

% define the step to do mc
deltaX = 20;

% define a width for counting data
boxwidth = 0.05;

% create vector for bins
binsX = limGlobalInf : boxwidth : limGlobalSup;

% number of bins
numBins = length(binsX) - 1;

% get the center of each bin
binCenter = (binsX(1:end-1) + binsX(2:end)) / 2;
% compute the value of potential energy per bin
V_bin = potentialEnergy(binCenter);

% Matrix to save data

% Rows: selected window, Cols: hist cell
allHistograms = zeros(size(lims, 1), numBins); 

numIntervals = size(lims, 1);

% define index of first histogram
idxHist_1 = (1 : numIntervals - 1)';

% index of the next histogram to compare
idxHist_2 = (2 : numIntervals)';

% unified index
histsIdx = [idxHist_1, idxHist_2];

% create a vector to store the ratio factor
ratioFactor = zeros(length(histsIdx),1);



% Create variable to store new position and inverse temp
xNew = NaN;

% create array of temperatures
T = 3:10;

meanEnergy = zeros(length(T),1);

for m=1:length(T)

    % fix a temperature
    beta = 1/T(m);
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

    % After observing the difference in superposed histograms
    % we have to fix it by adding a multiplicative factor on hists.

    % numIntervals give us the total number of intervals
    % we have to create the pairs of intervals to analyze

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

        ratioFactor(i) = mean(ratio);
        allHistograms(interval2_idx,:) = ratioFactor(i)*allHistograms(interval2_idx,:);

        idxStartH1 = find(allHistograms(interval1_idx,:),1,"first");
        idxFinH1 = find(allHistograms(interval1_idx,:),1,"last");

        idxStartH2 = find(allHistograms(interval2_idx,:),1,"first");
        idxFinH2 = find(allHistograms(interval2_idx,:),1,"last");

        combinedHist(idxStartH2:idxFinH2) = allHistograms(interval2_idx,idxStartH2:idxFinH2)';
        combinedHist(idxFinH2:idxFinH1) = allHistograms(interval1_idx,idxFinH2:idxFinH1)';


    end


    

    % get the mean potential energy
    meanPotential_T = sum(V_bin .* combinedHist')/sum(combinedHist);

    % get the mean energy by equipartition theorem
    meanEnergy(m) = meanPotential_T + 0.5*T(m);

end

%%

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

set(gca,'YScale', 'log');


normFactor = sum(combinedHist)*boxwidth;
figure
bar(binCenter, combinedHist/normFactor, 'LineWidth', 1.5, 'DisplayName', ['Interval ' num2str(i)]);
set(gca,'YScale', 'log');

save("dataEvsT_EJ2T2.mat","meanEnergy","T")
plot(T,meanEnergy,'*-')


