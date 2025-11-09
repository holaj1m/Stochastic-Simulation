tic
% Number of points considered in the sample

% the program runs with N up to 10^6 but the waiting time is longer

N = 1E5;

% select parameters for the gamma distribution
kShape = 2.5;
thetaScale = 5;

% parameter related to the max probability of the dist
xMax = (kShape-1)*thetaScale;
pMax = gampdf(xMax,kShape,thetaScale);

% vector of 0's to save the data
gammaData = zeros(N,1);

% counter with accepted data
acceptedCount = 1;

while acceptedCount < N+1
    % generate rand number within the limits of the dist
    xProposal    = 50*rand;
    probProposal = pMax*rand;

    % acceptance-rejection condition
    if probProposal < gampdf(xProposal,kShape,thetaScale)

        gammaData(acceptedCount) = xProposal;

        % increase the count
        acceptedCount = acceptedCount+ 1;
    end
end
toc

%% Plot section

% show a histogram wit the generated data
histogram(gammaData,"NumBins",200,Normalization="pdf")

% Compare the generated data with the analytical solution
xData = linspace(0,50,200);
hold on
plot(xData,gampdf(xData,kShape,thetaScale),'r-', 'LineWidth', 2)
yline(pMax,'r--', 'LineWidth', 2)
hold off
