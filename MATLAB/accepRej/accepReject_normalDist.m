tic
% Number of points considered in the sample
N = 1E7;

% select parameters for the normal distribution
mu = 0;
sigma = 1;

% parameter related to the max probability of the dist
xMax = 0;
pMax = normpdf(xMax,mu,sigma);

% vector of 0's to save the data
normalData = zeros(N,1);

% counter with accepted data
acceptedCount = 1;

while acceptedCount < N+1
    % generate rand number within the limits of the dist
    xProposal    = -5 + 2*5*rand;
    probProposal = pMax*rand;

    % acceptance-rejection condition
    if probProposal < normpdf(xProposal,mu,sigma)

        normalData(acceptedCount) = xProposal;

        % increase the count
        acceptedCount = acceptedCount+ 1;
    end
end
toc

%% Plot Section
% show a histogram wit the generated data
histogram(normalData,"NumBins",200,Normalization="pdf")

% Compare the generated data with the analytical solution
xData = linspace(-5,5,200);
hold on
plot(xData,normpdf(xData,mu,sigma),'r-', 'LineWidth', 2)
yline(pMax,'r--', 'LineWidth', 2)
hold off

ylabel('$P(x|I)$','FontSize',20,'interpreter','latex')
xlabel('$x$','FontSize',20,'interpreter','latex')


