% Code to generate data from exponential distirbution
% using inverse transform sampling

% Number of points considered in the sample
N = 1E7;

% define parameter exp distribution
lambda = 3;

% vector with uniform distribution
u = rand(N,1);

% generated data with the inverse CDF
distData = inverseCDF(u,lambda);

% show a histogram with data
histogram(distData,"NumBins",200,Normalization="pdf")

% compare the data with the distribution
xData = linspace(0,2.5,200);
hold on
plot(xData,lambda*exp(-lambda.*xData),'r-', 'LineWidth', 2)
hold off
xlim([0 2.5])

%%
function distData = inverseCDF(u, param)
% define the inverse of cdf to evaluate uniform dist
    distData = -(1/param)*log(1-u);
end