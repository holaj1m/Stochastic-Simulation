% Number of points considered in the sample
N = 1E7;

% select parameters for the normal distribution
mu = 0;
sigma = 1;

% create a vector sample of N random numbers 
normalData = normrnd(mu,sigma,N,1);

% show a histogram with data
histogram(normalData,"NumBins",200,Normalization="pdf")

% compare the generated data with the actual gaussian function
xData = linspace(-4.1,4.1,200);
yData = Gauss(xData,mu,sigma);

% Plot the normal pdf
hold on
plot(xData, yData, 'r-', 'LineWidth', 2);
hold off
%%
function normVar = Gauss(x,mu,sigma)
    normVar = exp(-0.5/(sigma^2) * (x-mu).^2)/(sigma*sqrt(2*pi));
end
%%
% We can also use the normPdf to avoid the function def

% Number of points considered in the sample
N = 1E7;

% select parameters for the normal distribution
mu = 0;
sigma = 1;

% create a vector sample of N random numbers 
normalData = normrnd(mu,sigma,N,1);

% show a histogram with data
histogram(normalData,"NumBins",200,Normalization="pdf")

% compare the generated data with the actual gaussian function
xData = linspace(-4.1,4.1,200);

% plot the normal pdf
hold on
plot(xData, normpdf(xData, mu, sigma), 'r-', 'LineWidth', 2);
hold off