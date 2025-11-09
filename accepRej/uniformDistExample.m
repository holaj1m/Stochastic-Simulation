% Number of points considered in the sample
N = 1E7;
% create the sample of N random numbers arreanged in a column
data = rand(N,1);

% show a histogram with data
histogram(data,"NumBins",200,Normalization="pdf")