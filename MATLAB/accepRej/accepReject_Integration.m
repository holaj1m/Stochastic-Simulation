% The idea of this code is to generate data from 
% P(x|I) = 1/Z exp(-0.5*x^2) using accept-reject method and 
% compute the expectation of <x^2> using the sampled values.

% set the total number of points to generate
N = 1E7;

% define the given probability distribution
distExp = @(x) exp(-0.5*x.^4);

% get the normalization factor
normFact = integral(distExp,-2,2);

% the max value of the pdf is given by
pMax = 1/normFact;

% definevector of zeros to save data
genData = zeros(N,1);

% counter with accepted data. I start form 1 because otherwise 
% i can't save the first value generated
acceptedCount = 1;

% define function to integrate
expectFun = @(x) x.^2;
% define variable to get the integral by limits theorem 
s = 0;
tic
while acceptedCount < N+1
    % generate a rand number within the limits of dist.
    xProposal = -2 + 4*rand;
    probProposal = pMax*rand;

    % acceptance-rejection condition
    if probProposal < distExp(xProposal)/normFact
        genData(acceptedCount) = xProposal;

        % increse the count
        acceptedCount = acceptedCount + 1;

        s = s + expectFun(xProposal);
    end
end
toc
% the value of the integral is obtained from s
aproxInt = s/N;
disp('La integral es '+ aproxInt)

% compare this value with the numerical integration
integrateFun = @(x) x.^2 .* exp(-0.5*x.^4)/normFact;
numValInt = integral(integrateFun,-2,2);

% If we set the numerical value as reference the relative
% error is given by:
errorInt = 100*abs(numValInt - aproxInt)/abs(numValInt);
disp(errorInt + "% de error en integración")
%% Plots comparing analytical pdf and the generated data

% show a histogram wit the generated data
histogram(genData,"NumBins",200,Normalization="pdf")

%compare with the analytical distribution
xData = linspace(-2,2,200);
dataDist = distExp(xData)/normFact;
hold on
plot(xData,dataDist,'r-', 'LineWidth', 2)
yline(pMax,'r--', 'LineWidth', 2)
hold off