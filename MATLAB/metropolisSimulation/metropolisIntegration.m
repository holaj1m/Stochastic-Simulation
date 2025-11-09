% Metropolis simulation for computing the value of a 1D integral.

% Number of points considered in the sample
N = 1E8;

% We want to use MC to integrate a function, thus we have 
% F(x) = P(x)A(x) = log(1+x)/sqrt(x)
intFun = @(x) log(1+x)./sqrt(x);
% the function chosen as P(x) must be normalized

% integration limits
xMin = 2; xMax = 10; 

% define P(x) normalized
probFun = @(x) log(1+x)/(-8 + 11*log(11) - log(27));

% verify that probFun is correctly normalized
integral(probFun,xMin,xMax);

% define A(x). The normalization factor must be added 
% in a way that F(x) is unchanged
expectFun = @(x) (-8 + 11*log(11) - log(27))/sqrt(x);

% seed for initial position and position
x0 = 6;
x = x0;

% vector of 0's to save the data
probData = zeros(N,1);
probData(1) = x;

% counter with rejections
rejectionsCount = 1;

% Create variable to store new position
xNew = NaN;
% variable to get the mean
S = 0;
% Loop to generate data 
tic
for i=2:N
    
    % generate a proposal of position 
    while true
        % uniform distribution wirhin the limts of dist
        xNew = x + (2*rand -1)*30;
        if (xNew < 10) && (xNew > 2)
            % break the while if P(xNew) exist
            break;
        end
    end

    probAccep = min(1,probFun(xNew)/probFun(x));

    if rand < probAccep
        x = xNew;
    else
        rejectionsCount = rejectionsCount + 1;
    end

    probData(i) = x;
    S = S + expectFun(x);
end
toc

disp(rejectionsCount*100/N + "% Rechazo")

%% Plot Section

% show a histogram wit the generated data
histogram(probData,"NumBins",200,Normalization="pdf")

% Compare the generated data with the analytical solution
xData = linspace(xMin,xMax,200);
hold on
plot(xData,probFun(xData),'r-', 'LineWidth', 2)
hold off

% Compute the integral
intMet = S/(N-1);
% get the numerical value
intNum = integral(intFun,xMin,xMax);

%compare them
errorInt = 100*abs(intNum - intMet)/abs(intNum);
disp(errorInt + "% de error en integración")

