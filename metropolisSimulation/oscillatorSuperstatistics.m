% Simulation of harmonic oscillator in superstatistics model.
% We assume that temperature follow a gamma distribution.

% number of points considered in the sample
N = 1E6;

% def. the coupling constant
k = 1;

% distributions parameters
deltaX = 4;
deltaBeta = 1;
bK = 2.5;
bTheta = 1/bK;

% energy of the oscilator
Energy = @(x) 0.5*k*x^2;

% joint prob. function
% here we assume that Temps. follow a gamma distribution
jointProb = @(x,beta) exp(-beta*Energy(x))*(beta^(bK-1))*exp(-beta/bTheta);

% seed for initial position and inverse temp
x0 = 0;
beta = 1;
x = x0;

% vector of 0's to save the data
superData = zeros(N,1);
superData(1) = x;

% counter with rejections
rejectionsCount = 0;

% Create variable to store new position and inverse temp
xNew = NaN;
betaNew = NaN;

% Loop to generate data
tic
for i=2:N
    % generate a proposal of position
    xNew = x + (2*rand - 1)*deltaX;
    % generate a proposal of inverse temp.
    while true
        betaNew = beta + (2*rand - 1)*deltaBeta;
        % break the while if betaNew can be evaluated
        if betaNew > 0
            break;
        end
    end

    % compute the probs. associated to old and new values
    pOld = jointProb(x,beta);
    pNew = jointProb(xNew,betaNew);

    if (pNew > pOld) || (rand < pNew/pOld)
        % accept the change
        x = xNew;
        beta = betaNew;
    else
        rejectionsCount = rejectionsCount + 1;
    end

    superData(i) = x;
end
toc
disp(rejectionsCount*100/N + "% Rechazo")

%% Plot Section
% show a histogram wit the generated data
histogram(superData,"NumBins",200,Normalization="pdf")
xlim([-6.5 6.5])
% Compare the generated data with the equilibrium solution
xData = linspace(-6.5,6.5,200);
hold on
plot(xData,normpdf(xData,mean(superData),std(superData)),'r-', 'LineWidth', 2)
hold off