% Monte Carlo Metropolis simulation for an anharmonic potential
% using canonical model. The code compute the average position.

% define the potential energy
potentialEnergy = @(x) 0.5*(x-1).^3;

% define the probability of a position given beta
canonicalProb = @(x,beta) exp(-beta.*potentialEnergy(x));

% Number of points considered in the sample
N = 1E7;

% define domain of the distribution
limInf = -5;
limSup =  5;

% define the step to do mc
deltaX = 0.6;

% Create variable to store new position and data
xNew = NaN;
canonicalData = zeros(N,1);

% initialize seed
x0 = limSup; 
x = x0;

% save the first value
canonicalData(1) = x;

% fix a temperatura
T = 5;

% variable to get the average position
S = 0;

% vector to store the mean energy per temperature
meanEnergy = zeros(length(T),1);

for j=1:length(T)
    
    % get the value of inverse temperature
    beta = 1/T(j);

    % initialize rejections count
    rejectionsCount = 0;

    tic
    for i=1:N
        % generate a new position
        while true
            xNew = x + (2*rand - 1)*deltaX;
            % break the while if the proposal is valid
            if xNew <= limSup && xNew >= limInf
                break
            end
        end
        % MC Step

        % compute the probs. associated to old and new values
        pOld = canonicalProb(x,beta);
        pNew = canonicalProb(xNew,beta);

        if (pNew > pOld) || (rand < pNew/pOld)
            % accept the change
            x = xNew;
        else
            rejectionsCount = rejectionsCount + 1;
        end
        % save the position
        S = S + x;
        % save the sampled value
        canonicalData(i) = x; 
    end
    
    % get the mean energy using equipartition theorem
    meanEnergy(j) = mean(potentialEnergy(canonicalData)) + 0.5*T(j);
    toc
    % display rejections %
    disp(rejectionsCount*100/N + "% Rechazo")
end

% Compute the expected position
intMet = S/N;

% get the numerical value
fun2Int = @(x,beta) (x).*(exp(-beta.*potentialEnergy(x)));

intNum = integral(@(x) fun2Int(x,beta),-5,5)/Z;

% compare numerical integration and metropolis
errorInt = 100*abs(intNum - intMet)/abs(intNum);
disp(errorInt + "% de error en integración")

%% PLOT SECTION
% show the trace of simulation
figure; 
plot(1:N, canonicalData, 'b-');
title(['Markov Trace at T = ' num2str(T)],'interpreter','latex');
ylabel('$x$','FontSize',25,'interpreter','latex')
xlabel('MC Step','FontSize',25,'interpreter','latex')

% show a histogram wit the generated data and analytical 
% distribution

% compute analytical distribution
% define points to evaluate the pdf
xTeo = linspace(limInf, limSup, 5000);

% compute unormalized distribution
pNoNorm = canonicalProb(xTeo, beta);

% get the normalization factor of distirbution
Z = integral(@(x) canonicalProb(x,beta),limInf,limSup);

% get the normalized distribution
pNormalized = pNoNorm / Z;

% create figure with sampled data
figure
histogram(canonicalData,"NumBins",200,Normalization="pdf")
xlim([(limInf-0.125) (limSup+0.125)])
hold on; 
% create the graph of analytical pdf
plot(xTeo, pNormalized, 'r-', 'LineWidth', 2);
set(gca,'YScale', 'log');
xlabel('x','FontSize',25,'interpreter','latex');
ylabel('PDF','FontSize',25,'interpreter','latex');
set(gca, 'FontSize', 18);
legend("Simulated Data", "Analytical PDF", 'interpreter', 'latex','FontSize',16,'location', 'northeast');
xlim([limInf limSup]);
hold off;

% Create figure of E vs T
figure
plot(T,meanEnergy,'*-','Color','r')
ylabel('$\langle E \rangle$','FontSize',25,'interpreter','latex')
xlabel('$T$','FontSize',25,'interpreter','latex')
set(gca, 'FontSize', 18);

