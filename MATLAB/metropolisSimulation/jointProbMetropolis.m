% Monte Carlo Metropolis for conditional probability

% The idea is getting the joint probability of P(\nu,\beta|\lambda)
% We know tha P(\beta | \nu, \lambda) \propto P(\nu|\beta)*P(\beta|\lambda)
% given.

% We have P(\beta|\lambda) = \lambda exp(-\lambda \beta)
% P(\nu|\beta) = sqrt(\beta / 2 \pi) exp(-0.5*\beta*\nu^2)
% So we can sample directly P(\beta, \nu | \lambda) by a 
% two step metropolis

% number of points to be sampled
N = 1E7;

% define the dist. of \nu
probNu = @(beta,nu) exp(-0.5*beta*nu^2);

% define P(\beta | \nu, \lambda)
probBetaNu = @(beta, nu, lambda) sqrt(beta)*exp(-beta*(0.5*nu^2+lambda));

% lambda is given by the problem
lambda = 2;

% the domain of \nu and beta are given by the problem
limInfNu = -15;
limSupNu =  15;

limInfBeta = 0;

% define delta for nu and beta to propose new values
deltaNu = 8.0;
deltaBeta = 15.0;

% initialize beta and nu
beta = 30;
nu = 10;

betaNew = NaN;
nuNew = NaN;

% define the joint distribution (nu,beta)
rawJointDist = zeros(N,2);

% initialize rejection count
rejectionCountNu = 0;
rejectionCountBeta = 0;

tic
for i=1:N

    % 1. Do metropolis to get a value of nu

    % propose a new value of nu
    while true
        nuNew = nu + (2*rand - 1)*deltaNu;
        if nuNew <= limSupNu && nuNew >= limInfNu
            break
        end
    end

    % check if the value is accepted or rejected
    % using the current value of beta

    probOldNu = probNu(beta,nu);
    probNewNu = probNu(beta,nuNew);

    if (probNewNu > probOldNu) || (rand < probNewNu/probOldNu)
        % accept the new value
        nu = nuNew;
    else
        rejectionCountNu = rejectionCountNu + 1;
    end

    % 2. now with the value of nu sampled we will get the 
    % new value of beta from de conditional prob

    % propose a value for beta
    while true
        betaNew = beta + (2*rand - 1)*deltaBeta;
        if betaNew >= limInfBeta
            break
        end
    end

    % check if the value is accepted or rejected
    probOldBeta = probBetaNu(beta,nu,lambda);
    probNewBeta = probBetaNu(betaNew,nu,lambda);

    if (probNewBeta > probOldBeta) || (rand < probNewBeta/probOldBeta)
        beta = betaNew;
    else
        rejectionCountBeta = rejectionCountBeta + 1;
    end
    
    % save the generated values
    rawJointDist(i,:) = [nu, beta];

end
toc
disp(rejectionCountNu*100/N + "% Rechazo Nu")
disp(rejectionCountBeta*100/N + "% Rechazo Beta")

%% Plot Section

% Here we compare the marginalized data with the analytical
% function.

figure
% show a histogram wit the generated data
histogram(rawJointDist(:,2),"NumBins",200,Normalization="pdf")
set(gca,'YScale', 'log');

% Compare with the theoretical curve 
teoFun = @(nu,lambda) lambda.*sqrt(1./(nu.^2.+2*lambda).^3);

nuTeo = linspace(-15,15,10000);

hold on
plot(nuTeo,teoFun(nuTeo,lambda),"Color",'r')



