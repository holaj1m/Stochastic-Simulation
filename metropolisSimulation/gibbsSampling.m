% Gibbs sampling

% we want samples from P(\nu,\beta|\lambda) given:
% P(\nu | \beta) = sqrt(\beta / 2 \pi) exp(-0.5*\beta*\nu^2) (GAUSSIAN)
% P(\beta | \nu, \lambda) \propto sqrt(beta)*exp(-beta*(0.5*nu^2+\lambda))
% (GAMMA DIST).
% so in P(\beta | \nu, \lambda) we have shape k = 3/2 and 
% scale \theta = (0.5*nu^2+\lambda)^(-1)

% points to be sampled
N = 1E7;

% value of lambda given
lambda = 2;

% domain of nu
limInfNu = -15;
limSupNu =  15;

% Initialization of nu and beta
nu = 10;
beta = 30;

% Fixed parameters of distributions
kShape = 1.5;
mu = 0;

% vector to save values (nu,beta)
rawJointDist = zeros(N,2);


tic 
for i=1:N

    % 1)
    % get samples from P(\beta | \nu, \lambda)
    thetaScale = 1/(0.5*nu^2 + lambda);

    % get the beta new from gamma random generator
    beta = gamrnd(kShape,thetaScale);

    % get samples from P(\nu|\beta)
    sigma = 1/sqrt(beta);
    
    % 2)
    % get nuNew from normal random generator
    while true
        nu = normrnd(mu,sigma);
        if (nu >= limInfNu) && (nu <= limSupNu)
            break
        end
    end

    % Save generated values
    rawJointDist(i,:) = [nu, beta];

end
toc

%%
% compare the values obtained with the theoretical dist
figure
% show a histogram wit the generated data
histogram(rawJointDist(:,1),"NumBins",200,Normalization="pdf")

% theoretical curve
teoFun = @(nu,lambda) lambda.*sqrt(1./(nu.^2.+2*lambda).^3);

nuTeo = linspace(-15,15,10000);

hold on
plot(nuTeo,teoFun(nuTeo,lambda),"Color",'r','LineWidth',2)



