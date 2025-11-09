% Code for a simulation of a one-dimensional Ising chain
% using montecarlo metropolis.

% set the couple constant
j = 1;
% number of states to consider
N = 100;
% create the array of states
s = ones(N,1);

% counter with rejections
rejectionsCount = 1;
% total steps simulation
nSteps = 1E4;
% total temperatures to visit
nTemps = 100;

% vectors to save energy and rejection rate for each T
energy = zeros(nSteps,1);
rejectionsRate = zeros(nSteps,1);

% create the array of temperatures
T = linspace(0.01,10,nTemps);

% vectors with final data for each T
meanE = zeros(nTemps,1);
stdE = zeros(nTemps,1);
meanRej = zeros(nTemps,1);
tic
for m=1:nTemps

    % inverse temperature
    beta = 1/T(m);

    % simulate nSteps with a fixed beta
    for l=1:nSteps

        % We do a MC sweep over the system
        for i=1:N
            % compute the current energy of sys
            energyOld = Energy(s, N, j);

            % change a random states in states
            k = randi(N);
            s(k) = -1*s(k);

            % compute the energy considering the change
            energyNew = Energy(s, N, j);

            if (energyNew < energyOld) || (rand < exp(-beta*(energyNew-energyOld)))
                % if the energy decrese or the probability increase accept the
                % change
                continue
            else
                % in the opposite case undo the change in the state
                rejectionsCount = rejectionsCount+1;
                s(k) = -1*s(k);
            end
        end

        energy(l) = Energy(s,N,j);
        rejectionsRate(l) = 100*rejectionsCount/N;

        % set rejectrions count to 0 again
        rejectionsCount = 0;

    end

    meanE(m) = mean(energy);
    stdE(m) = std(energy);
    meanRej(m) = mean(rejectionsRate);

end
toc
clear j N s rejectionsCount nSteps nTemps energy rejectionsRate
clear beta l m k i energyOld energyNew

errorbar(T,meanE,stdE)

%%
function e = Energy(states,numSites,j)
    e = 0;
    % compute the energy avoiding the las element in the chain
    for i=1:(numSites-1)
        e = e - j*states(i)*states(i+1);
    end
    % compute the energy for the las element in the chain
    e = e - j*states(1)*states(numSites);
end
