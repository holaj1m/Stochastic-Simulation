%% Simulated Annealing
% Simulated annealing optimization for the Rosenbrock function:
% f(x) = \sum_{i=1}^{n-1} \left[ 100(x_{i+1} - x_i^2)^2 + (1 - x_i)^2 \right]

% define dimension number
dims = 10;

% create a random initial state with values between -1 and 1
state = 2*rand(1,dims)-1;

% define the Rosenbrock function
evalFunction = @(s) sum(100*(s(2:end)-s(1:end-1).^2).^2 + ...
                            (1-s(1:end-1)).^2);

% initialize the montecarlo step
step = 0.5;

% initialize the temperature
initTemp = 20;
temp = initTemp;

% initialize best state and value
bestState = state;
bestVal = evalFunction(state);

% initialize global counter and rejection counter
% notice that global counter also gives the number of
% calls to the function
rejectionsCount = 0;
iterCount = 0;

% initialize variable to store difference between values
% obtained in each iteration
minDif = 100;

% initialize rejection rate
rejectRate = 1;

% set the threshold to adapt step over simulation
maxAdapt = 0.2; 


while (minDif > 1E-4) && (temp > 1E-12)

    % get the value of the function in the current state
    currentValue = evalFunction(state);

    % propose a new state
    newState = state + step.*(2*rand(1,dims)-1);
    
    % get the value 0^-12of the function in the new state
    proposedVal = evalFunction(newState);

    % do the MC update
    % notice that we verify if the new value minimize the function
    if (proposedVal < currentValue) || (rand < exp(-(proposedVal-currentValue)/temp))
        % update the state
        state = newState;

        % check if the new value is better than the best seen so far
        if proposedVal < bestVal
            minDif = bestVal-proposedVal;
            bestVal = proposedVal;
            bestState = state;
        end
    else
        rejectionsCount = rejectionsCount + 1;
    end

    % increse the global counter
    iterCount = iterCount + 1;
    
    % compute the rejections rate
    rejectRate = 100*rejectionsCount/iterCount;

    % cold the system
    temp = temp*0.999999;

    % adapt the step depending on temperature
    alpha = maxAdapt * (temp / initTemp);
    if rejectRate > 75
        step = step*(1 - alpha);
    elseif rejectRate < 60
        step = step*(1 + alpha);
    end


    % display min val, temperature and rejections rate
    % every 1E6 iterations
    if mod(iterCount,1E6)==0
        disp(bestVal+" "+temp+" "+rejectRate)
    end


end

% clear variables
clear state rejectRate rejectionsCount newState proposedVal 
clear maxAdapt currentValue 
