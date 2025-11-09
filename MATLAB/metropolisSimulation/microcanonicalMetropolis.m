% microcanonical simulation for a harmonic oscillator and plot 
% of temperature as a function of energy.


% define the potential energy
potentialEnergy = @(x) ((x+2).^2).*((x-2).^2);

% Number of points considered in the sample
N = 1E5;

% seed for initial position in the min potential
x0 = -2;
x = x0;

% vector of 0's to save the data
microData = zeros(N,1);

% counter with rejections
rejectionsCount = 0;

% Create variable to store new position
xNew = NaN;
% create varible to generate new position
delta = 3;

% total energy of the simulation
totalEnergy = linspace(1,25,100);
temps = zeros(length(totalEnergy),1);
tic
for j=1:length(totalEnergy)

    % Loop to generate data
    for i=1:N
        % generate a new position
        while true
            xNew = x + (2*rand - 1)*delta;

            % break the while if the proposal is valid
            if potentialEnergy(xNew) <= totalEnergy(j)
                break
            end
        end

        probAccep = min(1,sqrt(totalEnergy(j)-potentialEnergy(x))/sqrt(totalEnergy(j)-potentialEnergy(xNew)));

        if rand < probAccep
            x = xNew;
        else
            rejectionsCount = rejectionsCount + 1;
        end

        microData(i) = x;
    end
    temps(j) = 2*mean(totalEnergy(j)-potentialEnergy(microData));
end
toc
plot(totalEnergy,temps,'-o')
xlabel('$E$','Interpreter','latex')
ylabel('$T$','Interpreter','latex')