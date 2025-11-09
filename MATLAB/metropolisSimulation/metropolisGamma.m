% Monte Carlo simulation to sampling the gamma distribution.

% Number of points considered in the sample
N = 5E5;

% select parameters for the gamma distribution
kShape = 2.5;
thetaScale = 5;

% position where dist is max
xMax = (kShape-1)*thetaScale;

% seed for initial position and position
x0 = 30;
x = x0;

% vector of 0's to save the data
gammaData = zeros(N,1);
gammaData(1) = x;

% counter with rejections
rejectionsCount = 0;

% Create variable to store new position
xNew = NaN;

% Loop to generate data
tic
for i=2:N
    
    % generate a proposal of position 
    while true
        % uniform distribution wirhin the limts of dist
        xNew = x + (2*rand -1)*50;
        if (xNew < 50) && (xNew > 0)
            % break the while if P(xNew) exist
            break;
        end
    end

    probAccep = min(1,gampdf(xNew,kShape, thetaScale)/gampdf(x,kShape, thetaScale));

    if rand < probAccep
        x = xNew;
    else
        rejectionsCount = rejectionsCount + 1;
    end

    gammaData(i) = x;
end
toc
disp(rejectionsCount*100/N + "% Rechazo")

%% Plot Section

% show a histogram wit the generated data
histogram(gammaData,"NumBins",200,Normalization="pdf")

% Compare the generated data with the analytical solution
xData = linspace(0,50,200);
hold on
plot(xData,gampdf(xData,kShape,thetaScale),'r-', 'LineWidth', 2)
hold off

% trace plot
figure 
plot(gammaData,'-')