clear; clc; close all;
%% Exercise 3: Sequential Importance Sampling
disp("Sequential Importance Sampling")
estimate = zeros(1, 25);

for N = 2:26
    counter = 0;
    
    for size = 1:1e4
        X = zeros(1, N);
        Y = zeros(1, N);
        
        for i = 2:N
            r = randi([1, 4]);
            switch r 
                case 1 % Move one step to the right
                    X(i) = X(i - 1) + 1;
                    Y(i) = Y(i - 1);
                case 2 % Move one step to the left
                    X(i) = X(i - 1) - 1;
                    Y(i) = Y(i - 1);
                case 3 % Move one step up
                    Y(i) = Y(i - 1) + 1;
                    X(i) = X(i - 1);
                case 4 % Move one step down
                    Y(i) = Y(i - 1) - 1;
                    X(i) = X(i - 1);
            end  
        end

        selfAW = true; % Assume SAW is TRUE

        for i = 1:N-1
           for j = i+1:N
               if (X(i) == X(j) && Y(i) == Y(j))
                  selfAW = false; % SAW is FALSE
                  break
               end
           end

           if (~selfAW)
               counter = counter + 1;
               break
           end
        end
    end

    ratio = 1 - counter / size;
    steps = N - 1;
    estimate(steps) = 4^steps * ratio;
    
    if steps == 1
        disp("Possible SAW for " + steps + " step is " + estimate(steps))
    else
        disp("Possible SAW for " + steps + " steps is " + estimate(steps))
    end
end

constant = zeros(1, length(estimate) - 1);
gamma = 11 / 32;

for n = 1:length(constant)
    constant(n) = estimate(n+1) / estimate(n) * ((n / (n+1))^gamma);
end

green = [27, 158, 119]./255;
orange = [217, 95, 2]./255;
avgConstant = mean(constant);
showAvgConstant = zeros(1, length(constant)) + avgConstant;
figure()
p1 = plot(constant);
hold on;
p2 = plot(showAvgConstant);
hold off;
xlabel('Steps');
ylabel('Connectivity constant');
title('Estimating the Connectivity Constant with SIS');
legend('Step estimate', 'Average');
p1.Color = green;
p2.Color = orange;
%% Exercise 4: Improved Sequential Importance Sampling
disp("-----------------------------------------------------------------")
disp("Improved Sequential Importance Sampling")
size = 1e3; 
d = 2; % Dimension size
steps = 50; % Maximum number of stepsis, is changed to 100 for figure 4
XY = zeros(steps, d); 

for i = 2:size
    XY = cat(3, XY, zeros(steps, d)); % Stores data for all steps size times
end

w = zeros(steps + 1, size); % Importance weight function
w(1, :) = 1; % We set the first weight to one

for N = 2:steps + 1 
    for i = 1:size 
        z = 0; % Indicator function, initial value 0, changed to 1 if it is possible to take the next step
        preSteps = XY(1:N - 1, :, i); % Nxd matrix of previous steps
        currentPos = XY(N - 1, :, i); % Coordinates of last step
        possibleSteps = zeros(2, d);
        
        for dim = 1:d
            nextStep = zeros(1, d);
            nextStep(dim) = 1;
            possibleSteps(1, dim) = 1 - ismember(currentPos + nextStep, preSteps, 'rows'); % value at (1,dim) set to 1 if possible to walk, else 0
            possibleSteps(2, dim) = 1 - ismember(currentPos - nextStep, preSteps, 'rows');% value at (2,dim) set to 1 if possible to walk, else 0
        end
        
        nbrOfWalks = sum(sum(possibleSteps)); % Total number of possible walks
        
        if(nbrOfWalks == 0) 
            XY(N, :, i) = XY(N - 1, :, i); % If there are no possible walks, set X_k+1 = X_k
        else
            z = 1; % Set indicator function to 1 if it is possible to take one step
            g = 1 / nbrOfWalks; % Probability of our instrumental distribution
            [r, c] = find(possibleSteps, nbrOfWalks); % Returns index for row and cols for values that are nonzero (1)
            rand = randi(nbrOfWalks); % Generate uniform 
            row = r(rand); % row = 1 => positive step, row = 2 => negative step
            col = c(rand); % Which dimension to walk in [1, 2,..., d]
            temp = zeros(1, d);
            temp(col) = 1;
            
            if(row == 1)
                XY(N, :, i) = XY(N - 1, :, i) + temp; % Walk positive step in dimension temp
            else 
                XY(N, :, i) = XY(N - 1, :, i) - temp; % Walk negative step in dimension temp
            end
        end
        
        w(N, i) = (z / g) * w(N - 1, i); % Update the weighting function
    end
end

c_n = zeros(N, 1); 

for n = 1:N
    c_n(n) = mean(w(n, :)); % Estimate of c_n
end

c_n = c_n(2:end); % Removes the first value

for k = 1:length(c_n)
    if k == 1
        disp("Estimation of c_n for " + k + " step: " + c_n(k));
    else
        disp("Estimation of c_n for " + k + " steps: " + c_n(k));
    end
end

constant = zeros(1, length(c_n) - 1);
gamma = 11 / 32;

for n = 1:length(constant)
    constant(n) = c_n(n+1) / c_n(n) * ((n / (n+1))^gamma);
end

avgConstant = mean(constant);
showAvgConstant = zeros(1, length(constant)) + avgConstant;
figure()
p1 = plot(constant);
hold on;
p2 = plot(showAvgConstant);
xlabel('Steps');
ylabel('Connectivity constant');
title('Estimating the Connectivity Constant with the Improved SIS');
legend('Step estimate', 'Average');
p1.Color = green;
p2.Color = orange;
%% Exercise 5: Sequential Importance Sampling with Resampling
disp("-----------------------------------------------------------------")
disp("Sequential Importance Sampling with Resampling")
size = 1e3;
d = 2; % Dimension size
steps = 50; % Maximum number of steps, is changed to 100 for figure 4
XY = zeros(steps, d);

for i = 2:size
    XY = cat(3, XY, zeros(steps, d)); % Stores data for all steps size times
end

w = zeros(steps + 1, size); % Importance weight function
w(1, :) = 1; % We set the first weight to one
z_2 = zeros(size, 1);
z_1 = ones(size, 1);

for N = 2:steps + 1 
    for i = 1:size
        preSteps = XY(1:N - 1, :, i); % Nxd matrix of previous steps
        currentPos = XY(N - 1, :, i); % Coordinates of last step
        possibleSteps = zeros(2, d);
        
        for dim = 1:d
            nextStep = zeros(1, d);
            nextStep(dim) = 1;
            possibleSteps(1, dim) = 1 - ismember(currentPos + nextStep, preSteps, 'rows'); % value at (1,dim) set to 1 if possible to walk, else 0
            possibleSteps(2, dim) = 1 - ismember(currentPos - nextStep, preSteps, 'rows');% value at (2,dim) set to 1 if possible to walk, else 0
        end
        
        nbrOfWalks = sum(sum(possibleSteps)); % Total number of possible walks

        if(nbrOfWalks == 0) 
            XY(N, :, i) = XY(N - 1, :, i); % If there is no possible walks, set X_k+1 = X_k
            z_2(i) = 0;
        else
            z_2(i) = 1; % Set indicator function to 1 if it is possible to take one step
            g = 1 / nbrOfWalks; % Probability of our instrumental distribution
            [r, c] = find(possibleSteps, nbrOfWalks); % Returns index for row and cols for values that are nonzero (1)
            rand = randi(nbrOfWalks); % Generate uniform 
            row = r(rand); % row = 1 => positive step, row = 2 => negative step
            col = c(rand); % Which dimension to walk in [1, 2,..., d]
            temp = zeros(1, d);
            temp(col) = 1;
            
            if(row == 1)
                XY(N, :, i) = XY(N - 1, :, i) + temp; % Walk positive step in dimension temp
            else 
                XY(N, :, i) = XY(N - 1, :, i) - temp; % Walk negative step in dimension temp
            end
        end
        
        w(N, i) = z_2(i) / (z_1(i) * g); % Update the weighting function
    end
    
    ind = randsample(size, size, true, w(N, :));
    XY = XY(:, :, ind);
    z_1 = z_2(ind);
end

c_nMod = zeros(steps, 1); 
c_nMod(1) = sum(w(2,:))/size;

for s = 3:steps + 1
    c_nMod(s-1) = (sum(w(s, :))/size)*c_nMod(s-2);
end

for k = 1:length(c_nMod)
    if k == 1
        disp("Estimation of c_n for " + k + " step: " + c_nMod(k));
    else
        disp("Estimation of c_n for " + k + " steps: " + c_nMod(k));
    end
end

constantMod = zeros(1, length(c_nMod) - 1);
gamma = 11/32;

for n = 1:length(constantMod)
    constantMod(n) = c_nMod(n+1) / c_nMod(n) * ((n / (n+1))^gamma);
end

avgConstant = mean(constantMod);
showAvgConstant = zeros(1, length(constantMod)) + avgConstant;
figure()
p1 = plot(constantMod);
hold on;
p2 = plot(showAvgConstant);
hold off;
xlabel('Steps');
ylabel('Connectivity constant');
title('Estimating the Connectivity Constant with SISR');
legend('Step estimate', 'Average');
p1.Color = green;
p2.Color = orange;

% Figure 4
figure()
p1 = plot(constant);
hold on;
p2 = plot(constantMod);
hold off;
xlabel('Steps');
ylabel('Connectivity constant');
title('Comparing the Connectivity Constant for SIS and SISR');
legend('SIS', 'SISR');
p1.Color = green;
p2.Color = orange;
%% Exercise 6: Obtaining Estimates
disp("-----------------------------------------------------------------")
disp("Obtaining Estimates")
a2 = zeros(1, 10);
mu2 = zeros(1, 10);
gamma2 = zeros(1, 10);
size = 1e3; 
d = 2; % Dimension size
steps = 100; % Maximum number of steps

for iterations = 1:10
    XY = zeros(steps, d);
    
    for i = 2:size
        XY = cat(3, XY, zeros(steps, d)); % Stores data for all steps size times
    end
    
    w = zeros(steps + 1, size); % Importance weight function
    w(1, :) = 1; % We set the first weight to one
    z_2 = zeros(size, 1);
    z_1 = ones(size, 1);
    
    for N = 2:steps + 1 
        for i = 1:size
            preSteps = XY(1:N - 1, :, i); % Nxd matrix of previous steps
            currentPos = XY(N - 1, :, i); % Coordinates of last step
            possibleSteps = zeros(2, d);
            
            for dim = 1:d
                nextStep = zeros(1, d);
                nextStep(dim) = 1;
                possibleSteps(1, dim) = 1 - ismember(currentPos + nextStep, preSteps, 'rows'); % value at (1,dim) set to 1 if possible to walk, else 0
                possibleSteps(2, dim) = 1 - ismember(currentPos - nextStep, preSteps, 'rows');% value at (2,dim) set to 1 if possible to walk, else 0
            end
            
            nbrOfWalks = sum(sum(possibleSteps)); % Total number of possible walks

            if (nbrOfWalks == 0) 
                XY(N, :, i) = XY(N - 1, :, i); % If there is no possible walks, set X_k+1 = X_k
                z_2(i) = 0;
                w(N, i) = 0; 
            else
                z_2(i) = 1; % Set indicator function to 1 if it is possible to take one step
                g = 1 / nbrOfWalks; % Probability of our instrumental distribution
                [r, c] = find(possibleSteps, nbrOfWalks); % Returns index for row and cols for values that are nonzero (1)
                rand = randi(nbrOfWalks); % Generate uniform 
                row = r(rand); % row = 1 => positive step, row = 2 => negative step
                col = c(rand); % Which dimension to walk in [1, 2,..., d]
                temp = zeros(1, d);
                temp(col) = 1;
                w(N, i) = z_2(i) / (z_1(i) * g);
                
                if(row == 1)
                    XY(N, :, i) = XY(N - 1, :, i) + temp; % Walk positive step in dimension temp
                else 
                    XY(N, :, i) = XY(N - 1, :, i) - temp; % Walk negative step in dimension temp
                end
            end
        end
        
        ind = randsample(size, size, true, w(N, :));
        XY = XY(:, :, ind);
        z_1 = z_2(ind);
    end
    
    c_nEst = zeros(steps, 1); 
    c_nEst(1) = sum(w(2,:))/size;
    
    for s = 3:steps + 1
        c_nEst(s-1) = (sum(w(s, :))/size)*c_nEst(s-2);
    end
    
    q = 1:1:length(c_nEst);
    y = @(q) log(c_nEst(q));
    X1 = ones(length(c_nEst), 1);
    X2 = @(q) transpose(q * 1);
    X3 = @(q) transpose(1 * log(q));
    X = [X1(q), X2(q), X3(q)];
    b = regress(y(q), X);
    a2(iterations) = exp(b(1));
    mu2(iterations) = exp(b(2));
    gamma2(iterations) = 1 + b(3);
end

% Displaying the arrays containing the parameter estimates
disp(a2)
disp(mu2)
disp(gamma2)

% Displaying the parameters' variance
disp("Variance of A2" + var(a2))
disp("Variance of mu2" + var(mu2))
disp("Variance of gamma2" + var(gamma2))
%% Excercise 9: Obtaining Estimates for a Larger d
disp("-----------------------------------------------------------------")
disp("Obtaining Estimates for a Larger d")
size = 1e3; 
d = 5; % Dimension size
steps = 50; % Maximum number of steps
XY = zeros(steps, d);

for i = 2:size
    XY = cat(3, XY, zeros(steps, d)); % Stores data for all steps size times
end

w = zeros(steps + 1, size); % Importance weight function
w(1, :) = 1; % We set the first weight to one
z_2 = zeros(size, 1);
z_1 = ones(size, 1);

for N = 2:steps + 1 
    for i = 1:size
        preSteps = XY(1:N - 1, :, i); % Nxd matrix of previous steps
        currentPos = XY(N - 1, :, i); % Coordinates of last step
        possibleSteps = zeros(2, d);
        
        for dim = 1:d
            nextStep = zeros(1, d);
            nextStep(dim) = 1;
            possibleSteps(1, dim) = 1 - ismember(currentPos + nextStep, preSteps, 'rows'); % value at (1,dim) set to 1 if possible to walk, else 0
            possibleSteps(2, dim) = 1 - ismember(currentPos - nextStep, preSteps, 'rows');% value at (2,dim) set to 1 if possible to walk, else 0
        end
        
        nbrOfWalks = sum(sum(possibleSteps)); % Total number of possible walks

        if(nbrOfWalks == 0) 
            XY(N, :, i) = XY(N - 1, :, i); % If there is no possible walks, set X_k+1 = X_k
            z_2(i) = 0;
            w(N, i) = 0; 
        else
            z_2(i) = 1; % Set indicator function to 1 if it is possible to take one step
            g = 1 / nbrOfWalks; % Probability of our instrumental distribution
            [r, c] = find(possibleSteps, nbrOfWalks); % Returns index for row and cols for values that are nonzero (1)
            rand = randi(nbrOfWalks); % Generate uniform 
            row = r(rand); % row = 1 => positive step, row = 2 => negative step
            col = c(rand); % Which dimension to walk in [1, 2,..., d]
            temp = zeros(1, d);
            temp(col) = 1;
            w(N, i) = z_2(i) / (z_1(i) * g);
            
            if(row == 1)
                XY(N, :, i) = XY(N - 1, :, i) + temp; % Walk positive step in dimension temp
            else 
                XY(N, :, i) = XY(N - 1, :, i) - temp; % Walk negative step in dimension temp
            end
        end
    end
    
    ind = randsample(size, size, true, w(N, :));
    XY = XY(:, :, ind);
    z_1 = z_2(ind);
end

c_nEstMod = zeros(steps, 1); 
c_nEstMod(1) = sum(w(2,:))/size;

for s = 3:steps + 1
    c_nEstMod(s-1) = (sum(w(s, :))/size)*c_nEstMod(s-2);
end

q = 1:1:length(c_nEstMod);
y = @(q) log(c_nEstMod(q));
X1 = ones(length(c_nEstMod), 1);
X2 = @(q) transpose(q * 1);
X3 = @(q) transpose(1 * log(q));
X = [X1(q), X2(q), X3(q)];
b = regress(y(q), X);
a_d = exp(b(1));
mu_d = exp(b(2));
gamma_d = 1+ b(3);
disp(a_d)
disp(mu_d)
disp(gamma_d)
%% Exercise 10: Filter Estimation
% Gathering all values
clear;
load population;

N = 1e3;
n = 50;
tau = zeros(1, n+1); % Filter means array
tauLow = zeros(1, n+1);
tauUp = zeros(1, n+1);
p = @(x, y) unifpdf(y, .5*x, x); % Observation density

for k = 1:n+1
    if k == 1
        part = .2 + rand(N, 1)*.4; % Initialization
    else
        part = (.5 + rand(N, 1)*2.5).*part.*(1-part); % Mutation
    end
    
    w = p(part, Y(k)); % Weighting
    tau(k) = sum(part.*w) / sum(w); % Estimation
    [xx, I] = sort(part);
    cum = cumsum(w(I)) / sum(w);
    ILow = find(cum >= .025, 1);
    IUp = find(cum >= .975, 1);
    tauLow(k) = xx(ILow);
    tauUp(k) = xx(IUp);
    index = randsample(N, N, true, w); % Selection
    part = part(index);
end

counter = 0;

for i = 1:n+1
    if X(i) >= tauLow(i) && X(i) <= tauUp(i)
        counter = counter + 1;
    else
        disp(i) % Display the generations not covered by the CI
    end
end

disp("Within CI: " + counter);

% Plotting tau
green = [27, 158, 119]./255;
orange = [217, 95, 2]./255;
figure()
x = 1:1:51;
p1 = plot(x, tau);
grid on;
xlim([0, 52]);
ylim([0, .9]);
title("Filter Expectation of X")
xlabel("Generation");
ylabel("Relative population size");
p1.Color = green;
p1.Marker = 'o';

% Plotting CI and X_k
figure()
p2 = plot(x, tauLow, '--r');
hold on;
p3 = plot(x, tauUp, '--r');
hold on;
p4 = plot(x, X, '-k');
grid on;
xlim([0, 52]);
ylim([0, .9]);
title("95% Point Wise Confidence Interval and True Values of X")
xlabel("Generation");
ylabel("Relative population size");
hold off;

red = [211 94 96]./255;
black = [128 133 133]./255;
p2.Marker = 'x';
p3.Marker = 'x';
p4.Marker = 'sq';
p2.Color = green;
p3.Color = green;
p4.Color = orange;