clear; clc; close all;
%% 1c: The Behavior of the Chain for Different Breakpoints
disp("The Behavior of the Chain for Different Breakpoints")

% # of disasters
load("coal_mine_disasters.mat")

% Defining start & end points 
tStart = 1658;
tEnd = 1980;

% Defining the initial breakpoint selections
tFirst = 1690;
tLast = 1960;

d = 6;
samples = 1e4;
burnIn = 1e3;
psi = 20;
tau = T;

for count = 2:d
    step = (tLast - tFirst) / count;
    tMid = tFirst:step:tLast;
    t = [tStart, tMid(2:end - 1), tEnd];
    
    rho = .01 * ones(count, 1);
    
    total = zeros(samples, length(t));
    theta = gamrnd(2, 1 / psi);
    lambda = gamrnd(2, 1 / theta, 1, count);
    
    for i = 1:burnIn
        theta = gamrnd(2 * length(lambda) + 2, 1 ./ (psi + sum(lambda)));
        
        samplesTemp = zeros(1, length(t) - 1);
        
        for j = 1:length(t) - 1
            samplesTemp(j) = sum((t(j) <= tau) & (tau < t(j + 1)));
        end
        
        lambda = gamrnd(samplesTemp' + 2, 1./(theta + (t(2:end) - t(1:end - 1))'));
        
        accepted = zeros(1, length(t) - 2);
        
        for k = 2:length(t) - 1
            R = rho(1) * (t(k + 1) - t(k - 1));
            Xstar =  t(k) - R + 2 * R * rand;
            
            while (Xstar < t(k - 1) || Xstar >  t(k + 1))
                Xstar =  t(k) - R + 2 * R * rand;
            end
            
            num = formula(lambda, [t(1:k - 1), Xstar, t(k + 1:end)], tau);
            den = formula(lambda, t, tau);
            alpha = min(1, num / den);
            U = rand(1);

            if (U <= alpha)
                t(k) = Xstar;
                accepted(k - 1) = accepted(k - 1) + 1;
            end
        end
    end
    
    for i = 1:samples
        theta = gamrnd(2 * length(lambda) + 2, 1 ./ (psi + sum(lambda)));
        
        samplesTemp = zeros(1, length(t) - 1);
        
        for j = 1:length(t)-1
            samplesTemp(j) = sum((t(j) <= tau) & (tau < t(j + 1)));
        end
        
        lambda = gamrnd(samplesTemp' + 2, 1./(theta + (t(2:end) - t(1:end - 1))'));
        
        accepted = zeros(1, length(t) - 2);
        
        for k = 2:length(t) - 1
            R = rho(1) * (t(k + 1) - t(k - 1));
            Xstar =  t(k) - R + 2 * R * rand;
            
            while (Xstar < t(k-1) || Xstar >  t(k+1))
                Xstar =  t(k) - R + 2 * R * rand;
            end

            num = formula(lambda, [t(1:k - 1), Xstar, t(k + 1:end)], tau);
            den = formula(lambda, t, tau);
            alpha = min(1, num / den);
            U = rand(1);

            if (U <= alpha)
                t(k) = Xstar;
                accepted(k - 1) = accepted(k - 1) + 1;
            end
        end
        
        total(i, :) = t;
    end
    
    figure();
    plot(total);
    
    if (count == 2)
        title("The Behavior of the Chain For " + (count - 1) + " Breakpoint");
    else
        title("The Behavior of the Chain For " + (count - 1) + " Breakpoints");
    end
    
    xlabel("Samples");
    ylabel("Year");
    
    % Plotting 5 breakpoint markings in the coal-mine-disasters-plot
    if (count == 6)
        figure();
        plot(T, 1:length(T));
        hold on;

        for i = 2:d
            line([mean(total(:, i)) mean(total(:, i))], [0, length(T)], "Color", [rand, rand, rand]);
        end

        hold off;
        title("Total Number of Accidents From 1658-1980");
        xlabel("Year");
        ylabel("Accidents");
        axis([tStart, T(end), 0, length(T)]);
    end
end
%% 1d: How Psi Affects Theta, Lambda & t
disp("-----------------------------------------------------------------")
disp("How Psi Affects Theta, Lambda & t")

step = (tLast - tFirst) / d;
tMid = tFirst:step:tLast;
t =[tStart, tMid(2:end - 1), tEnd];

psi = 50;
rho = .01 * ones(d, 1);

thetaMean = zeros(psi, 1);
thetaVar = zeros(psi, 1);
lambdaMean = zeros(psi, d);
lambdaVar = zeros(psi, d);
tArr = zeros(samples, length(t));
tMean = zeros(psi, length(t));
tVar = zeros(psi, length(t));

for p = 1:psi
    theta = gamrnd(2, 1 / p);
    lambda = gamrnd(2, 1 / theta, 1, d);
    
    thetaTemp = zeros(samples, 1);
    lambdaTemp = zeros(samples, d);
    
    for i = 1:burnIn
        theta = gamrnd(2 * length(lambda) + 2, 1 ./ (p + sum(lambda)));
        
        samplesTemp = zeros(1, length(t) - 1);
        
        for j = 1:length(t) - 1
            samplesTemp(j) = sum((t(j) <= tau) & (tau < t(j + 1)));
        end
        
        lambda = gamrnd(samplesTemp' + 2, 1./(theta + (t(2:end) - t(1:end - 1))'));
        
        accepted = zeros(1, length(t) - 2);
        
        for k = 2:length(t) - 1
            R = rho(1) * (t(k + 1) - t(k - 1));
            Xstar =  t(k) - R + 2 * R * rand;
            
            while (Xstar < t(k - 1) || Xstar >  t(k + 1))
                Xstar =  t(k) - R + 2 * R * rand;
            end

            num = formula(lambda, [t(1:k - 1), Xstar, t(k + 1:end)], tau);
            den = formula(lambda, t, tau);
            alpha = min(1, num / den);
            U = rand(1);

            if (U <= alpha)
                t(k) = Xstar;
                accepted(k - 1) = accepted(k - 1) + 1;
            end
        end
    end
    
    for i = 1:samples
        theta = gamrnd(2 * length(lambda) + 2, 1 ./ (p + sum(lambda)));
        
        samplesTemp = zeros(1, length(t) - 1);
        
        for j = 1:length(t) - 1
            samplesTemp(j) = sum((t(j) <= tau) & (tau < t(j + 1)));
        end
        
        lambda = gamrnd(samplesTemp' + 2, 1./(theta + (t(2:end) - t(1:end - 1))'));
        
        accepted = zeros(1, length(t) - 2);
        
        for k = 2:length(t) - 1
            R = rho(1) * (t(k + 1) - t(k - 1));
            Xstar =  t(k) - R + 2 * R * rand;
            
            while (Xstar < t(k - 1) || Xstar >  t(k + 1))
                Xstar =  t(k) - R + 2 * R * rand;
            end

            num = formula(lambda, [t(1:k - 1), Xstar, t(k + 1:end)], tau);
            den = formula(lambda, t, tau);
            alpha = min(1, num / den);
            U = rand(1);

            if (U <= alpha)
                t(k) = Xstar;
                accepted(k - 1) = accepted(k - 1) + 1;
            end
        end
        
        thetaTemp(i) = theta;
        lambdaTemp(i, :) = lambda';
        tArr(i, :) = t;
    end
    
    thetaMean(p) = mean(thetaTemp);
    thetaVar(p) = var(thetaTemp);
    lambdaMean(p, :) = mean(lambdaTemp);
    lambdaVar(p, :) = var(lambdaTemp);
    tMean(p, :) = mean(tArr);
    tVar(p, :) = var(tArr);
end

figure();
plot(thetaMean, "-x");
title("The Mean of the Parameter \theta's Dependency on \psi");
xlabel("\psi");
ylabel("Mean");
grid on;

figure();
plot(thetaVar, "-x");
title("The Variance of the Parameter \theta's Dependency on \psi");
xlabel("\psi");
ylabel("Variance");
grid on;

figure();
plot(lambdaMean, "-x");
title("The Mean of the Intensities \lambda's Dependency on \psi");
xlabel("\psi");
ylabel("Mean");
legend("\lambda_1", "\lambda_2", "\lambda_3", "\lambda_4", "\lambda_5", "\lambda_6");

figure();
plot(lambdaVar, "-x");
title("The Variance of the Intensities \lambda's Dependency on \psi");
xlabel("\psi");
ylabel("Variance");
legend("\lambda_1", "\lambda_2", "\lambda_3", "\lambda_4", "\lambda_5", "\lambda_6");

figure();
plot(tMean(:, 2:end - 1), "-x");
title("The Mean of the Breakpoints t's Dependency on \psi");
xlabel("\psi");
ylabel("Mean");
legend("t_1", "t_2", "t_3", "t_4", "t_5");

figure();
plot(tVar(:, 2:end - 1), "-x");
title("The Variance of the Breakpoints t's Dependency of \psi");
xlabel("\psi");
ylabel("Variance");
legend("t_1", "t_2", "t_3", "t_4", "t_5");
%% 1e: How Rho Affects Theta, Lambda & t
disp("-----------------------------------------------------------------")
disp("How Rho Affects Theta, Lambda & t")
% How Rho Affects Theta & Lambda

nor = 30;
psi = 20;

rhoTemp = (1:nor) * .01;
rho = zeros(d, nor);

for i = 1:d
   rho(i, :) = rhoTemp; 
end

thetaMean = zeros(nor, 1);
lambdaMean = zeros(nor, d);

for n = 1:nor
    theta = gamrnd(2, 1 / psi);
    lambda = gamrnd(2, 1 / theta, 1, d);
    thetaTemp = zeros(samples, 1);
    lambdaTemp = zeros(samples, d);
    
    for i = 1:burnIn
        theta = gamrnd(2 * length(lambda) + 2, 1 ./ (psi + sum(lambda)));
        
        samplesTemp = zeros(1, length(t) - 1);
        
        for j = 1:length(t) - 1
            samplesTemp(j) = sum((t(j) <= tau) & (tau < t(j + 1)));
        end
        
        lambda = gamrnd(samplesTemp' + 2, 1./(theta + (t(2:end) - t(1:end-1))'));
        
        accepted = zeros(1, length(t) - 2);
        
        for k = 2:length(t) - 1
            R = rho(1, n) * (t(k + 1) - t(k - 1));
            Xstar =  t(k) - R + 2 * R * rand;
            
            while (Xstar < t(k - 1) || Xstar >  t(k + 1))
                Xstar =  t(k) - R + 2 * R * rand;
            end

            num = formula(lambda, [t(1:k - 1), Xstar, t(k + 1:end)], tau);
            den = formula(lambda, t, tau);
            alpha = min(1, num / den);
            U = rand(1);

            if (U <= alpha)
                t(k) = Xstar;
                accepted(k - 1) = accepted(k - 1) + 1;
            end
        end
    end
    
    for i = 1:samples
        theta = gamrnd(2 * length(lambda) + 2, 1 ./ (psi + sum(lambda)));
        
        samplesTemp = zeros(1, length(t) - 1);
        
        for j = 1:length(t) - 1
            samplesTemp(j) = sum((t(j) <= tau) & (tau < t(j + 1)));
        end
        
        lambda = gamrnd(samplesTemp' + 2, 1./(theta + (t(2:end) - t(1:end-1))'));
        
        accepted = zeros(1, length(t) - 2);
        
        for k = 2:length(t) - 1
            R = rho(1, n) * (t(k + 1) - t(k - 1));
            Xstar =  t(k) - R + 2 * R * rand;
            
            while (Xstar < t(k - 1) || Xstar >  t(k + 1))
                Xstar =  t(k) - R + 2 * R * rand;
            end

            num = formula(lambda, [t(1:k - 1), Xstar, t(k + 1:end)], tau);
            den = formula(lambda, t, tau);
            alpha = min(1, num / den);
            U = rand(1);

            if (U <= alpha)
                t(k) = Xstar;
                accepted(k - 1) = accepted(k - 1) + 1;
            end
        end
        
        thetaTemp(i) = theta;
        lambdaTemp(i, :) = lambda';
    end
    
    thetaMean(n) = mean(thetaTemp);
    lambdaMean(n, :) = mean(lambdaTemp);
end

figure();
plot(rho(1, :), thetaMean, "-x");
title("The Mean of the Parameter \theta's Dependency on \rho");
xlabel("\rho");
ylabel("Mean");

figure();
plot(rho(1, :), lambdaMean, "-x");
title("The Mean of the Intensities \lambda's Dependency on \rho");
xlabel("\rho");
ylabel("Mean");
legend("\lambda_1", "\lambda_2", "\lambda_3", "\lambda_4", "\lambda_5", "\lambda_6");

% How Rho Affects t
rArr = [.01, .02, .03, .04];
nor = length(rArr);
rho = zeros(d, nor);

for i = 1:d
   rho(i, :) = rArr; 
end

tArr = zeros(samples, length(t));

for n = 1:nor
    tArr = zeros(samples, length(t));
    theta = gamrnd(2, 1 / psi);
    lambda = gamrnd(2, 1 / theta, 1, d);
    
    for i = 1:burnIn
        theta = gamrnd(2 * length(lambda) + 2, 1 ./ (psi + sum(lambda)));
        
        samplesTemp = zeros(1, length(t) - 1);
        
        for j = 1:length(t) - 1
            samplesTemp(j) = sum((t(j) <= tau) & (tau < t(j + 1)));
        end
        
        lambda = gamrnd(samplesTemp' + 2, 1./(theta + (t(2:end) - t(1:end-1))'));
        
        accepted = zeros(1, length(t) - 2);
        
        for k = 2:length(t) - 1
            R = rho(1, n) * (t(k + 1) - t(k - 1));
            Xstar =  t(k) - R + 2 * R * rand;
            
            while (Xstar < t(k - 1) || Xstar >  t(k + 1))
                Xstar =  t(k) - R + 2 * R * rand;
            end

            num = formula(lambda, [t(1:k - 1), Xstar, t(k + 1:end)], tau);
            den = formula(lambda, t, tau);
            alpha = min(1, num / den);
            U = rand(1);

            if (U <= alpha)
                t(k) = Xstar;
                accepted(k - 1) = accepted(k - 1) + 1;
            end
        end
    end
    
    for i = 1:samples
        theta = gamrnd(2 * length(lambda) + 2, 1 ./ (psi + sum(lambda)));
        
        samplesTemp = zeros(1, length(t) - 1);
        
        for j = 1:length(t) - 1
            samplesTemp(j) = sum((t(j) <= tau) & (tau < t(j + 1)));
        end
        
        lambda = gamrnd(samplesTemp' + 2, 1./(theta + (t(2:end) - t(1:end-1))'));
        
        accepted = zeros(1, length(t) - 2);
        
        for k = 2:length(t) - 1
            R = rho(1, n) * (t(k + 1) - t(k - 1));
            Xstar =  t(k) - R + 2 * R * rand;
            
            while (Xstar < t(k - 1) || Xstar >  t(k + 1))
                Xstar =  t(k) - R + 2 * R * rand;
            end

            num = formula(lambda, [t(1:k - 1), Xstar, t(k + 1:end)], tau);
            den = formula(lambda, t, tau);
            alpha = min(1, num / den);
            U = rand(1);

            if (U <= alpha)
                t(k) = Xstar;
                accepted(k - 1) = accepted(k - 1) + 1;
            end
        end
        
        tArr(i, :) = t;
    end
    
    figure();
    subplot(2, 2, 1);
    autocorr(tArr(:, 2), 5e2)
    title("The Correlation Function for t_1");
    xlabel("Time lag");
    ylabel("Dependency on \rho = " + rho(1, n));
    
    subplot(2, 2, 2)
    autocorr(tArr(:, 3), 5e2)
    title("The Correlation Function for t_2");
    xlabel("Time lag");
    ylabel("Dependency on \rho = " + rho(1, n));
    
    subplot(2, 2, 3)
    autocorr(tArr(:, 4), 5e2)
    title("The Correlation Function for t_3");
    xlabel("Time lag");
    ylabel("Dependency on \rho = " + rho(1, n));
    
    subplot(2, 2, 4)
    autocorr(tArr(:, 5), 5e2)
    title("The Correlation Function for t_4");
    xlabel("Time lag");
    ylabel("Dependency on \rho = " + rho(1, n));
    
    if (n == 3)
        figure();
        plot(tArr);
        title("The Behavior of the Chain For " + (d - 1) + " Breakpoints with \rho = " + rho(1, n));
        xlabel("Samples");
        ylabel("Year");
    end
end

% Acceptance Ratio Plot
rhoPlot = linspace(.001, .05);
acceptedPlot = zeros(length(rhoPlot), d - 1);
accSamp = 100;

for r = 1:length(rhoPlot)
   for i = 1:accSamp
        theta = gamrnd(2 * length(lambda) + 2, 1 ./ (psi + sum(lambda)));
        
        samplesTemp = zeros(1, length(t) - 1);
        
        for j = 1:length(t) - 1
            samplesTemp(j) = sum((t(j) <= tau) & (tau < t(j + 1)));
        end
        
        lambda = gamrnd(samplesTemp' + 2, 1./(theta + (t(2:end) - t(1:end-1))'));
        
        accepted = zeros(1, length(t) - 2);
        
        for k = 2:length(t) - 1
            R = rhoPlot(r) * (t(k + 1) - t(k - 1));
            Xstar =  t(k) - R + 2 * R * rand;
            
            while (Xstar < t(k - 1) || Xstar >  t(k + 1))
                Xstar =  t(k) - R + 2 * R * rand;
            end

            num = formula(lambda, [t(1:k - 1), Xstar, t(k + 1:end)], tau);
            den = formula(lambda, t, tau);
            alpha = min(1, num / den);
            U = rand(1);

            if (U <= alpha)
                t(k) = Xstar;
                accepted(k - 1) = accepted(k - 1) + 1;
            end
        end
       
       acceptedPlot(r, :) = acceptedPlot(r, :) + accepted;
   end
end

ratio = sum(acceptedPlot, 2) / (accSamp * (d - 1));

figure();
plot(rhoPlot, ratio, "-x");
hold on;
plot([.001, .05],[.3, .3]);
title("Acceptance Ratio as a Function of \rho");
xlabel("\rho");
ylabel("Acceptance ratio");
hold off;

%% 2b: Bootstrapped 95% Confidence Intervals
disp("-----------------------------------------------------------------")
disp("Bootstrapped 95% Confidence Intervals")

load("atlantic.txt")

fInv = @(u, beta, mu) mu - beta * log(-log(u));
n = length(atlantic);
B = 1e3;
[betaHat, muHat] = est_gumbel(atlantic);
boot = zeros(2, B);

for b = 1:B % Bootstrap
    yBoot =  fInv(rand(n, 1), betaHat, muHat);
    [betaTemp, muTemp] = est_gumbel(yBoot);
    boot(1, b) = betaTemp;
    boot(2, b) = muTemp;
end

deltaBeta = sort(boot(1, :) - betaHat);
deltaMu = sort(boot(2, :) - muHat);
alpha = .05;
lowBeta = betaHat - deltaBeta(ceil((1 - alpha / 2) * B));
upperBeta = betaHat - deltaBeta(ceil(alpha * B / 2));
lowMu = muHat - deltaMu(ceil((1 - alpha / 2) * B));
upperMu = muHat - deltaMu(ceil(alpha * B / 2));
disp("Beta confidence interval LB: " + lowBeta + " UB: " + upperBeta)
disp("Mu confidence interval LB: " + lowMu + " UB: " + upperMu)
%% 2c: One-Sided 95% Confidence Interval
disp("-----------------------------------------------------------------")
disp("One-Sided 95% Confidence Interval")

T = 3 * 14 * 100;
waves = zeros(1, B);

for i = 1:B
    waves(1, i) = fInv(1 - 1 / T, boot(1, i), boot(2, i));
end

averageWave = fInv(1 - 1 / T, betaHat, muHat);
deltaWave = sort(waves(1, :) - averageWave);
disp(averageWave - deltaWave(ceil(alpha * B))) % Upper bound
