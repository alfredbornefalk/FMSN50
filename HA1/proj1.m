%% Part 1 of 2a: Monte Carlo Sampling
clear; clc; % Good practice

disp("Standard Monte Carlo")
load powercurve_D236;

month = ["January", "February", "March", "April", "May", "June", "July", "August", "September", "October", "November", "December"];
lambda = [11.7, 10.7, 10.1, 8.8, 8.6, 8.9, 8.6, 8.9, 10.0, 10.9, 11.7, 11.7];
k = [2.0, 2.0, 2.0, 1.9, 1.9, 1.9, 1.9, 1.9, 2.0, 1.9, 2.0, 2.0];
lower = .025;
upper = 1 - lower;
size = 1e6;
sizeSq = sqrt(size);
pVal = 1.96;

for c = 1:12
    v = wblrnd(lambda(c), k(c), size, 1);
    power = P(v);
    tau = mean(power); % Law of Large Numbers
    sd = std(power);
    LB = tau - (pVal * sd / sizeSq);
    UB = tau + (pVal * sd / sizeSq);
    width = UB - LB;
    disp(month(c) + ": LB = " + LB + ", UB = " + UB + "; Width = " + width)
end
%% Part 2 of 2a: Weibull Distribution
disp("-----------------------------------------------------------------")
disp("Weibull Distribution")
a = 3;
b = 30;

for c = 1:12
    fxA = wblcdf(a, lambda(c), k(c));
    fxB = wblcdf(b, lambda(c), k(c));
    u = rand(size, 1);
    par = fxA + u * (fxB - fxA);
    X = wblinv(par, lambda(c), k(c));
    power = P(X) * (fxB - fxA);
    tau = mean(power); % Law of Large Numbers
    sd = std(power);
    LB = tau - (pVal * sd / sizeSq);
    UB = tau + (pVal * sd / sizeSq);
    width = UB - LB;
    disp(month(c) + ": LB = " + LB + ", UB = " + UB + "; Width = " + width)
end
%% 2b: Control Variate
disp("-----------------------------------------------------------------")
disp("Control Variate")

for c = 1:12
   v = wblrnd(lambda(c), k(c), size, 1); % Control variate
   power = P(v); % phi(v)
   temp = -cov([power, v]) / var(v);
   beta = temp(1, 2);
   zeta = power; % zeta = r since E[v] = m
   tauCV = mean(zeta);
   tempMod = cov([power, v]);
   sdCV = sqrt(var(power) + (beta^2) * var(v) + 2*beta*tempMod(1, 2));
   LB = tauCV - pVal * sdCV / sizeSq;
   UB = tauCV + pVal * sdCV / sizeSq;
   width = UB - LB;
   % Calculating the correlation
   tempRho = cov([power, v]) / (sqrt(var(power))*sqrt(var(v)));
   rho = tempRho(1, 2);
   disp(month(c) + ": LB = " + LB + ", UB = " + UB + "; Width = " + width + "; Corr = " + rho)
end
%% 2c: Importance Sampling 
disp("-----------------------------------------------------------------")
disp("Importance Sampling")

for c = 1:12
    x = linspace(0, 35, size);
    sigma = 5.08; % standard deviation of g(x)
    f = @(x) wblpdf(x, lambda(c), k(c));
    phi = @(x) P(x);
    phiTimesF = phi(x) .* f(x)';
    xValue = find(phiTimesF == max(phiTimesF), 1, 'first');
    index = x(xValue); % x-value corresponding to the maximum value of y
    g = @(x) normpdf(x, index, sigma);
    X = sigma.*randn(1,size) + index; % generating samples from g(x)
    target_func = @(x) (f(x) .* P(x)') ./ g(x); % phi(x) * f(x)/g(x) (slide 16)
    Y = target_func(X);
    tau = mean(Y);
    sd = std(Y);
    LB = tau - pVal * sd / sizeSq;
    UB = tau + pVal * sd / sizeSq;
    width = UB - LB;
    disp(month(c) + ": LB = " + LB + ", UB = " + UB + "; Width = " + width)
end

figure()
plot(x, phiTimesF) % Figure 1 in report
xlabel('X')
ylabel('Y')
%% 2d: Antithesis Sampling
disp("-----------------------------------------------------------------")
disp("Antithesis Sampling")

for c = 1:12
    u = rand(size / 2, 1); % Half the sample size
    power = P(wblinv(u, lambda(c), k(c)));
    powerTilde = P(wblinv(1 - u, lambda(c), k(c)));
    W = (power + powerTilde) ./ 2;
    tau = mean(W);
    LB = tau - pVal * std(W) / sqrt(size / 2);
    UB = tau + pVal * std(W) / sqrt(size / 2);
    width = UB - LB;
    disp(month(c) + ": LB = " + LB + ", UB = " + UB + "; Width = " + width)
end
%% 2e: Probability of Power Generation
disp("-----------------------------------------------------------------")
disp("Probability of Power Generation")

for c = 1:12
    v = wblrnd(lambda(c), k(c), size, 1);
    power = P(v);
    probPower = nnz(power) / numel(power); % Probability that the turbine generates power 
    disp(month(c) + ": Probability of power = " + probPower)
end
%% 2f: Power Coefficient
disp("-----------------------------------------------------------------")
disp("Power Coefficient")

d = 236;
rho = 1.225;
m = 3;

for c = 1:12
    v = wblrnd(lambda(c), k(c), size, 1);
    power = P(v);
    powerExp = mean(power);
    powerTotExp = .5 * rho * pi * .25 * d^2 * gamma(1 + m / k(c))*lambda(c)^m;
    powerCoeff = powerExp / powerTotExp;
    sdPowerCoeff = std(power ./ powerTotExp);
    LB = powerCoeff - (pVal * sdPowerCoeff / sizeSq);
    UB = powerCoeff + (pVal * sdPowerCoeff / sizeSq);
    width = UB - LB;
    disp(month(c) + ": LB = " + LB + ", UB = " + UB + "; Width = " + width)
end
%% 2g: Capacity & Availability Factor
disp("-----------------------------------------------------------------")
disp("Capacity & Availability Factor")

accCapCoeff = 0;
accProb = 0;

for c = 1:12
    v = wblrnd(lambda(c), k(c), size, 1);
    power = P(v);
    powerExp = mean(power);
    capCoeff = powerExp / 15e6; % 15e6 = maximum possible output per second
    accCapCoeff = accCapCoeff + capCoeff;
    probPower = nnz(power) / numel(power); % Probability that the turbine generates power
    accProb = accProb + probPower;
end

meanCapCoeff = accCapCoeff / 12;
disp("Average Capacity Factor: " + meanCapCoeff)
meanProb = accProb / 12;
disp("Average Availability Factor: " + meanProb)
%% 3a: Expected Combined Power Generation
disp("-----------------------------------------------------------------")
disp("Expected Combined Power Generation")

% Given parameters
alpha = .638;
p = 3;
q = 1.5;
kVal = 1.95;
lambdaVal = 10.05;

% PDF, CDF, Joint PDF & Joint CDF
pdf = @(x) wblpdf(x, lambdaVal, kVal);
cdf = @(x) wblcdf(x, lambdaVal, kVal);
pdfJoint = @(x, y) pdf(x).*pdf(y) .*(1+alpha*(1-cdf(x).^p).^(q-1).*(1-cdf(y).^p).^(q-1).*((cdf(x).^p).*(1+p*q)-1).*((cdf(y).^p).*(1+p*q)-1));
cdfJoint = @(x, y) cdf(x).*cdf(y) .*(1+alpha*(1-cdf(x).^p).^q.*(1-cdf(y).^p).^q);

x = linspace(0, 35, size);
sigma = 5.08; % standard deviation of g(x), same value as in 2c
phi = @(x) P(x);
phiTimesPDF = phi(x) .* pdf(x)';
xValue = find(phiTimesPDF == max(phiTimesPDF), 1, 'first');
index = x(xValue); % x-value corresponding to the maximum value of the function
g = @(x) normpdf(x, index, sigma);
X1 = sigma.*randn(1,size) + index; % generating samples from g(x1)
X2 = sigma.*randn(1,size) + index; % generating samples from g(x2)
target_func = @(x) (pdf(x) .* P(x)') ./ g(x); % phi(x) * f(x)/g(x) (slide 16)
Y1 = target_func(X1);
Y2 = target_func(X2);
tau1 = mean(Y1); % E(Pv1)
tau2 = mean(Y2); % E(Pv2)
tauSum = tau1 + tau2; % E(Pv1) + E(Pv2)
disp("Expected Value: " + tauSum)
%% 3b: Covariance of Combined Power Generation
disp("-----------------------------------------------------------------")
disp("Covariance of Combined Power Generation")

mu = [10.960 10.960];
sigma = [20 1; 1 20];
R = mvnrnd(mu, sigma, size);
X = R(:,1);
Y = R(:, 2);
omega = @(x, y) (P(x).*P(y)) .* pdfJoint(x, y) ./ mvnpdf([x y], mu, sigma);

res = mean(omega(X, Y));
cov = res - tau1 * tau2;
disp("Covariance: " + cov)
%% 3c: Variance & SD of Combined Power Generation
disp("-----------------------------------------------------------------")
disp("Variance & SD of Combined Power Generation")
V = wblrnd(lambdaVal, kVal, size, 1);
power = P(V);
varPhi = var(power);
varTot = varPhi * 2 + 2*cov;
disp("Variance: " + varTot)

sdTot = sqrt(varTot);
disp("Standard Deviation: " + sdTot)
%% 3d: 95 % Confidence Interval for Combined Power Generation Outputs
disp("-----------------------------------------------------------------")
disp("95 % Confidence Interval for Combined Power Generation Outputs")
muOver = [9.00, 9.00];
sigmaOver = [30, 1;1, 30];
muUnder = [6.00, 6.00];
sigmaUnder = [40, 1; 1, 40];
ROver = mvnrnd(muOver, sigmaOver, size);
XOver = ROver(:, 1);
YOver = ROver(:, 2);
RUnder = mvnrnd(muUnder, sigmaUnder, size);
XUnder = RUnder(:, 1);
YUnder = RUnder(:, 2);
omegaUnder = @(x, y) ((P(x) + P(y)) < 15e6) .* pdfJoint(x, y) ./ mvnpdf([x y], muUnder, sigmaUnder);
omegaOver = @(x, y) ((P(x) + P(y)) > 15e6) .* pdfJoint(x, y) ./ mvnpdf([x y], muOver, sigmaOver);
underMean = mean(omegaUnder(XUnder, YUnder));
overMean = mean(omegaOver(XOver, YOver));
underSD = std(omegaUnder(XUnder, YUnder));
overSD = std(omegaOver(XOver, YOver));
sumUnderOver = underMean + overMean;
underLB = underMean - pVal * underSD / sqrt(size);
underUB = underMean + pVal * underSD / sqrt(size);
overLB = overMean - pVal * overSD / sqrt(size);
overUB = overMean + pVal * overSD / sqrt(size);
disp("Over: LB = " + overLB + ", UB = " + overUB)
disp("Under: LB = " + underLB + ", UB = " + underUB)
disp("Sum over and under 15 MW: " + sumUnderOver)
%% Figure 2
omega = @(x, y) (P(x) .* P(y)) .* pdfJoint(x, y);
V1 = linspace(0,35,100);
V2 = linspace(0,35,100);
Z = zeros(100,100);

for r=1:100
    for c=1:100
       Z(r,c) =  omega(V1(r),V2(c));
    end
end

figure()
surf(V1,V2,Z);
colorbar;
xlabel('X')
ylabel('Y')
%% Figure 3
mu = [10.960 10.960];
sigma = [20 1; 1 20];
omega = @(x, y) (P(x) .* P(y)) .* pdfJoint(x, y) / mvnpdf([x y], mu, sigma);
V1 = linspace(-10,50,100);
V2 = linspace(-10,50,100);
Z = zeros(100,100);

for r=1:100
    for c=1:100
       Z(r,c) =  omega(V1(r),V2(c));
    end
end

figure()
surf(V1,V2,Z);
colorbar;
xlabel('X')
ylabel('Y')
%% Figure 4
mu = [9 9];
sigma = [30 1; 1 30];
omega = @(x, y) ((P(x) + P(y)) > 15e6) .* pdfJoint(x, y);
V1 = linspace(0,35,100);
V2 = linspace(0,35,100);
Z = zeros(100,100);

for r=1:100
    for c=1:100
       Z(r,c) =  omega(V1(r),V2(c));
    end
end

figure() 
surf(V1,V2,Z);
colorbar;
xlabel('X')
ylabel('Y')
%% Figure 5
mu = [6 6];
sigma = [40 1; 1 40];
omega = @(x, y) ((P(x) + P(y)) < 15e6) .* pdfJoint(x, y);
V1 = linspace(0,35,100);
V2 = linspace(0,35,100);
Z = zeros(100,100);

for r=1:100
    for c=1:100
       Z(r,c) =  omega(V1(r),V2(c));
    end
end

figure()
surf(V1,V2,Z);
colorbar;
xlabel('X')
ylabel('Y')
%% Figure 6
mu = [9 9];
sigma = [30 1; 1 30];
omega = @(x, y) ((P(x) + P(y)) > 15e6) .* pdfJoint(x, y) ./ mvnpdf([x y], mu, sigma);

V1 = linspace(-10, 50,100);
V2 = linspace(-10, 50,100);
Z = zeros(100,100);

for r=1:100
    for c=1:100
       Z(r,c) =  omega(V1(r),V2(c));
    end
end

figure() 
surf(V1,V2,Z);
colorbar;
xlabel('X')
ylabel('Y')
%% Figure 7
mu = [6 6];
sigma = [60 1; 1 60];
omega = @(x, y) ((P(x) + P(y)) < 15e6) .* pdfJoint(x, y) ./ mvnpdf([x y], mu, sigma);

V1 = linspace(-10, 50,100);
V2 = linspace(-10, 50,100);
Z = zeros(100,100);

for r=1:100
    for c=1:100
       Z(r,c) =  omega(V1(r),V2(c));
    end
end

figure() 
surf(V1,V2,Z);
colorbar;
xlabel('X')
ylabel('Y')
%% Figure 8
omega = @(x, y) ((P(x) + P(y)));
V1 = linspace(0, 35,100);
V2 = linspace(0, 35,100);
Z = zeros(100,100);

for r=1:100
    for c=1:100
       Z(r,c) =  omega(V1(r),V2(c));
    end
end

figure()
surf(V1,V2,Z);
colorbar;
xlabel('X')
ylabel('Y')