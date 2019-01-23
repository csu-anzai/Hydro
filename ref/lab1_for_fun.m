clc, clear all, clf

% Part 1:

% define size of vector
n = 100000

% create n random numbers
RNG = rand(1,n);

% Part 2:

% produce a histogram of random data
figure(1)
histogram( RNG )

% compute the mean of the data
meanRNG = mean( RNG )

% compute the standard deviation
stdRNG = std( RNG )

% Part 3:

% mean & variance of gaussian
mu = 0;
sigma = 1.0;

% compute PDF of normal distribution
u = sort(RNG);
PDF = norminv( u, mu, sigma );

% Part 4:

% plot gaussian PDF
figure(2)
histogram( PDF )

% compute mean and variance
meanPDF = mean(PDF)
varPDF = std(PDF)^2

% Part 5:

% create Q-Q plot
figure(3)
qqplot( PDF );

% Part 6:

% perform the kolmogorov-smirnov test
[H,P] = kstest( PDF )
