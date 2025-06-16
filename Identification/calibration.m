
% Raw data (use dots instead of commas for decimals in MATLAB)
data = [0.9, -0.7, 0, -0.5, -0.7, -0.35, -0.5, -0.2, -1.4, -1.4, ...
        0.4, 0.4, -0.7, -0.9, 0];

%data = [1.9, 1.2, -0.2, 0.1, -0.4, 0, 0.2, 0.4, 0.6, ...
%         -0.2, -0.2, -0.4, -0.9, 0.9, 0];

% Parameters
mu = -0.3766666667;               % Mean
sigma2 = 0.4188809524;            % Variance
%mu = 0.2;               % Mean
%sigma2 = 0.51;            % Variance
sigma = sqrt(sigma2);             % Standard deviation

% Histogram
figure;
histogram(data, 'Normalization', 'pdf','NumBins', 6);  % Normalize to get PDF
hold on;

% X range for normal curve
x = linspace(min(data)-1, max(data)+1, 1000);

% Normal distribution PDF
y = normpdf(x, mu, sigma);

% Plot normal distribution
plot(x, y, 'r', 'LineWidth', 2);
xline([mu-2*sigma,mu+2*sigma])

% Labels and legend
xlabel('x');
ylabel('Probability Density');
title('Histogram with Normal Distribution Overlay');
legend('Data Histogram', ['Normal PDF (\mu = ' num2str(mu) ', \sigma^2 = ' num2str(sigma2) ')']);
xlim([-5,5])
grid on;
