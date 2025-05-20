function plot_mean_ci(time, data, color)
% time: 1 x T vector
% data: N x T matrix (N time series of length T)
% color: RGB triplet or color character (e.g., 'b', [0.2 0.4 0.6])

    % Compute mean and 95% CI
    mu = mean(data, 1);
    SEM = std(data, 0, 1) / sqrt(size(data,1));
    CI95 = 1.96 * SEM;  % for 95% CI

    % Plot shaded area
    hold on;
    fill([time, fliplr(time)], ...
         [mu + CI95, fliplr(mu - CI95)], ...
         color, 'FaceAlpha', 0.5, 'EdgeColor', 'none');

    % Plot mean line
    plot(time, mu, 'Color', color, 'LineWidth', 2);
end

% Example
% Simulated data
% T = 100;
% time = linspace(0, 10, T);
% data = sin(time) + 0.2*randn(30, T);  % 30 trials with noise
% 
% % Plot
% figure;
% plot_mean_ci(time, data, [0 0.447 0.741]);  % blue shade/line
% xlabel('Time'); ylabel('Signal');
% title('Mean with 95% Confidence Interval');
