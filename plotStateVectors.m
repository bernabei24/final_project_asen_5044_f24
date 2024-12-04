function plotStateVectors(x_lin, dt)
    % This function plots all 6 state vectors over time.
    % Inputs:
    % - x_lin: A 6-by-N matrix where each row represents a state variable and
    %          each column represents a timestep.
    % - dt: Time step interval (in seconds).

    % Define time vector
    timesteps = size(x_lin, 2);
    time = 0.1:dt:(timesteps) * dt;

    % Plot each state vector
    figure;
    for i = 1:6
        subplot(6, 1, i); % Create a 6x1 grid of subplots
        plot(time, x_lin(i, :), 'LineWidth', 1.5); % Plot the state vector
        grid on;
        title(['State Vector ', num2str(i)]);
        xlabel('Time (s)');
        ylabel(['x_', num2str(i)]);
    end

    % Adjust layout
    sgtitle('State Vectors Over Time'); % Add a super-title to the figure
end
