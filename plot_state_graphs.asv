function plot_state_graphs(delx_plus, P_plus, tvec)
    % Validate input dimensions
    if size(delx_plus, 2) ~= length(tvec)
        error('The number of columns in delx_plus must match the length of tvec.');
    end
    
    % Prepare figure
    figure;
    
    % Number of states
    num_states = size(delx_plus, 1);
    
    % Plot each state with its error bounds
    for i = 1:num_states
        subplot(num_states, 1, i);
        
        % Extract the state estimate and its standard deviation
        state_est = delx_plus(i, :); % Directly use as a row vector
        state_std = sqrt(squeeze(P_plus(i, i, :)))'; % Standard deviations from diagonal of P_plus
        
        % Time to plot
        time = tvec;
        
        % Plot the state estimate
        plot(time, state_est, 'LineWidth', 1.5);
        hold on;
        
        % Plot the error bounds at +-2 standard deviations
        fill([time, fliplr(time)], [state_est + 2 * state_std, fliplr(state_est - 2 * state_std)], ...
            'r', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
        
        % Labels and titles
        title(sprintf('State %d', i));
        ylabel(sprintf('State %d', i));
        xlabel('Time (s)');
        
        % Grid and formatting
        grid on;
        hold off;
    end
    
    % Super title for the figure
    sgtitle('State Estimates with \pm2 Standard Deviations Error Bounds');
end
