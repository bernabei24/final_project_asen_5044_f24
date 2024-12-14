function plot_states_and_errors(delx_plus, x_nominal, x_true, P_plus, tvec)
    % Calculate the total state by adding the deviations from the nominal trajectory to the nominal state
    total_state = x_nominal + delx_plus';  % Transpose x_nominal to match dimensions
    
    % Calculate the errors for each state
    state_errors = total_state - x_true;
    
    % Initialize the figure
    figure;
    
    % Loop over each state to plot
    for i = 1:6
        sigma = sqrt(squeeze(P_plus(i, i, :)));  % Standard deviations for state i
        
        % Plotting the error with 2-sigma bounds
        subplot(2, 6, i);
        hold on;
        plot(tvec, state_errors(:, i), 'b', 'LineWidth', 1.5);
        plot(tvec, 2*sigma, 'r--', tvec, -2*sigma, 'r--');
        title(['Error for state ', num2str(i)]);
        xlabel('Time');
        ylabel('Error');
        grid on;
        
        % Plotting total state and true state
        subplot(2, 6, i+6);
        hold on;
        plot(tvec, total_state(:, i), 'b', 'LineWidth', 1.5);
        plot(tvec, x_true(:, i), 'g', 'LineWidth', 1.5);
        legend('Total State', 'True State');
        title(['State ', num2str(i), ' Comparison']);
        xlabel('Time');
        ylabel('State value');
        grid on;
    end
    
    % Enhance layout
    sgtitle('State Errors and Comparisons');
end
