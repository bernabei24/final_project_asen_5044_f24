function plot_states_and_errorsEKF(total_state, x_true, P_plus, tvec)
    % Calculate the total state by adding the deviations from the nominal trajectory to the nominal state
    
    % Calculate the errors for each state
    state_errors = total_state - x_true;
    
    state_errors(:,3) = wrappedAngleDiff(total_state(:,3), x_true(:,3));
    state_errors(:,6) = wrappedAngleDiff(total_state(:,6), x_true(:,6));
    % Initialize the figure
    figure;
    
    % Loop over each state to plot
    for i = 1:6
        sigma = sqrt(squeeze(P_plus(i, i, :)));  % Standard deviations for state i
        
        % Plotting the error with 2-sigma bounds
        subplot(6, 2, 2*i-1);  % Adjusted to be in a 6x2 grid
        hold on;
        p1 = plot(tvec, state_errors(:, i), 'b', 'LineWidth', 1.5);  % State errors
        p2 = plot(tvec, 2*sigma, 'r--', 'LineWidth', 1);            % Positive 2-sigma
        plot(tvec, -2*sigma, 'r--', 'LineWidth', 1);                % Negative 2-sigma
        
        % Custom titles for each state error plot with LaTeX formatting
        switch i
            case 1
                title('$\xi$ (Easting of ground)', 'Interpreter', 'latex');
                ylabel('$\xi$ (m)', 'Interpreter', 'latex');
            case 2
                title('$\eta$ (Northing of ground)', 'Interpreter', 'latex');
                ylabel('$\eta$ (m)', 'Interpreter', 'latex');
            case 3
                title('$\theta$ (Heading of ground)', 'Interpreter', 'latex');
                ylabel('$\theta$ (rad)', 'Interpreter', 'latex');
            case 4
                title('$\xi$ (Easting of air)', 'Interpreter', 'latex');
                ylabel('$\xi$ (m)', 'Interpreter', 'latex');
            case 5
                title('$\eta$ (Northing of air)', 'Interpreter', 'latex');
                ylabel('$\eta$ (m)', 'Interpreter', 'latex');
            case 6
                title('$\theta$ (Heading of air)', 'Interpreter', 'latex');
                ylabel('$\theta$ (rad)', 'Interpreter', 'latex');
        end
        grid on;
        legend([p1 p2], {'State Errors', '$2\sigma$ bounds'}, 'Interpreter', 'latex');

        % Plotting total state and true state
        subplot(6, 2, 2*i);  % Adjusted to be in a 6x2 grid
        hold on;
        plot(tvec, total_state(:, i), 'b', 'LineWidth', 1.5);
        plot(tvec, x_true(:, i), 'g', 'LineWidth', 1.5);
        legend('Total State', 'True State');
        
        % Custom titles for each state comparison plot with LaTeX formatting
        switch i
            case 1
                title('Easting of ground comparison', 'Interpreter', 'latex');
            case 2
                title('Northing of ground comparison', 'Interpreter', 'latex');
            case 3
                title('Heading of ground comparison', 'Interpreter', 'latex');
            case 4
                title('Easting of air comparison', 'Interpreter', 'latex');
            case 5
                title('Northing of air comparison', 'Interpreter', 'latex');
            case 6
                title('Heading of air comparison', 'Interpreter', 'latex');
        end
        
        xlabel('Time');
        ylabel('State value');
        grid on;
    end
    
    % Enhance layout
    sgtitle('State Errors and Comparisons');
end
