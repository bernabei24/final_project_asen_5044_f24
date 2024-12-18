function plotInnovationsAndMeasurements(dely_diff, dely, dely_predict)
    % Validate input dimensions
    if size(dely_diff, 1) ~= 5 || size(dely_diff, 2) ~= 1001 || ...
       size(dely, 1) ~= 1001 || size(dely, 2) ~= 5 || ...
       size(dely_predict, 1) ~= 5 || size(dely_predict, 2) ~= 1001
        error('Incorrect matrix dimensions.');
    end

    % Create a single figure for all plots
    figure;

    % Iterate through each measurement
    for i = 1:5
        % Plot innovations in the 1st column of the ith row
        subplot(5, 2, 2*i-1); % Positions subplots in columns 1 of rows 1 to 5
        plot(1:1001, dely_diff(i, :), 'b');
        title(sprintf('Innovations for Measurement %d', i));
        xlabel('Time Step');
        ylabel('Innovation');
        grid on; % Optionally add a grid

        % Plot true and predicted measurements in the 2nd column of the ith row
        subplot(5, 2, 2*i); % Positions subplots in columns 2 of rows 1 to 5
        plot(1:1001, dely(:, i), 'b-', 1:1001, dely_predict(i, :), 'r--');
        legend('True', 'Predicted');
        title(sprintf('True vs. Predicted for Measurement %d', i));
        xlabel('Time Step');
        ylabel('Measurement Value');
        grid on; % Optionally add a grid
    end
end
