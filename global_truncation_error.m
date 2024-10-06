function global_truncation_error()
    % function that will plot and solve the global error for different
    % values of h_ref, and plot them against each other. Also uses loglog
    % regression to solve for k and p values of error

    % controlled variables for rate_func01
    tspan1 = [0, 50];
    X01 = 1;

    % controlled variables for rate_func02
    tspan2 = [0, 50];
    X02 = [1;0];
    
    % the analytical X(t_f)
    X_f = solution01(tspan1(2));

    % a vector of 100 equally spaced h_ref values from 10e-5 to 10e1 we
    % will solve our integration functions
    h_refs = logspace(-5, 1, 100);

    % establishing the size of the vectors we will store the error for a
    % given h value
    midpoint_errors = zeros(length(h_refs),1);
    euler_errors = zeros(length(h_refs),1);

    % for loop to the final X values for the three methods to solve for the
    % global error, epsilon = |X_f - X(t_f)| for different values of h
    for i = 1:length(h_refs)
        % solving for the X_f values
        [~, X_midpoint_list, ~, ~] = explicit_midpoint_fixed_step_integration(@rate_func01, tspan1, X01, h_refs(i));
        [~, X_euler_list, ~, ~] = forward_euler_fixed_step_integration(@rate_func01, tspan1, X01, h_refs(i));
        
        % pulling out the final X value
        X_f_midpoint = X_midpoint_list(end);
        X_f_euler = X_euler_list(end);

        % calcualting the error
        midpoint_errors(i) = norm(X_f_midpoint-X_f);
        euler_errors(i) = norm(X_f_euler- X_f);
    end
    
    % using provided log regression function to solve for each method's p
    % and k values
    [p_midpoint, k_midpoint] = loglog_fit(h_refs, midpoint_errors);
    [p_euler, k_euler] = loglog_fit(h_refs, euler_errors);

     % and now solving lines of best fit
    fit_line_x = 10e-5:0.1:10e1;
    fit_line_midpoint = k_midpoint*fit_line_x.^p_midpoint;
    fit_line_euler = k_euler*fit_line_x.^p_euler;

    
    % plotting errors on a log scale
    clf
    
    % plotting calculated data
    loglog(h_refs, midpoint_errors, 'ro','markerfacecolor','r','markersize',1)
    hold on
    loglog(h_refs, euler_errors, 'bo','markerfacecolor','r','markersize',1)
    
    % plotting lines of best fit
    loglog(fit_line_x, fit_line_midpoint, 'c-','markerfacecolor','r','markersize',1)
    loglog(fit_line_x, fit_line_euler, 'm-','markerfacecolor','r','markersize',1)
    title("Error Plots")
    xlabel("h values [-]")
    ylabel("error [-]")
    legend("Midpoint Data", "Euler Data", "Midpoint Fit",...
        "Euler Fit", 'location', 'nw')
end