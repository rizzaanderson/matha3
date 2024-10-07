function local_truncation_error()
    % function that will plot and solve the local error for different
    % values of h_ref, and plot them against each other. Also uses loglog
    % regression to solve for k and p values of error

    % controlled variables for rate_func01
    t1 = 5;
    X01 = 3;

    % solve for analytical solution for rate_func01
    X01_sol = solution01(t1);

    % controlled variables for rate_func02
    t2 = 5;
    X02 = [1;0];


    % a vector of 100 equally spaced h_ref values from 10e-5 to 10e1 we
    % will solve our integration functions
    h_refs = logspace(-5, 1, 100);

    % establishing the size of the vectors we will store the error for a
    % given h value
    midpoint_errors = zeros(length(h_refs),1);
    euler_errors = zeros(length(h_refs),1);
    X_sol = zeros(length(h_refs), 1);
    analytical = zeros(length(h_refs), 1);

    % for loop to solve for X_n+1 and X(t+h) with midpoint, euler, and
    % analytical methods to calculate the error associated with that 
    % particular h value for epsilon = |G(t, X(t), h) - X(t+h)|
    for i = 1:length(h_refs)
        % X_next with the midpoint and euer method
        [X_next_midpoint, ~] = explicit_midpoint_step(@rate_func01, t1, X01_sol, h_refs(i));
        [X_next_euler, ~] = forward_euler_step(@rate_func01, t1, X01_sol, h_refs(i));
        
        % solving for the analytical solution
        X_sol(i) = solution01(t1 + h_refs(i));
        
        % solving for the error of the three methods
        midpoint_errors(i) = norm(X_next_midpoint - X_sol(i));
        euler_errors(i) = norm(X_next_euler - X_sol(i));
        analytical(i) = norm(X_sol(i) - X01_sol); 
            % ^^ shouldn't this error always be zero because this is what we're comparing to?
    end
        
    % using provided log regression function to solve for each method's p
    % and k values
    [p_midpoint, k_midpoint] = loglog_fit(h_refs, midpoint_errors);
    [p_euler, k_euler] = loglog_fit(h_refs, euler_errors);
    [p_analytical, k_analytical] = loglog_fit(h_refs, analytical);

     % and now solving lines of best fit
    fit_line_x = 10e-7:0.1:10e1;
    fit_line_midpoint = k_midpoint*fit_line_x.^p_midpoint;
    fit_line_euler = k_euler*fit_line_x.^p_euler;
    fit_line_analytical = k_analytical*fit_line_x.^p_analytical;

    
    % plotting errors on a log scale
    clf;
    
    % plotting calculated data
    loglog(h_refs, midpoint_errors, 'ro','markerfacecolor','r','markersize',2)
    hold on
    loglog(h_refs, euler_errors, 'bo','markerfacecolor','r','markersize',2)
    loglog(h_refs, analytical, 'go', 'markerfacecolor', 'g', 'markersize', 2)
    
    % plotting lines of best fit
    loglog(fit_line_x, fit_line_midpoint, 'c-','markerfacecolor','r','markersize',1)
    loglog(fit_line_x, fit_line_euler, 'm-','markerfacecolor','r','markersize',1)
    loglog(fit_line_x, fit_line_analytical, 'y-','markerfacecolor','r','markersize',1)
    title("Error Plots")
    xlabel("h values [-]")
    ylabel("error [-]")
    legend("Midpoint Data", "Euler Data", "Analytical Data", "Midpoint Fit",...
        "Euler Fit", "Analytical Fit", 'location', 'nw')
    
    
end