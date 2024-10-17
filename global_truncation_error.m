function global_truncation_error()
    % function that will plot and solve the global error for different
    % values of h_ref, and plot them against each other. Also uses loglog
    % regression to solve for k and p values of error
    tspan = [0, 10];

    % controlled variables for rate_func01
    X01 = 1;

    % controlled variables for rate_func02
    X02 = [1;0];
    
    % set funciton number
    X0 = X01;
    rate_func = @rate_func01;
    X_f = solution01(tspan(2));

    % a vector of 100 equally spaced h_ref values from 10e-5 to 10e1 we
    % will solve our integration functions
    num_h_vals = 50;
    h_refs = logspace(-5, 1, num_h_vals);

    % establishing the size of the vectors we will store the error for a
    % given h value
    exp_midpoint_errors = zeros(length(h_refs),1);
    exp_euler_errors = zeros(length(h_refs),1);
    imp_midpoint_errors = zeros(length(h_refs),1);
    imp_euler_errors = zeros(length(h_refs),1);
    
    total_num_evals = zeros(length(h_refs));

    % for loop to the final X values for the three methods to solve for the
    % global error, epsilon = |X_f - X(t_f)| for different values of h
    for i = 1:length(h_refs)
        % solving for the X_f values
        [~, X_exp_midpoint_list, ~, num_evals_exp_midpoint] = explicit_midpoint_fixed_step_integration(rate_func, tspan, X0, h_refs(i));
        [~, X_exp_euler_list, ~, num_evals_exp_euler] = forward_euler_fixed_step_integration(rate_func, tspan, X0, h_refs(i));
        [~, X_imp_midpoint_list, ~, num_evals_imp_midpoint] = implicit_midpoint_fixed_step_integration(rate_func, tspan, X0, h_refs(i));
        [~, X_imp_euler_list, ~, num_evals_imp_euler] = backward_euler_fixed_step_integration(rate_func, tspan, X0, h_refs(i));
        
        % pulling out the final X value
        X_f_exp_midpoint = X_exp_midpoint_list(:, end);
        X_f_exp_euler = X_exp_euler_list(:, end);
        X_f_imp_midpoint = X_imp_midpoint_list(:, end);
        X_f_imp_euler = X_imp_euler_list(:, end);

        % calcualting the error
        exp_midpoint_errors(i, :) = norm(X_f_exp_midpoint-X_f);
        exp_euler_errors(i, :) = norm(X_f_exp_euler- X_f);
        imp_midpoint_errors(i, :) = norm(X_f_imp_midpoint-X_f);
        imp_euler_errors(i, :) = norm(X_f_imp_euler- X_f);
        
        total_num_evals(i) = num_evals_exp_midpoint;
        
    end
    
    % using provided log regression function to solve for each method's p
    % and k values
    filter_params.min_xval = 10e-4;
    filter_params.max_xval = 10e-2;
    
    [p_exp_midpoint, k_exp_midpoint] = loglog_fit(h_refs, exp_midpoint_errors, filter_params)
    [p_exp_euler, k_exp_euler] = loglog_fit(h_refs, exp_euler_errors, filter_params, filter_params)
    [p_imp_midpoint, k_imp_midpoint] = loglog_fit(h_refs, imp_midpoint_errors, filter_params)
    [p_imp_euler, k_imp_euler] = loglog_fit(h_refs, imp_euler_errors, filter_params)

     % and now solving lines of best fit
    fit_line_x = 10e-5:0.1:10e1;
    fit_exp_line_midpoint = k_exp_midpoint*fit_line_x.^p_exp_midpoint;
    fit_exp_line_euler = k_exp_euler*fit_line_x.^p_exp_euler;

%     Graphs for Assignemnt Day 1
        % plotting errors on a log scale
        clf; figure(1)

        % plotting calculated data
        loglog(h_refs, exp_midpoint_errors, 'ro','markerfacecolor','r','markersize',2)
        hold on;
        loglog(h_refs, exp_euler_errors, 'bo','markerfacecolor','b','markersize',2)

        % plotting lines of best fit
        loglog(fit_line_x, fit_exp_line_midpoint, 'm-')
        loglog(fit_line_x, fit_exp_line_euler, 'c-')
        title("Global Truncation Error for Explicit Methods: Rate Func 1")
        xlabel("h values [-]")
        ylabel("error [-]")
        legend("Explicit Midpoint Data", "Explicit Euler Data", "Midpoint Fit",...
            "Euler Fit", 'location', 'nw')

    % Graphs for Assignment Day 2
        % plotting errors on a log scale
        figure(2);

        % plotting calculated data
        loglog(h_refs, exp_midpoint_errors, 'ro','markerfacecolor','r','markersize',2)
        hold on
        loglog(h_refs, exp_euler_errors, 'bo','markerfacecolor','b','markersize',2)
        loglog(h_refs, imp_midpoint_errors, 'go','markerfacecolor','g','markersize',2)
        loglog(h_refs, imp_euler_errors, 'ko','markerfacecolor','k','markersize',2)

        title("Global Truncation Error Plot: Explicit vs Implicit Methods");
        subtitle("Function 1: 0-10 seconds");
        xlabel("h values [-]")
        ylabel("error [-]")
        legend("Explicit Midpoint", "Explicit Euler", "Implicit Midpoint", "Implicit Euler",...
            'location', 'nw');
        
    % Graphs for Day 3:
        figure(3)
        
        % plotting calculated data
        loglog(total_num_evals(:,1), exp_midpoint_errors, 'ro','markerfacecolor','r','markersize',2)
        hold on;
        loglog(total_num_evals(:,1), exp_euler_errors, 'bo','markerfacecolor','b','markersize',2)

        % plotting lines of best fit
        filter_params.min_xval = 10e2;
        filter_params.max_xval = 10e4;
        
        [p_midpoint_evals_exp, k_midpoint_evals_exp] = loglog_fit(total_num_evals(:,1), exp_midpoint_errors, filter_params)
        [p_euler_evals_exp, k_euler_evals_exp] = loglog_fit(total_num_evals(:,1), exp_euler_errors, filter_params)
        [p_midpoint_evals_imp, k_midpoint_evals_imp] = loglog_fit(total_num_evals(:,1), imp_midpoint_errors, filter_params)
        [p_euler_evals_imp, k_euler_evals_imp] = loglog_fit(total_num_evals(:,1), imp_euler_errors, filter_params)

        fit_line_x = 10e-5:0.1:10e1;
        fit_line_exp_midpoint_evals = k_midpoint_evals_exp*fit_line_x.^p_midpoint_evals_exp;
        fit_line_exp_euler_evals = k_euler_evals_exp*fit_line_x.^p_euler_evals_exp;
        
        loglog(fit_line_x, fit_line_exp_midpoint_evals, 'm-')
        loglog(fit_line_x, fit_line_exp_euler_evals, 'c-')
        
        title("Global Truncation Error Plot vs Number of Function Calls: Explicit Methods");
        xlabel("number of evaluations [-]")
        ylabel("error [-]")
        legend("Explicit Midpoint Data", "Explicit Euler Data", "Midpoint Fit",...
            "Euler Fit", 'location', 'nw')
    
    
        % comparing global errors with number of funciton calls
        figure(4);
    
        % plotting calculated data
        loglog(total_num_evals(:,1), exp_midpoint_errors, 'ro','markerfacecolor','r','markersize',2)
        hold on;
        loglog(total_num_evals(:,1), exp_euler_errors, 'bo','markerfacecolor','b','markersize',2)
        loglog(total_num_evals(:,1), imp_midpoint_errors, 'go','markerfacecolor','g','markersize',2)
        loglog(total_num_evals(:,1), imp_euler_errors, 'ko','markerfacecolor','k','markersize',2)


        title("Global Truncation Error Plot vs Number of Function Calls");
        subtitle("Function 1: 0-10 seconds");
        xlabel("number of evaluations [-]")
        ylabel("error [-]")
        legend("Explicit Midpoint", "Explicit Euler", "Implicit Midpoint", "Implicit Euler",...
            'location', 'nw');
    
end