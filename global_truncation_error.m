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
    num_h_vals = 100;
    h_refs = logspace(-5, 1, num_h_vals);

    % establishing the size of the vectors we will store the error for a
    % given h value
    exp_midpoint_errors = zeros(length(h_refs),1);
    exp_euler_errors = zeros(length(h_refs),1);
    imp_midpoint_errors = zeros(length(h_refs),1);
    imp_euler_errors = zeros(length(h_refs),1);

    % for loop to the final X values for the three methods to solve for the
    % global error, epsilon = |X_f - X(t_f)| for different values of h
    for i = 1:length(h_refs)
        % solving for the X_f values
        [~, X_exp_midpoint_list, ~, num_evals_exp_midpoint] = explicit_midpoint_fixed_step_integration(@rate_func01, tspan1, X01, h_refs(i));
        [~, X_exp_euler_list, ~, num_evals_exp_euler] = forward_euler_fixed_step_integration(@rate_func01, tspan1, X01, h_refs(i));
        [~, X_imp_midpoint_list, ~, num_evals_imp_midpoint] = implicit_midpoint_fixed_step_integration(@rate_func01, tspan1, X01, h_refs(i));
        [~, X_imp_euler_list, ~, num_evals_imp_euler] = backward_euler_fixed_step_integration(@rate_func01, tspan1, X01, h_refs(i));
        
        % pulling out the final X value
        X_f_exp_midpoint = X_exp_midpoint_list(end);
        X_f_exp_euler = X_exp_euler_list(end);
        X_f_imp_midpoint = X_imp_midpoint_list(end);
        X_f_imp_euler = X_imp_euler_list(end);

        % calcualting the error
        exp_midpoint_errors(i) = norm(X_f_exp_midpoint-X_f);
        exp_euler_errors(i) = norm(X_f_exp_euler- X_f);
        imp_midpoint_errors(i) = norm(X_f_imp_midpoint-X_f);
        imp_euler_errors(i) = norm(X_f_imp_euler- X_f);
    end
    
    % using provided log regression function to solve for each method's p
    % and k values
    [p_exp_midpoint, k_exp_midpoint] = loglog_fit(h_refs, exp_midpoint_errors);
    [p_exp_euler, exp_k_euler] = loglog_fit(h_refs, exp_euler_errors);
    [p_imp_midpoint, k_imp_midpoint] = loglog_fit(h_refs, imp_midpoint_errors);
    [p_imp_euler, imp_k_euler] = loglog_fit(h_refs, imp_euler_errors);

     % and now solving lines of best fit
%     fit_line_x = 10e-5:0.1:10e1;
%     fit_exp_line_midpoint = k_exp_midpoint*fit_line_x.^p_exp_midpoint;
%     fit_exp_line_euler = exp_k_euler*fit_line_x.^p_exp_euler;

    % Graphs for Assignemnt Day 1
%         % plotting errors on a log scale
%         clf
% 
%         % plotting calculated data
%         loglog(h_refs, midpoint_errors, 'ro','markerfacecolor','r','markersize',1)
%         hold on
%         loglog(h_refs, euler_errors, 'bo','markerfacecolor','r','markersize',1)
% 
%         % plotting lines of best fit
%         loglog(fit_line_x, fit_line_midpoint, 'c-','markerfacecolor','r','markersize',1)
%         loglog(fit_line_x, fit_line_euler, 'm-','markerfacecolor','r','markersize',1)
%         title("Error Plots")
%         xlabel("h values [-]")
%         ylabel("error [-]")
%         legend("Midpoint Data", "Euler Data", "Midpoint Fit",...
%             "Euler Fit", 'location', 'nw')

    % Graphs for Assignment Day 2
        % plotting errors on a log scale
        clf; figure(1);

        % plotting calculated data
        loglog(h_refs, exp_midpoint_errors, 'ro','markerfacecolor','r','markersize',3)
        hold on
        loglog(h_refs, exp_euler_errors, 'bo','markerfacecolor','b','markersize',3)
        loglog(h_refs, imp_euler_errors, 'go','markerfacecolor','g','markersize',3)
        loglog(h_refs, imp_euler_errors, 'ko','markerfacecolor','k','markersize',3)

        title("Global Truncation Error Plot: Explicit vs Implicit Methods");
        xlabel("h values [-]")
        ylabel("error [-]")
        legend("Explicit Midpoint", "Explicit Euler", "Implicit Midpoint", "Implicit Euler",...
            "Analytical", 'location', 'nw');
        
        % plotting errors on a log scale
        figure(2);

        % plotting calculated data
        loglog(h_refs, exp_midpoint_errors*num_evals_exp_midpoint, 'ro','markerfacecolor','r','markersize',3)
        hold on
        loglog(h_refs, exp_euler_errors*num_evals_exp_euler, 'bo','markerfacecolor','b','markersize',3)
        loglog(h_refs, imp_midpoint_errors*num_evals_imp_midpoint, 'go','markerfacecolor','g','markersize',3)
        loglog(h_refs, imp_euler_errors*num_evals_imp_euler, 'ko','markerfacecolor','k','markersize',3)

        title("Scaled Global Truncation Error Plot");
        xlabel("h values [-]")
        ylabel("error [-]")
        legend("Explicit Midpoint", "Explicit Euler", "Implicit Midpoint", "Implicit Euler",...
            "Analytical", 'location', 'nw');
        
    
end