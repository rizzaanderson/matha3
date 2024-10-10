function implicit_solver()
    % function that solves implicit_midpoint_step and
    % implicit_midpoint_fixed_step_integration. Plots
    % implicit_midpoint_fixed_step_integrations for different values of h
    
    % constants
    t = 0;
    tspan = [0, 2*pi];
    h_ref = 0.1;
    X0 = 1;
    
    % solving step functions for the same point
    [XB_mid, ~] = implicit_midpoint_step(@rate_func01, t, X0, h_ref);
    [XB_eul, ~] = backward_euler_step(@rate_func01, t, X0, h_ref);
    
    % solving over a time span
    [t_list,X_list,h_avg, num_evals] = ...
    implicit_midpoint_fixed_step_integration(@rate_func02,tspan,X0,h_ref);

    % looping through different values of h
    h = linspace(10e-5, 10e-2, 7)

    % setting up plot for multiple h values
    clf;
    
    legendStrings = "h = " + string(h)
    for i = 1:length(h)
        
        [t_plot, X_plot, ~, ~] = implicit_midpoint_fixed_step_integration(@rate_func01,tspan,X0,h(i));
        hold on;
        plot(t_plot, X_plot)
        title("Rate Func 01: X vs t (for X0)");
        xlabel("time (seconds)"); ylabel("x (meters)");
        legend(legendStrings)
    end
    
    
end
