function midpoint_euler_test()
    t = 0;
    tspan = [0, 100];
    h_ref = 0.1;
    X0 = 0;
    
    % solving step functions for the same point
    [XBm, ~] = explicit_midpoint_step(@rate_func01, t, X0, h_ref);
    [XBe, ~] = forward_euler_step(@rate_func01, t, X0, h_ref);
    
    % solving the analytical solution
    sols = length(tspan);
    ts = linspace(tspan(1), tspan(2));
    for i=1:length(ts)
        sols(i) = solution01(ts(i));
    end
    
    [t_listm, X_listm, ~,~] = explicit_midpoint_fixed_step_integration(@rate_func01, tspan, X0, h_ref);
    [t_liste, X_liste, ~, ~] = forward_euler_fixed_step_integration(@rate_func01, tspan, X0, h_ref);
    clf; hold on;
    title("Rate Func 01: X vs t (for X0)");
    xlabel("time (seconds)"); ylabel("x (meters)");
    plot(t_listm, X_listm);
    hold on
    plot(t_liste, X_liste);
    plot(ts, sols)
    legend("midpoint", "euler", "analytical")
end