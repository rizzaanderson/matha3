function test_euler_fixed_integration()
    t = 0;
    tspan = [0, 50];
    h_ref = 5;
    X0 = 0;
    [XB, ~] = forward_euler_step(@rate_func01, t, X0, h_ref)
    [t_list,X_list,h_avg, num_evals] = forward_euler_fixed_step_integration(@rate_func01,tspan,X0,h_ref);
    
    clf; hold on;
    title("Rate Func 01: X vs t (for X0)");
    xlabel("time (seconds)"); ylabel("x (meters)");
    plot(t_list, X_list);
    
end