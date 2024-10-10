function test_euler_fixed_integration()
    tspan = [0, 50];
    h_ref = 0.01;
    X0 = [1;0];
%     t = 0;
%     [XB, ~] = forward_euler_step(@rate_func01, t, X0, h_ref)

    [t_list,X_list,h_avg, num_evals] = backward_euler_fixed_step_integration(@rate_func02,tspan,X0,h_ref);
    
    clf; hold on; 
    title("Rate Func 01: X vs t (for X0)");
    xlabel("time (seconds)"); ylabel("x (meters)");
    plot(t_list, X_list);
    
end