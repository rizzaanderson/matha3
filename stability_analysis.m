function stability_analysis()
    % function that compares the stability of back/forward euler and
    % im/explicity midpoint for eq 5, rate_func01, to the numerical
    % solution

    % recommended conditions
    tspan = [0, 20];
    h38 = 0.38;
    h45 = 0.45;
    X0 = 1;

    % solving for the numerical solution
    t_sol = linspace(0, 20);
    X_sol = zeros(length(t_sol),1);
    
    for i = 1:length(t_sol)
        X_sol(i) = solution01(t_sol(i));
    end

    % calculating integration using the euler methods for h_ref = 0.38
    [forward_t38, forward_X38, ~, ~] = forward_euler_fixed_step_integration...
        (@rate_func01, tspan, X0, h38);
    [backward_t38, backward_X38, ~, ~] = backward_euler_fixed_step_integration...
         (@rate_func01, tspan, X0, h38);
    
    
     % calculating integration using the midpoint methods for h_ref = 0.38
     [immid_t38, immid_X38, ~, ~] = implicit_midpoint_fixed_step_integration...
         (@rate_func01, tspan, X0, h38);
     [exmid_t38, exmid_X38, ~, ~] = explicit_midpoint_fixed_step_integration...
         (@rate_func01, tspan, X0, h38);

     % plotting each method against the solution
     tile_plot38 = tiledlayout(2,2);
     title(tile_plot38, 'Comparing Stabiliy for h_{ref} = 0.38')
     xlabel(tile_plot38, 't [s]')
     ylabel(tile_plot38, 'X [m]')
 
     
     nexttile
     axis([0 20 -1.5 1.5])
     plot(t_sol, X_sol)
     hold on
     plot(forward_t38, forward_X38)
     hold off
     title('Foward Euler')
     legend('numerical', 'analytical', 'location', 'nw')

     nexttile
     axis([0 20 -1.5 1.5])
     plot(t_sol, X_sol)
     hold on
     plot(backward_t38, backward_X38)
     hold off
     title('Backward Euler')
     legend('numerical', 'analytical', 'location', 'nw')

     nexttile
     axis([0 20 -1.5 1.5])
     plot(t_sol, X_sol)
     hold on
     plot(exmid_t38, exmid_X38)
     hold off
     title('Explicit Midpoint')
     legend('numerical', 'analytical', 'location', 'nw')

     nexttile
     axis([0 20 -1.5 1.5])
     plot(t_sol, X_sol)
     hold on
     plot(immid_t38, immid_X38)
     hold off
     title('Implicit Midpoint')
     legend('numerical', 'analytical', 'location', 'nw')


     % calculating integration using the euler methods for h_ref = 0.45
    [forward_t45, forward_X45, ~, ~] = forward_euler_fixed_step_integration...
        (@rate_func01, tspan, X0, h45);
    [backward_t45, backward_X45, ~, ~] = backward_euler_fixed_step_integration...
         (@rate_func01, tspan, X0, h45);


     % calculating integration using the midpoint methods for h_ref = 0.45
     [immid_t45, immid_X45, ~, ~] = implicit_midpoint_fixed_step_integration...
         (@rate_func01, tspan, X0, h45);
     [exmid_t45, exmid_X45, ~, ~] = explicit_midpoint_fixed_step_integration...
         (@rate_func01, tspan, X0, h45);

     % plotting each method against the solution
     tile_plot45 = tiledlayout(2,2);
     title(tile_plot45, 'Comparing Stabiliy for h_{ref} = 0.45')
     xlabel(tile_plot45, 't [s]')
     ylabel(tile_plot45, 'X [m]')


     nexttile
     axis([0 20 -1.5 1.5])
     plot(t_sol, X_sol)
     hold on
     plot(forward_t45, forward_X45)
     hold off
     title('Foward Euler')
     legend('numerical', 'analytical', 'location', 'nw')

     nexttile
     axis([0 20 -1.5 1.5])
     plot(t_sol, X_sol)
     hold on
     plot(backward_t45, backward_X45)
     hold off
     title('Backward Euler')
     legend('numerical', 'analytical', 'location', 'nw')

     nexttile
     axis([0 20 -1.5 1.5])
     plot(t_sol, X_sol)
     hold on
     plot(exmid_t45, exmid_X45)
     hold off
     title('Explicit Midpoint')
     legend('numerical', 'analytical', 'location', 'nw')

     nexttile
     axis([0 20 -1.5 1.5])
     plot(t_sol, X_sol)
     hold on
     plot(immid_t45, immid_X45)
     hold off
     title('Implicit Midpoint')
     legend('numerical', 'analytical', 'location', 'nw')
end