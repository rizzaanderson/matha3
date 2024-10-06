function [t_list,X_list,h_avg, num_evals] = forward_euler_fixed_step_integration(rate_func_in,tspan,X0,h_ref)
%{
Runs numerical integration using forward Euler approximation
    INPUTS:
        rate_func_in: the function used to compute dXdt. rate_func_in will
            have the form: dXdt = rate_func_in(t,X) (t is before X)
        tspan: a two element vector [t_start,t_end] that denotes the integration endpoints
        X0: the vector describing the initial conditions, X(t_start)
        h_ref: the desired value of the average step size (not the actual value)
    OUTPUTS:
        t_list: the vector of times, [t_start;t_1;t_2;...;.t_end] that X is approximated at
        X_list: the vector of X, [X0’;X1’;X2’;...;(X_end)’] at each time step
        h_avg: the average step size
        num_evals: total number of calls made to rate_func_in during the integration
%}
    
    N = ceil((tspan(2) - tspan(1))/h_ref);
    h_avg = (tspan(2) - tspan(1))/N;
        
    t_list = linspace(tspan(1), tspan(2), N);
    X_list = zeros(length(t_list),length(X0));
    X_list(1,:) = transpose(X0);

    num_evals = 0;
    
    for i = 2:length(t_list)
        [X_next, step_evals] = forward_euler_step(rate_func_in,t_list(i),X_list(i-1),h_avg);
        X_list(i) = X_next;
        num_evals = num_evals + step_evals;
    end
    
end