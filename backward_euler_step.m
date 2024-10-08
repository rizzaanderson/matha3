
function [XB,num_evals] = backward_euler_step(rate_func_in,t,XA,h)
%{
    This function computes the value of X at the next time step
    using the Backward Euler approximation
   
INPUTS:
    rate_func_in: the function used to compute dXdt. rate_func_in will
                  have the form: dXdt = rate_func_in(t,X) (t is before X)
    t: the value of time at the current step
    XA: the value of X(t)
    h: the time increment for a single step i.e. delta_t = t_{n+1} - t_{n}
OUTPUTS:
    XB: the approximate value for X(t+h) (the next step)
        formula depends on the integration method used
    num_evals: A count of the number of times that you called
    rate_func_in when computing the next step
%}  
    G_p = @(XB) rootFunc(rate_func_in,t,XA,h, XB);
    XB = multiNewton(G_p, XA);
    num_evals = 1;
    
end

function G = rootFunc(rate_func_in,t,XA,h, XB)
    G = XA + h*rate_func_in(t+h, XB) - XB;
end