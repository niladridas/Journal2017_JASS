%
function dx = norm_equinoctial_dyn(t,x,norm_matrix,t_norm,equinoctial_dyn)
    % First convert the normalised x and t back to normal x and t
    x_true =  norm_matrix*x;
    t_true =  t_norm*t;
    dx = equinoctial_dyn(t_true,x_true);
    dx = t_norm*((norm_matrix)\dx);
end