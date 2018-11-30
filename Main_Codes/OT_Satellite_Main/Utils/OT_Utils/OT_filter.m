function fltrd_data =  OT_filter(theta_samples,measured_output,cost,weight,OT_constants)
    M = size(theta_samples,2);
    [Aeq,Aeq_1] = OT_constants(M);
    P = Optimal_Transport(theta_samples,measured_output,cost,Aeq,Aeq_1,weight);
    fltrd_data = theta_samples*P;
end