function C = bessel_C(x,N)
% Applying Simpsons Rule
    t = linspace(0,(1-1e-2),N);
    fun2 = @(k) fun1(k,x);
%     x = e(i);
    val_list = fun2(t); % main list
    list1 = val_list(1:(end-1));
    list2 = val_list(2:end);
    arg1 = t(1:(end-1));
    arg2 = t(2:end);
    arg3 = 0.5*(arg1+arg2);
    arg3 = arg3./(ones(1,length(arg3))-arg3);
    list3 = fun2(arg3);
    int_val = ((1)/(N*6))*sum(list1+list2+4*list3);
    C = (int_val);
    function f = fun1(t,x)
    % This function returns the vales evaluated at the given k
    f1 = exp((t./(1-t)).*cos(x));
    f2 = besseli(0,(t./(1-t)));
    f = f1./((1-t).^2.*f2);
    end 
end