function P_dot = legendre_diff(n,val) % First Derivative
    % Input: n = order of the polynomial
    %        val = value at which the derivative needs to be calculated
    % Legendre Polynomial Derivative
    syms id x
    expr = legendreP(id,x); % i th degree at x
    P_dot = subs(diff(expr,x,1),[id,x],[n,val]); % first derivative at x
end