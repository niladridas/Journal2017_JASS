function dEl = equinoctial_dyn(t,El,params)
    % Constant Inputs:
    dEl = zeros(6,1);
    mu = params.mu;
    R_e = params.R_e;
    J = params.J; % [0,J2,J3,...];
    % El : Array of Equinoctial Element as defined in the above paper
    % El : {p,f,g,h,k,L}
    p = El(1);f = El(2);g = El(3);h = El(4);k = El(5);L = El(6);
    if(El(1)<0)
        disp('Error1')
    end
    % We need the calculate the following 
    s = sqrt(1 + h^2 + k^2);
    w = 1 + f*cos(L) + g*sin(L);
    r = p/w;
    
    sin_phi = 2*(h*sin(L)-k*cos(L))/(s^2);
    D1 = 0;
    temp_a = nil_legendreP(sin_phi);
    temp_b = nil_legendreP_diff(sin_phi);
    
    for i = 2:length(J) % count starts from 2
        D1 = D1 + (i+1)*J(i)*((R_e/r)^i)*temp_a(1,i);
    end
    C1 = (-mu*cos(L)/(w*r))*D1; 
    C2 = (-mu*sin(L)/(w*r))*D1;
    %
    D2 = 0;
    for i = 2:length(J) % count starts from 2
        D2 = D2 + J(i)*((R_e/r)^i)*temp_b(1,i);
    end
    C3 = (-2*mu/(r*(s^2)))*(h*cos(L)+k*sin(L))*D2 - (mu/(r*w))*(g*cos(L)-f*sin(L))*D1;
    C4 = (mu/(w*(r^2)))*D1;
    C5 = (-2*mu/(r*(s^4)))*((1-h^2+k^2)*sin(L)+2*h*k*cos(L))*D2;
    C6 = (2*mu/(r*(s^4)))*((1+h^2-k^2)*cos(L)+2*h*k*sin(L))*D2;
    
    dEl(1) = 2*sqrt(p/mu)* (-g*C1 + f*C2 + C3);
    dEl(2) = (1/sqrt(mu*p))*(2*p*g*C4 - (1-f^2-g^2)*C2 - (g*(s^2)/2)*(h*C5+k*C6) + ( f + (1+w)*cos(L))*C3);
    dEl(3) = (1/sqrt(mu*p))*(-2*p*f*C4 + (1-f^2-g^2)*C1 + (f*s^2/2)*(h*C5+k*C6) + ( g + (1+w)*sin(L))*C3);
    dEl(4) = ((s^2)/(2*sqrt(mu*p)))*(h* (g*C1 - f*C2 - C3) - ((s^2)/2)*C6);
    dEl(5) = ((s^2)/(2*sqrt(mu*p)))*(k* (g*C1 - f*C2 - C3) + ((s^2)/2)*C5);
    dEl(6) = sqrt(mu*p)*(w/p)^2 + ((s^2)/(2*sqrt(mu*p)))*(h*C5 + k*C6);
end