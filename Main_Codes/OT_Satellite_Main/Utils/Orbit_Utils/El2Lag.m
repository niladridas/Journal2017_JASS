function Lag_element  = El2Lag(El)
    p = El(1); f = El(2); g = El(3); h = El(4); k = El(5); L = El(6) ;
    
    e = sqrt(f^2+g^2);
    a = p/(1-e^2);
    Ohm = atan2(k,h);
    w = atan2(g,f)-Ohm;
    i = 2*atan(h/cos(Ohm));
    nu = L-w-Ohm;
    Lag_element = [a;e;i;w;Ohm;nu];
end