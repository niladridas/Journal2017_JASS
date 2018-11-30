function El_val = Lag2El(Lag_val)
    a = Lag_val(1); e = Lag_val(2); i = Lag_val(3); w = Lag_val(4); Ohm = Lag_val(5); nu = Lag_val(6);
    
    p = a*(1-(e^2));
    f = e*cos(w+Ohm);
    g = e*sin(w+Ohm);
    h = tan(i/2)*cos(Ohm);
    k = tan(i/2)*sin(Ohm);
    L = Ohm + w + nu; 
    El_val = [p;f;g;h;k;L];
end