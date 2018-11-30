% Position and velocity in cartesian coordinates from TLE elements using
% SGP4 model (Spacetrack Report No.3)
% Author: Niladri Das
% Email: niladri@tamu.edu
% Affiliation: Laboratory for Uncertainty Quantification
%              Aerospace Engineering Department, TAMU, TX, USA
% Date: 20 April 2017
% NOTE: the notation are kept as in the report
% Units: Kg, Km, mins

function pos_vel = sgp4(tle_parsed_elements,sgp4_const_elements)
    % Step 1:
    % Input elements from the TLE are first parsed
    % Check for perigee : effects the expression for s_star 
%     PERIGEE_FLAG = 0;
    % Llimit1 : Lower limit settting for solving Kepler Equation in section 10
    Llimit1 = 1e-12;
    % Step 2:
    % Original Mean motion and semimajor axis are first recovered from the
    % Input elements:
    n_o = tle_parsed_elements.n_o;
    e_o = tle_parsed_elements.e_o;
    i_o = tle_parsed_elements.i_o;
    M_o = tle_parsed_elements.M_o;
    w_o = tle_parsed_elements.w_o; 
    Ohm_o = tle_parsed_elements.Ohm_o; 
    B_star = tle_parsed_elements.B_star;
    t_o = tle_parsed_elements.t_o;
    t = tle_parsed_elements.t;
    % (t-t_o) is the time since epoch
    % 1: n_o: SGP type "mean" mean motion at epoch
    %   L2C53-63 rev/day
    % 2: e_o: the “mean” eccentricity at epoch
    %   L2C27-33
    % 3: i_o: the "mean" inclination at epoch
    %   L2C9-16
    % 4: M_o: the "mean" mean anomaly at epoch
    %   L2C44-51 (degrees)
    % 5: w_o: the "mean" argument of perigee at epoch
    %   L2C35-42
    % 6: Ohm_o: the "mean" logitude of ascending node at epoch
    %   L2C18-25
    % 7*: n_o_dot: the time rate of change of “mean” mean motion at epoch
    %   L1C34-43 n_o_dot/2
    % 8*: n_o_ddot: the second time rate of change of “mean” mean motion at epoch
    %   L1C45-52 n_o_ddot/6
    % 9*: B_star: the SGP4 type drag coefficient
    %   L1C54-61
    % 10: k_e: sqrt(GM), where G is Newton’s universal gravitational constant and M is the mass
    % of the Earth
    % 11: a_E:the equatorial radius of the Earth (unit:Km)
    % 12*: J_2: the second gravitational zonal harmonic of the Earth
    % 13*: J_3: the third gravitational zonal harmonic of the Earth
    % 14*: J_4: the fourth gravitational zonal harmonic of the Earth
    %   [J_2,J_3,J_4] = [ 0.0010826269 -0.0000025323 -0.0000016204 ];
    % 15*: (t-t_o): time since epoch
    % 16: k_2: 0.5*J_2*a_E^2
    % 17: k_4: (-3/8)*J_4*a_E^4
    % 18: A_3_0: -J_3*a_E^3
    % 19: q_o: parameter for the SGP4/SGP8 density function
    % 20: s: parameter for the SGP4/SGP8 density function 
    % 21: B: 0.5*C_D*A/m, where the ballistic coefficient for SGP8 where C D is a dimensionless drag
    % coefficient and A is the average cross-sectional area of the satellite of mass m 
    XKMPER = sgp4_const_elements.XKMPER;
    a_E = sgp4_const_elements.a_E;
    k_2 = sgp4_const_elements.k_2;
    k_e = sgp4_const_elements.k_e; % First error here sqrt(GM)
    A_3_0= sgp4_const_elements.A_3_0;
    s = sgp4_const_elements.s;
    qo_s4 = sgp4_const_elements.qo_s4;
    B = B_star;%sgp4_const_elements.B; % Big error
    k_4 = sgp4_const_elements.k_4;
    % Original mean motion n_o_pp (pp: prime prime)
    % Semi-major axis
    a_1 = (k_e/n_o)^(2/3);
    del_1 = (3/2)*(k_2/(a_1^2))*((3*((cos(i_o))^2)-1)/((1-(e_o^2))^(3/2)));
    a_o = a_1*(1-(1/3)*del_1-(del_1^2)-(134/81)*(del_1^3));
    delta_o = (3/2)*(k_2/(a_o^2))*((3*((cos(i_o)^2))-1)/((1-e_o^2)^(3/2)));
    n_o_pp = n_o/(1+delta_o);
    a_o_pp = a_o/(1-delta_o);

    ISIMP = 0;
    if (a_o_pp*(1-e_o)/a_E <= 220/XKMPER+a_E)
        ISIMP = 1;
    end

%   Calculate the Perigee
    Perigee = (a_o_pp*(1-e_o)-a_E)*XKMPER;
    % If Perigee >= 156, no changes needed for s and qo_s4
    if Perigee <= 156 
        s = Perigee - 78;
        if Perigee <= 98
            s = 20;
        end
        qo_s4 = ((120-s)*a_E/XKMPER)^4;
        s = s/XKMPER + a_E;
    end
    
    % Step 5:
    theta = cos(i_o);
    zeta = 1/(a_o_pp-s);
    beta_o = (1-e_o^2)^0.5;
    eta = a_o_pp*e_o*zeta;

    % Step 6:
    C_2 = qo_s4*(zeta^4)*n_o_pp*((1-eta^2)^(-7/2))*(a_o_pp*(1+(3/2)*(eta^2)+4*e_o*eta+e_o*eta^3)...
        + (3/2)*(k_2*zeta/(1-eta^2))*(-0.5+1.5*(theta^2))*(8+24*(eta^2)+3*(eta^4)));
    C_1 = B*C_2;
    C_3 = qo_s4*(zeta^5)*A_3_0*n_o_pp*a_E*sin(i_o)/(k_2*e_o);
    C_4 = 2*n_o_pp*qo_s4*(zeta^4)*a_o_pp*(beta_o^2)*((1-eta^2)^(-7/2))*((2*eta*(1+e_o*eta)+0.5*e_o+0.5*(eta^3))-...
        (2*k_2*zeta/(a_o_pp*(1-eta^2)))*(3*(1-3*(theta^2))*(1+1.5*(eta^2)-2*e_o*eta-0.5*e_o*(eta^3))+...
        (3/4)*(1-theta^2)*(2*eta^2-e_o*eta-e_o*(eta^3))*cos(2*w_o)));
    C_5 = 2*qo_s4*(zeta^4)*a_o_pp*(beta_o^2)*((1-eta^2)^(-7/2))*(1+(11/4)*eta*(eta+e_o)+e_o*(eta^3));

    
    % In the actual code implementation the part independent of D2, D3,D4
    % is done later
    % First
    M_DF = M_o + (1 + (3*k_2*(-1+3*theta^2)/(2*(a_o_pp^2)*beta_o^3)) +...
        (3*k_2^2*(13-78*theta^2+137*theta^4)/(16*(a_o_pp^4)*beta_o^7)) )*n_o_pp*(t-t_o);
    w_DF = w_o + ( -(3*k_2*(1-5*theta^2)/(2*(a_o_pp^2)*beta_o^4)) +...
        (3*k_2^2*(7-114*theta^2+395*theta^4)/(16*(a_o_pp^4)*beta_o^8)) +...
        (5*k_4*(3-36*(theta^2)+49*theta^4)/(4*(a_o_pp^4)*beta_o^8)) )*n_o_pp*(t-t_o);
    Ohm_DF = Ohm_o + (-(3*k_2*theta/(a_o_pp^2*beta_o^4))+...
        (3*k_2^2*(4*theta-19*theta^3)/(2*a_o_pp^4*beta_o^8))+...
        (5*k_4*theta*(3-7*theta^2)/(2*a_o_pp^4*beta_o^8)))*n_o_pp*(t-t_o);
    delw = B*C_3*cos(w_o)*(t-t_o);
    delM = -(2/3)*qo_s4*B*zeta^4*(a_E/(e_o*eta))*((1+eta*cos(M_DF))^3 - (1+eta*cos(M_o))^3);
    M_p = M_DF + delw +delM;
    w = w_DF - delw - delM;
    Ohm = Ohm_DF - (21/2)*(n_o_pp*k_2*theta/(a_o_pp^2*beta_o^2))*C_1*(t-t_o)^2;
    temp_a = 1 - C_1*(t-t_o);
    temp_e = B*C_4*(t-t_o);
    temp_L = (3/2)*C_1*(t-t_o)^2;
%     a = a_o_pp*(1-C_1*(t-t_o)-D_2*(t-t_o)^2-D_3*(t-t_o)^3-D_4*(t-t_o)^4)^2;
    % Step 7:
    if ISIMP==0
        D_2 = 4*a_o_pp*zeta*C_1^2;
        D_3 = (4/3)*a_o_pp*(zeta^2)*(17*a_o_pp+s)*C_1^3;
        D_4 = (2/3)*a_o_pp*(zeta^3)*(221*a_o_pp+31*s)*C_1^4;
%         e = e_o - B*C_4*(t-t_o)-B*C_5*(sin(M_p)-sin(M_o));
        temp_a = (temp_a - D_2*((t-t_o)^2)-D_3*((t-t_o)^3)-D_4*((t-t_o)^4));
        temp_e = temp_e  + B*C_5*(sin(M_p)-sin(M_o));
        temp_L  =  temp_L + (D_2+2*C_1^2)*(t-t_o)^3+...
            (1/4)*(3*D_3+12*C_1*D_2+10*C_1^3)*(t-t_o)^4+...
        (1/5)*(3*D_4+12*C_1*D_3+6*D_2^2+30*C_1^2*D_2+15*C_1^4)*(t-t_o)^5;
    end
    % Step 8:
    % Secular effects of atmospheric drag and gravitation are included through
    % the equations

    a = a_o_pp*(temp_a)^2;
    e =  e_o - temp_e;
    L = M_p + w + Ohm + n_o_pp*temp_L;
    beta = sqrt(1-e^2);
    n = k_e/(a^(3/2));

    % TODO 1:
    % when epoch perigee height is less than 220 kilometers, 
    % the equations for a and L are truncated after the C_1 term, 
    % and the terms involving C_5 , delw, and delM are dropped.

    % Step 9:
    % Adding the long-period periodic terms
    a_xN = e*cos(w);
    L_L = (A_3_0*sin(i_o)/(8*k_2*a*beta^2))*(e*cos(w))*((3+5*theta)/(1+theta));
    a_yNL = A_3_0*sin(i_o)/(4*k_2*a*beta^2);
    L_T = L + L_L;
    a_yN = e*sin(w) + a_yNL;

    % Step 10:
    U = mod(L_T - Ohm,2*pi);
    % Solve Kepler's equation for (E+w) = Ew by using the iteration 
    Ew = U;
    del_Ew = 10; % some random number greater than Llimit1
    while (abs(del_Ew)>=Llimit1)
    del_Ew = (U-a_yN*cos(Ew)+a_xN*sin(Ew)-Ew)/(-a_yN*sin(Ew)-a_xN*cos(Ew)+1);
    Ew = Ew + del_Ew;
    end
%     keyboard
    % Step 11:
    % calculate preliminary quantities needed for short-period periodics
    ecosE = a_xN*cos(Ew) + a_yN*sin(Ew);
    esinE = a_xN*sin(Ew) - a_yN*cos(Ew);
    e_L = (a_xN^2+a_yN^2)^(1/2);
    p_L = a*(1-e_L^2);
    r = a*(1-ecosE);
    r_dot = k_e*sqrt(a)*esinE/r;
    r_fdot = k_e*sqrt(p_L)/r;

    cosu = (a/r)*(cos(Ew)-a_xN+(a_yN*esinE/(1+sqrt(1-e_L^2))));
    sinu = (a/r)*(sin(Ew)-a_yN+(a_xN*esinE/(1+sqrt(1-e_L^2))));
    % TODO 3: source of error
    u = mod(atan2(sinu,cosu),2*pi);
    r_del = (k_2/(2*p_L))*(1-theta^2)*cos(2*u);
    u_del = -(k_2/(4*(p_L^2)))*(7*theta^2-1)*sin(2*u);
    Ohm_del = 3*k_2*theta*sin(2*u)/(2*(p_L^2));
    i_del = 3*k_2*theta*sin(i_o)*cos(2*u)/(2*(p_L^2));
    rdot_del = -(k_2*n/(p_L))*(1-theta^2)*sin(2*u);
    r_del_fdot = (k_2*n/p_L)*((1-theta^2)*cos(2*u)-(3/2)*(1-3*theta^2));

    % Step 12:
    % short-period periodics are added to give the osculating quantities
    r_k = r*(1-(3/2)*k_2*sqrt(1-e_L^2)*(3*theta^2-1)/p_L^2)+r_del;
    u_k = u + u_del;
    Ohm_k = Ohm + Ohm_del;
    i_k = i_o + i_del;
    rdot_k = r_dot + rdot_del;
    r_fdot_k = r_fdot + r_del_fdot;

    % Step 13:
    M_x = -sin(Ohm_k)*cos(i_k);
    M_y = cos(Ohm_k)*cos(i_k);
    M_z = sin(i_k);
    M = [M_x;M_y;M_z];

    N_x = cos(Ohm_k);
    N_y = sin(Ohm_k);
    N_z = 0;
    N = [N_x;N_y;N_z];

    U_vec = M.*sin(u_k) + N.*cos(u_k);
    V_vec = M.*cos(u_k) - N.*sin(u_k);

    % Step 14:
    % POSITION
    r_vec = r_k.*U_vec;
    % Position in Km
    r_vec = XKMPER.*r_vec;
    % VELOCITY
    rdot_vec = rdot_k.*U_vec +  r_fdot_k.*V_vec; 
    pos_vel = [r_vec',rdot_vec'];
    % Velocity in Km/sec
    rdot_vec = rdot_vec.*(XKMPER/60);
    pos_vel = [r_vec',rdot_vec'];
end



