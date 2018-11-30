function [date,dot_n,ddot_n,bstar,inclination,rightascension,...
    essentricity,argperi,meanano,meanmotion, totalrev,fid] = tleparser(fid,FLAG)
% function tle_set = tleparser(tle_line)
% Author: Niladri Das
% Email: niladri@tamu.edu
% Affiliation: Laboratory for Uncertainty Quantification
%              Aerospace Engineering Department, TAMU, TX, USA
% Date: 18th April 2017

% Input: A single two line element
% FLAG: 0 means one read
%       1 keep the file open

% Output: Epoch, n, n_dot, n_ddot, B_star, e, M, i, Ohm, w
% Units will be in radian, mins, kg, ER (earths radius)
% Default values
date = zeros(1,6);
dot_n = 0;
ddot_n = 0;
bstar = 0;
inclination = 0;
rightascension = 0;
essentricity = 0;
argperi = 0;
meanano = 0;
meanmotion = 0;
totalrev = 0;

tline = fgetl(fid);
if ischar(tline)~=0
    % One atomic module has three parts:
    % The name 
    % First Line of TLE, starts with 1
    % Second Line of TLE, starts with 2
    while ischar(tline)
        if str2double(tline(1))==1
            line_1 = tline;
            tline = fgetl(fid);
            line_2 = tline;
            break
        else
         tline = fgetl(fid);
         continue;
        end
    end
    if FLAG == 0
        fclose(fid);
    end
    % We have line 1 and line 2
    % Read the size of each one of them.
    % Sanity check: If the size of each line is not 69,
    % Error messege is "BAD FORMATING OF THE TLE FILE"
    if (length(line_1)~=69 || length(line_2)~=69)
        error('BAD FORMATING OF THE TLE FILE');
    end
    % Main parsing section

    % Step 1: The precious UTC time at which this epoch was generated
    % L1 19-20 : Year
    date = parse_epoch(str2double(line_1(1,19:32)));
    % [Year, Month, Day, Time, Mins, Secs];

    % Step 2: Mean Motion Derivative in 
    % Original Unit: (rev/day^2)/2
    % Transformed Unit: radians/(min^2)
%     date(1,6)
    dotn = str2double(line_1(1,34:43));
    dot_n = 2*dotn*2*pi/((24*60)^2);

    % Step 3: Second Time Derivative of Mean Motion (decimal point assumed) 
    % The last two characters are exponent
    % Original Unit: (rev/day^3)/6
    % Transformed Unit: radians/(min^3)
    % First read the whole set
    ddot_n_read = line_1(1,45:52);
    % Clip until you reach the first character
    while ddot_n_read(1)==' '
        ddot_n_read = ddot_n_read(1,2:end);
    end
    temp_ddotn = strcat('.',ddot_n_read(1,1:end-2));
    temp_ddotn_exp = str2double(ddot_n_read(1,(end-1):end)); % the exponential term
    temp_ddotn = str2double(temp_ddotn);
    ddotn = temp_ddotn*10^(temp_ddotn_exp);
    ddot_n = 6*ddotn*2*pi/((24*60)^3);

    % Step 4:  BSTAR drag term (decimal point assumed)
    % L1 54-61
    % Similar to that if Second time derivative
    % Last two elemets are for eponential part
    % If the SGP4 model was used this field (characters 53 to 60) is the Bstar Drag Parameter or Pseudo Ballistic Coefficient. 
    % Otherwise it is the radiation pressure coefficient.
    % Original Unit : /ER
    % Unit remains unchanged
    % First read the whole set
    bstar_read = line_1(1,54:61);
    % Clip until you reach the first character
    %
    % In some TLE is data set the 54th element is '+' whereas in others it
    % is blank
    %
    while bstar_read(1)==' ' || bstar_read(1)=='+'
        bstar_read = bstar_read(1,2:end);
    end
    temp_bstar = strcat('.',bstar_read(1,1:end-2));
    temp_bstar_exp = str2double(bstar_read(1,(end-1):end)); % the exponential term
    temp_bstar = str2double(temp_bstar);
    bstar = temp_bstar*10^(temp_bstar_exp);

    % Step 5: Ephemeris type
    % L1 63
    % Indicates the model used to calculate the data. 
    % If the SGP or SGP4 models was used (almos allways) this field's 
    % value is 0 (zero). SGP stands for "Simplified General Perturbations",
    % it is the orbital decay model used by the US Space Command 
    % to calculate the orbital elements.
    if str2double(line_1(1,63))~= 0
        error('SGP or SGP4 was not used to generate this TLE');
    end

    %------------------------------------------------------------------%
    % LINE 2 starts
    % Step 6: Inclination in degrees
    % L2 09-16
    % Original Unit: degrees
    % Final Unit: Radians
    inclination_read = line_2(1,9:16);
    inclination = str2double(inclination_read)*pi/180;

    % Step 7: Right Ascension of the Ascending Node [Degrees]
    % L2 18-25
    % Original Unit: degrees
    % Final Unit: Radians
    rightascension_read = line_2(1,18:25);
    rightascension = str2double(rightascension_read)*pi/180;

    % Step 8: Eccentricity (decimal point assumed)
    % L2 27-33
    e_read = line_2(1,27:33);
    temp_e = strcat('.',e_read);
    essentricity = str2double(temp_e);

    % Step 9: Argument of Perigee [Degrees]
    % L2 35-42	
    % Original Unit: degrees
    % Final Unit: Radians
    argperi_read = line_2(1,35:42);
    argperi = str2double(argperi_read)*pi/180;

    % Step 10: Mean Anomaly [Degrees]
    % L2 44-51	
    % Original Unit: degrees
    % Final Unit: Radians
    meanano_read = line_2(1,44:51);
    meanano = str2double(meanano_read)*pi/180;

    % Step 11: Mean Motion [Revs per day]
    % L2 53-63
    % Original Unit: rev/day
    % Final Unit : rad/min
    meanmotion_read = line_2(1,53:63);
    meanmotion = str2double(meanmotion_read)*2*pi/(24*60);

    % Step 12: Revolution number at epoch [Revs]
    % L2 64-68
    % Original Unit: Revolutions
    % Final Unit: Radians
    totalrev_read = line_2(1,64:68);
    totalrev = str2double(totalrev_read)*2*pi;
else
    fid = -1;
end
% % One atomic module has three parts:
% % The name 
% % First Line of TLE, starts with 1
% % Second Line of TLE, starts with 2
% while ischar(tline)
%     if str2double(tline(1))==1
%         line_1 = tline;
%         tline = fgetl(fid);
%         line_2 = tline;
%         break
%     else
%      tline = fgetl(fid);
%      continue;
%     end
% end
% if FLAG == 0
%     fclose(fid);
% end
% % We have line 1 and line 2
% % Read the size of each one of them.
% % Sanity check: If the size of each line is not 69,
% % Error messege is "BAD FORMATING OF THE TLE FILE"
% if (length(line_1)~=69 || length(line_2)~=69)
%     error('BAD FORMATING OF THE TLE FILE');
% end
% % Main parsing section
% 
% % Step 1: The precious UTC time at which this epoch was generated
% % L1 19-20 : Year
% date = parse_epoch(str2double(line_1(1,19:32)));
% % [Year, Month, Day, Time, Mins, Secs];
% 
% % Step 2: Mean Motion Derivative in 
% % Original Unit: (rev/day^2)/2
% % Transformed Unit: radians/(min^2)
% date(1,6)
% dotn = str2double(line_1(1,34:43));
% dot_n = 2*dotn*2*pi/((24*60)^2);
% 
% % Step 3: Second Time Derivative of Mean Motion (decimal point assumed) 
% % The last two characters are exponent
% % Original Unit: (rev/day^3)/6
% % Transformed Unit: radians/(min^3)
% % First read the whole set
% ddot_n_read = line_1(1,45:52);
% % Clip until you reach the first character
% while ddot_n_read(1)==' '
%     ddot_n_read = ddot_n_read(1,2:end);
% end
% temp_ddotn = strcat('.',ddot_n_read(1,1:end-2));
% temp_ddotn_exp = str2double(ddot_n_read(1,(end-1):end)); % the exponential term
% temp_ddotn = str2double(temp_ddotn);
% ddotn = temp_ddotn*10^(temp_ddotn_exp);
% ddot_n = 6*ddotn*2*pi/((24*60)^3);
% 
% % Step 4:  BSTAR drag term (decimal point assumed)
% % L1 54-61
% % Similar to that if Second time derivative
% % Last two elemets are for eponential part
% % If the SGP4 model was used this field (characters 53 to 60) is the Bstar Drag Parameter or Pseudo Ballistic Coefficient. 
% % Otherwise it is the radiation pressure coefficient.
% % Original Unit : /ER
% % Unit remains unchanged
% % First read the whole set
% bstar_read = line_1(1,54:61);
% % Clip until you reach the first character
% while bstar_read(1)==' '
%     bstar_read = bstar_read(1,2:end);
% end
% temp_bstar = strcat('.',bstar_read(1,1:end-2));
% temp_bstar_exp = str2double(bstar_read(1,(end-1):end)); % the exponential term
% temp_bstar = str2double(temp_bstar);
% bstar = temp_bstar*10^(temp_bstar_exp);
% 
% % Step 5: Ephemeris type
% % L1 63
% % Indicates the model used to calculate the data. 
% % If the SGP or SGP4 models was used (almos allways) this field's 
% % value is 0 (zero). SGP stands for "Simplified General Perturbations",
% % it is the orbital decay model used by the US Space Command 
% % to calculate the orbital elements.
% if str2double(line_1(1,63))~= 0
%     error('SGP or SGP4 was not used to generate this TLE');
% end
% 
% %------------------------------------------------------------------%
% % LINE 2 starts
% % Step 6: Inclination in degrees
% % L2 09-16
% % Original Unit: degrees
% % Final Unit: Radians
% inclination_read = line_2(1,9:16);
% inclination = str2double(inclination_read)*pi/180;
% 
% % Step 7: Right Ascension of the Ascending Node [Degrees]
% % L2 18-25
% % Original Unit: degrees
% % Final Unit: Radians
% rightascension_read = line_2(1,18:25);
% rightascension = str2double(rightascension_read)*pi/180;
% 
% % Step 8: Eccentricity (decimal point assumed)
% % L2 27-33
% e_read = line_2(1,27:33);
% temp_e = strcat('.',e_read);
% essentricity = str2double(temp_e);
% 
% % Step 9: Argument of Perigee [Degrees]
% % L2 35-42	
% % Original Unit: degrees
% % Final Unit: Radians
% argperi_read = line_2(1,35:42);
% argperi = str2double(argperi_read)*pi/180;
% 
% % Step 10: Mean Anomaly [Degrees]
% % L2 44-51	
% % Original Unit: degrees
% % Final Unit: Radians
% meanano_read = line_2(1,44:51);
% meanano = str2double(meanano_read)*pi/180;
% 
% % Step 11: Mean Motion [Revs per day]
% % L2 53-63
% % Original Unit: rev/day
% % Final Unit : rad/min
% meanmotion_read = line_2(1,53:63);
% meanmotion = str2double(meanmotion_read)*2*pi/(24*60);
% 
% % Step 12: Revolution number at epoch [Revs]
% % L2 64-68
% % Original Unit: Revolutions
% % Final Unit: Radians
% totalrev_read = line_2(1,64:68);
% totalrev = str2double(totalrev_read)*2*pi;
end








