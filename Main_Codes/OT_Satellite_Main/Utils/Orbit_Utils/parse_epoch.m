function Date = parse_epoch(epoch_val)
% Author: Niladri Das
% Email: niladri@tamu.edu
% Affiliation: Laboratory for Uncertainty Quantification
%              Aerospace Engineering Department, TAMU, TX, USA
% Date: 18th April 2017

% Input: Epoch from TLE pasrsed 
%       1. 5 digits before decimal
%           1.1 First 2 digits gives the year 19**
%           1.2 Gives the day of the year
%       8 digits after decimal
% Output: [Yr,Mn,Dy,Hr,Mn,Sc]

% Step 1:
% Seperate into integer and decimal
epoch_int = floor(epoch_val); % 5 digits
epoch_dec = epoch_val-epoch_int; % 8 digits
% [year_tmp,day_tmp] = quorem(epoch_int,1000); % Extract the year number and the day
year_tmp = floor(epoch_int/1000);
day_tmp = epoch_int - year_tmp*1000;
Yr = 2000+year_tmp;
% Determine if 'Yr' is a leap year or not
if mod(Yr,4)==0
    LYR_FLAG = 1;
else
    LYR_FLAG =0;
end
% Determine the month and the day
% Initialize day_rem
day_rem = day_tmp;
for i = 1:12
    if mod(i,2)==0 % even month
        if i==2
            day_rem = day_rem - (28+LYR_FLAG);
            if day_rem <= 31 % MArch has 31 days
                Mn = i+1;
                break;
            end
        else
            day_rem = day_rem - 30; 
            if day_rem <= 31 % Odd month has 31 days
                Mn = i+1;
                break;
            end
        end
    else
       day_rem = day_rem - 31; 
       if day_rem <= 30 % Even month has 30 days
            Mn = i+1;
            break;
       end
    end
end
Dy = day_rem;
% Decimal Part
% 
Hr = floor(epoch_dec*24);
dec_Hr = epoch_dec*24 - Hr;
%
Min = floor(dec_Hr*60);
dec_Min = dec_Hr*60 - Min;
%
Sec = dec_Min*60;

Date = [Yr,Mn,Dy,Hr,Min,Sec];
end