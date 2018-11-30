function Tsince = time_diff(e1,e2)
% e1 and e2 is the 1X6 element data set representing the date
% date: Year, Month, Day, Hour, Mins, Secs
% Output: Time difference in seconds
% Assumption: e2 > e1
% Measurements taken in the same year
if e1(1,1)~=e2(1,1)
    error('MEASUREMENTS NOT TAKEN IN THE SAME YEAR');
end
% Check if the year is a leap year
if mod(e1(1,1),4)==0
    LEAP_YEAR = 29;
else
    LEAP_YEAR = 28;
end
% For e1 calculate how many seconds has passed since 1st Jan
% 1. Number of months passed e1(1,2)-1
% 2. Number of days passed for these months
i = 1;
days_paste1 = 0;
while i < e1(1,2)
    if mod(i,2)==1
    days_paste1 = days_paste1+31;
    else
        if i == 2
            days_paste1 = days_paste1 + LEAP_YEAR;
        else
            days_paste1 = days_paste1 + 30;
        end
    end
    i = i+1;
end
exact_dayse1 = days_paste1 + e1(1,3);
exact_hourse1 = exact_dayse1*24+e1(1,4);
exact_minse1 = exact_hourse1*60+e1(1,5);
exact_secs1 = exact_minse1*60+e1(1,6);
%%
% For e2 calculate how many seconds has passed since 1st Jan
% 1. Number of months passed e2(1,2)-1
% 2. Number of days passed for these months
i = 1;
days_paste2 = 0;
while i < e2(1,2)
    if mod(i,2)==1
    days_paste2 = days_paste2+31;
    else
        if i == 2
            days_paste2 = days_paste2 + LEAP_YEAR;
        else
            days_paste2 = days_paste2 + 30;
        end
    end
    i = i+1;
end
exact_dayse2 = days_paste2 + e2(1,3);
exact_hourse2 = exact_dayse2*24+e2(1,4);
exact_minse2 = exact_hourse2*60+e2(1,5);
exact_secs2 = exact_minse2*60+e2(1,6);
%%
Tsince = exact_secs2 - exact_secs1; % e2 > e1
end