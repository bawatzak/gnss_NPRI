% This function is used to check the results of GNSS-IR at Newport against 
% the known water levels taken from the co-located tidal gauge.
% First you must import the hourly tidal gauge data as a numerical array 
% (only import the column of verified measurements)

% NPRI data have been downloaded from https://tidesandcurrents.noaa.gov/stationhome.html?id=8452660
% Be sure that the measurements are recorded in Meters! AND with respect to the correct datum!

daily = zeros(31,1);
total = 0; avg = 0;
i = 1; % since matlab starts at index 1
med_daily = zeros(31,1);
month = input('which month do you want to use? (2-digit): ');

% make sure we pick up the right month's data from newport
% You will want to, of course, name your own tidal gauge datasets 
% This is really only practical for the NPRI setting because we had less
% than 1 year of data, so it was easy enough to just do things manually


% if month == 1
%     hourly = TidalJan;
% elseif month == 2
%     hourly = TidalFeb;
% elseif month == 3
%     hourly = TidalMar;
% elseif month == 4
%     hourly = TidalApr;
% elseif month == 5
%     hourly = TidalMay;
% elseif month == 6
%     hourly = TidalJune;
% elseif month == 7
%     hourly = TidalJuly;
% elseif month == 11
%     hourly = TidalNov;
% elseif month == 12
%     hourly = TidalDec;
% 
% else
%     return;
% end

for day = 1:(length(hourly)/24) %i.e. run through only as many days of the month as there are
    total = 0;
    med_daily(day) = median(hourly(i:i+23));
    for hr = 0:23 %this is the timing convention of the tidal gauge (instead of 1:24)
        total = total + hourly(i);
        i = i + 1; 
    end
    avg = total/24;
    daily(day) = avg;
    
end
