% This is an adaptation of the tools that accompany the publication below.
% Roesler, C.J. and K. M. Larson, Software Tools for GNSS Interferometric Reflectometry, 
% GPS Solutions Vol 22:80, doi:10.1007/s10291-018-0744-8, 2018 
%
% Adapted for NPRI by Ben Watzak
% June 2019
% adapted from sample_gnss_ir.m
%
% This is a copy of npt_monthly_RH but then modified with another script to get an hourly track of the month
% The purpose of this code is to generate a month-long plot of subdaily average RH of the Newport site.
% This code asks the user for the year, the month, and the GPS frequency which they would like to use.

year = input('4-digit year: ');
month = input('2-digit month: ');
daysInMonth = input ('number of days in this month: ');
% med_avg = input('0 for median, 1 for average: '); %not necessary for this script
% start = input('start day of month: ');
% endD = input('end day of month: ');

  wave1  = 0.1902936; 

  wave2 = 0.2442102134245683;

  wave5 = 0.254828048;
 
full_daily = zeros(daysInMonth*24,3);
errors = zeros(daysInMonth*24,3);
combined_hourly = zeros(daysInMonth*24,1);
combined_hourly_err = zeros(daysInMonth*24,1);
full_daily_1 = []; 
full_daily_2= []; 
full_daily_5= []; 
 errors_1= []; 
 errors_2= []; 
 errors_5 = [];
 sizes_1 = [];
 sizes_2 = [];
 sizes_5 = [];
std_dev = zeros(daysInMonth,3);
combined_daily = zeros(daysInMonth,1);
l1_avg = 0;


for freqtype = [1 2 5]
    if freqtype == 5
        freq_ind = 3; %column index of the matrix where the frequency's data is stored
    else
        freq_ind = freqtype;
    end
    
  
    
    for day = 1:daysInMonth %THIS WILL NEED TO BE MANUALLY ADJUSTED FOR EACH PARTICULAR MONTH
        % run through every day to get that day's median RH
        [today, today_err, today_size] = npt_hourly_RH(year, month, day, freqtype, true); 
        % day_index = ((day -1) * 24) + 1;
        % This is for a sanity check. It just tells you how the hours are
        % indexed in this particular function. i.e. we are not creating 2D
        % arrays with entries (day, hour); rather it is a 1D vector of hourly values.
        
        if freqtype ==1
            full_daily_1= [full_daily_1 ; today];
            errors_1 = [errors_1  ; today_err];
            sizes_1 = [sizes_1; today_size];
        elseif freqtype ==2
            full_daily_2= [full_daily_2 ; today];
            errors_2 = [errors_2  ; today_err];
            sizes_2 = [sizes_2; today_size];
        else
            full_daily_5= [full_daily_5 ; today];
            errors_5 = [errors_5  ; today_err];
            sizes_5 = [sizes_5; today_size];
        end
        
      

    end
    %done cycling through every day of the month that we wanted
end
%done cycling through all three frequencies
full_daily = [full_daily_1, full_daily_2, full_daily_5];

data = full_daily;
std_dev = [errors_1,errors_2,errors_5];

    ylab = sprintf('Average Reflector Height');

data(data==0) = NaN;
sizes_1(sizes_1==0) = NaN;
sizes_2(sizes_2==0) = NaN;
sizes_5(sizes_5==0) = NaN;
%remove the days where there either was no data that passed QC, or the month didn't have that day (e.g. February 31st)

for day = 1:daysInMonth
    for hour = 1:24
        index = (day-1)*24 + hour;
       combined_hourly(index) = nanmean(data(index,:));
       w_i = 1 - (std_dev(index,:) / nansum(std_dev(index,:)));
       combined_hourly(index) = nansum(w_i(~isnan(w_i)) * data(index,~isnan(w_i)).') / nansum(w_i);

       combined_hourly_err(index) = sqrt(  ((w_i(~isnan(w_i)).^2) * (std_dev(index,~isnan(w_i)).^2).') / ((nansum(w_i))^2)  ) ;
       %wtd var = {sum (w^2 * var) }/ {(sum(w))^2}
    end
end






% Create figure
figure1 = figure;

% Create axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');

hold on;
plot(data(:,1), 'b o', 'MarkerFaceColor','b','MarkerSize', 12);
hold on;
errorbar(data(:,1),std_dev(:,1),'LineStyle', 'none', 'Color', 'b')
hold on;
labels = num2str(sizes_1);
%text(1:length(data(:,1)) , data(:,1), labels, 'Color', 'b');
hold on;
plot(data(:,2), 'r v', 'MarkerFaceColor','r','MarkerSize', 12);
hold on;
errorbar(data(:,2),std_dev(:,2),'LineStyle', 'none', 'Color', 'r')
hold on;
labels = num2str(sizes_2);
%text(1:length(data(:,2)) , data(:,2), labels, 'Color', 'r');
hold on;
plot(data(:,3), 'm h', 'MarkerFaceColor','m','MarkerSize', 12);
hold on;
errorbar(data(:,3),std_dev(:,3),'LineStyle', 'none', 'Color', 'm')
hold on;
labels = num2str(sizes_5);
%text(1:length(data(:,3)) , data(:,3), labels, 'Color', 'm');
hold on;
%plot(1:31, combined_daily(:), '-k');

axis(axes1,'ij');
% Set the remaining axes properties
set(axes1,'XGrid','on','XTick',[0:24:(31*24) ],...
    'XTickLabel',...
    { '1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','23','24','25','26','27','28','29','30','31','Next Month'});


title(['Hourly Reflector Height for ' sprintf('%04d',year) '-' sprintf('%02d', month)]);
xlabel('Day of Month')
ylabel(ylab)
legend('L1 Frequency','L1 Error', 'L2 Frequency','L2 Error', 'L5 Frequency','L5 Error')
set(gca, 'YDir','reverse')