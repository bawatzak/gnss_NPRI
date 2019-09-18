% Author: Ben Watzak
% July 2019

% This code is meant to accompany the tools adapted from the publication below.
% Roesler, C.J. and K. M. Larson, Software Tools for GNSS Interferometric Reflectometry, 
% GPS Solutions Vol 22:80, doi:10.1007/s10291-018-0744-8, 2018 

% The purpose of this code is to generate a month-long plot of the daily median RH of the Newport site.
% This code asks the user for the year, the month, and the GPS frequency which they would like to use.
% Then it creates a matrix with day of month in column 1, and the median RH for that day in col 2, and creates a plot from this. 

% clear

year = input('4-digit year: ');
month = input('2-digit month: ');
med_avg = input('0 for median, 1 for average: ');

  wave1  = 0.1902936; 

  wave2 = 0.2442102134245683;

  wave5 = 0.254828048;
 


med_monthly_RH = zeros(31,3);
avg_monthly_RH = zeros(31,3);
std_dev = zeros(31,3);
combined_daily = zeros(31,1);
combined_daily_err = zeros(31,1);

l1_avg = 0;


for freqtype = [1 2 5]
    if freqtype == 5
        freq_ind = 3; %column index of the matrix where the frequency's data is stored
    else
        freq_ind = freqtype;
    end
    
  
    
    for day = 1:31
        % run through every day to get that day's median RH
        [today, lsp_amps] = npt_RH(year, month, day, freqtype, true); 
        
         %The following comment block was used when I was originally trying to account for L2 double peaks
         %Now, using the better sp3 files, it's not so much of an issue at NPRI site
        
%         if freqtype ==2
%             %check to make sure there are no double peaks
%             if ~isempty(today)
%                 if max(today) > (wave2/wave1) * monthly_RH(day,1) - std_dev(day, 1) 
%                     %L2 double peaks: one should occur at roughly the same spot as the L1 peak, the other should
%                     %occur at roughly (.2445/.1905) * the L1 peak ... this coefficient is the ratio of wavelengths: L2/L1
%                     % So if max of L2 peaks > (.2445/.1905) * the L1 peak - 1 std dev of L1 peaks, then we're likely dealing with double peaks
%                     fprintf(1,'%s%d\n','You likely have a double peak for L2 on day ',day)
%                     
%                     % now to actually do something about it
%                     % exclude all RH at the second peak
%                     % want to find standard deviation of mean of L1 RH, then set +2 sd as window outside of which to exclude
%                     % max_good = maximum allowed RH = median(L1 RH) + 
%                     max_good =  monthly_RH(day,1) + 2*std_dev(day, 1);
%                     today( today > max_good) = []; %removes values > max_good
%                     L2_new_sd(day) = std(today);
%                 end
%                 if L2_new_sd(day) > 1.5 * std_dev(day, 1) 
%                     fprintf(1,'%s%d\n','You still likely have a double peak for L2 on day ',day)
%                 end
%             end
%         end
        if~isempty(today)
            if med_avg == 0
                %taking medians; we don't care about error
                med_monthly_RH(day,freq_ind) = median(today);
                
            else
                %looking at averages
                if freqtype == 2
                    l1_avg = avg_monthly_RH(day,1);
                    %pass this into peakcheck_auto for comparison purposes
                    std_dev(day,freq_ind) = 0;
                else
                    l1_avg = 0; %not useful if not analyzing L2 peaks
                end
                
                
                %
                %Look to make sure that we can actually take an average
                if length(today) < 3
                    if length(today) < 2
                        %means there's only one satellite track; don't want an average of one thing
                        avg_monthly_RH(day,freq_ind) = NaN;
                        std_dev(day, freq_ind) = NaN;
                    else 
                        avg_monthly_RH(day,freq_ind) = mean(today);
                        std_dev(day, freq_ind) = std(today);
                    end
                    %end
                        
                else
                    %we have enough data to take an average
                    [~, avg_monthly_RH(day,freq_ind),std_dev(day,freq_ind)] = peakcheck_auto(today, lsp_amps, 7, freqtype, l1_avg);
                    %avg_monthly_RH(day,freq_ind) = mean(today);
                    %    std_dev(day, freq_ind) = std(today);
                end
                %
              
            end %Now we're done taking either the median or the average and sd
        end %end of if(~isempty(today))... done with the data set for that day

    end
    %done cycling through every day of the month
end
%done cycling through all three frequencies



if med_avg == 0
    data = med_monthly_RH;
    ylab = sprintf('Median Reflector Height');
else
    data = avg_monthly_RH;
    ylab = sprintf('Average Reflector Height');
end
data(data==0) = NaN;
std_dev(std_dev==0) = NaN;
%remove the days where there either was no data that passed QC, or the month didn't have that day (e.g. February 31st)


for day = 1:31
   combined_daily(day) = nanmean(data(day,:));
   w_i = 1 - (std_dev(day,:) / nansum(std_dev(day,:)));
   combined_daily(day) = nansum(w_i(~isnan(w_i)) * data(day,~isnan(w_i)).') / nansum(w_i);
   
   combined_daily_err(day) = sqrt(  ((w_i(~isnan(w_i)).^2) * (std_dev(day,~isnan(w_i)).^2).') / (nansum(w_i))^2  ) ; %wtd var = {sum (w^2 * var) }/ {(sum(w))^2}
   
end


figure;
hold on;
plot(1:31,data(:,1), 'b o');
hold on;
errorbar(data(:,1),std_dev(:,1),'LineStyle', 'none', 'Color', 'b')
hold on;
plot(1:31,data(:,2), 'r v');
hold on;
errorbar(data(:,2),std_dev(:,2),'LineStyle', 'none', 'Color', 'r')
hold on;
plot(1:31,data(:,3), 'm h');
hold on;
errorbar(data(:,3),std_dev(:,3),'LineStyle', 'none', 'Color', 'm')
hold on;
plot(1:31, combined_daily(:), 'LineStyle', '-', 'Color', 'k', 'LineWidth', 2);

title(['Daily Reflector Height for ' sprintf('%04d',year) '-' sprintf('%02d', month)]);
xlabel('Day of Month')
ylabel(ylab)
legend('L1 Frequency','L1 Error', 'L2 Frequency','L2 Error', 'L5 Frequency','L5 Error', 'Combined Daily Average')
set(gca, 'YDir','reverse')