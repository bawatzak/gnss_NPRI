function [doublePk, avg_RH, std_dev] = peakcheck(RH_peaks, LSP_amps, maxHeight, freqtype, l1_avg)
% Author: Ben Watzak
% July 2019

% This code is meant to accompany the tools adapted from the publication below.
% Roesler, C.J. and K. M. Larson, Software Tools for GNSS Interferometric Reflectometry, 
% GPS Solutions Vol 22:80, doi:10.1007/s10291-018-0744-8, 2018 

% This function is used to analyze RH peak values for possible bimodal distribution
% It is set to find a Gaussian mixed distribution model(GMM) with 2 peaks

%fit a bimodal distribution to the peak RH data for the day
doublePk = fitgmdist(RH_peaks, 2,  'RegularizationValue', .0001);
%.Sigma is the covariance matrix (i.e. it returns the variances) 
%so let's get the standard deviations of both peaks
sd1 = sqrt(doublePk.Sigma(:,:,1));
sd2 = sqrt(doublePk.Sigma(:,:,2));

    %at each respective peak, mark the mixing proportion
    %Mixing proportion j determines the proportion of the population composed by component j,
    pk1 = doublePk.mu(1); 
    mx1 = doublePk.ComponentProportion(1);


    cmx1 = sprintf('%4.2f ',mx1);

    pk2 = doublePk.mu(2);
    mx2 = doublePk.ComponentProportion(2);

    cmx2 = sprintf('%4.2f ',mx2);


%Here is where it gets interesting
  %if the two overlap, then they're not statistically different
  %Use the method of assuming unequal variances... I don't know enough yet to assume equal variance of the possibly two peaks


n1 = doublePk.ComponentProportion(1) * length(RH_peaks);
n2 = doublePk.ComponentProportion(2) * length(RH_peaks);

sd1 = sqrt(doublePk.Sigma(:,:,1) *n1 / (n1 - 1));
sd2 = sqrt(doublePk.Sigma(:,:,2) *n2 / (n2 - 1));
df = ( ((sd1^2/n1) + (sd2^2/n2))^2 ) / ( ((sd1^2/n1)^2/(n1-1)) + ((sd2^2/n2)^2/(n2-1))  );
  
tcrit = tinv(.975,df);
%Get a confidence interval for the difference of the two means
ci = [((pk1 - pk2) - tcrit * sqrt((sd1^2/n1) + (sd2^2/n2)) ) , (pk1 - pk2) + tcrit * sqrt((sd1^2/n1) + (sd2^2/n2)) ];
testT = (pk1 - pk2) / sqrt( (sd1^2/n1) + (sd2^2/n2)  ); %test t-value, which can be compared with tcrit as well (instead of seeing if conf. interval contains 0

%if freqtype ~= 2 
    if (n1 <= 1 || n2 <=1)
        %means that it's not even worth investigating; there's at most one point with a peak that isn't with the rest
        %i.e. the data are very highly clustered around a single peak, and double peaks probably aren't a concern
       % fprintf('%s\n', 'There is only one peak!');
        avg_RH = mean(RH_peaks); %just take the average of the peak values for now

        std_dev = std(RH_peaks);
        return
    
    
    elseif max(ci) > 0 && min(ci) < 0
        %Then the difference is not statistically significant 
      %  fprintf('%s\n', 'There is only one peak!');
        avg_RH = mean(RH_peaks); %return the simple average
        std_dev = std(RH_peaks);
        return
    else
       % fprintf('%s\n', 'You have two statistically different peaks. A weighted average is calculated.');
        avg_RH = mx1 * pk1  +  mx2 * pk2 ;
        %weighted average means weighted variance
        
        std_dev = sqrt( (sd1^2)*(mx1^2) + (sd2^2)*(mx2^2) ); % ( Sum of weights )^2 =  ((mx1 + mx2)^2) = 1^2 = 1
        return
    end
    
%end

% THE REST OF THIS IS FOR DEALING WITH L2 DATA, AND HAS BEEN COMMENTED OUT
% BECAUSE THE ISSUE HAS BEEN RESOLVED IN THE NPRI DATA USING SP3 FILES TO
% CREATE SNR DATA FILES. 
% Nonetheless, it remains here for you to use at your site (recommended, at
% least for the first looks at your data), should you need it.

% if max(ci) > 0 && min(ci) < 0
%     %Then the difference is not statistically significant (see may 27 L1 data for example)
%     %If we have taken 2 maxima per track for L2, and STILL don't have different peaks, then this is great,
%     %It would certainly be ok to then remove the smaller peaks, since they're all centered around the same value in this case
%     fprintf('%s\n', 'There is only one peak!');
%     %RETURN THE AVERAGE OF ONLY THE HIGHEST L2 TRACK PEAKS
%     max_RH_peaks = zeros(length(RH_peaks)/2,1);
%     for i = 2:2:length(RH_peaks)
%         %since L2 is told to take the two maximum peaks for each track, we will now tell it to only take one lsp maximum per satellite track
%         j = find( LSP_amps == max([LSP_amps(i) LSP_amps(i-1)]) ); %find the index of the maximum lsp amplitude
%         max_RH_peaks(i/2) = RH_peaks(j);
%     end
%     avg_RH = mean(max_RH_peaks) %now we can return an average of only the highest LSP peaks
%     std_dev = std(RH_peaks);
% 
% else
%     fprintf('%s\n', 'You have two statistically different peaks in this L2 data. Double-peak values are removed to form an average.');
%         %Based off Roesler & Larson 2018, we'll just remove the peak where the RH is greater
%         %"The second peak is located at the L1 peak multiplied by the ratio of the L2 and L1 wavelengths (0.244 and 0.19 m, respectively). 
%         %This second peak is to be expected given that geodetic receivers must cross correlate to extract the L2P data. 
%         %While the secondary peak is straight- forward to observe and exclude, we caution that the peaks will be much more difficult to separate at smaller values of HR. "

%         %This code will effectively remove peak data that are from the second peak. 
%         %So we'll analyze RH peak data that may not necessarily correspond to the highest lsp value, but it will correspond to the correct RH.
%         keep_refining = true;
%         
%         %set maximum RH to be the (ratio of two wavelengths * the smaller peak) - (2 std. dev of the larger peak)
%         min_ind = find( doublePk.mu(:) == min(doublePk.mu(:)) ); %index for the smaller RH peak (most likely ~the same as the L1 peak)
%         max_ind = mod(min_ind, 2) + 1; %i.e. the index of the larger RH-valued peak, which is interpolated from L1
%         mean2 = (0.2442102134245683/0.1902936) * l1_avg %will replace doublePk.mu with L1 mean
% 
%         max_good = mean2 - (2* sqrt(doublePk.Sigma(:,:,max_ind)) )
%         %i.e. we're removing all RH values w/in +/- 2 standard deviations of mean2, as well as any RH greater than 2 std. deviations above mean2
%         RH_peaks
%         while(keep_refining)
%             
%             for i = 1:length(RH_peaks)
%                 if RH_peaks(i) > max_good
%                     RH_peaks(i) = 0; %zero out the bad RH values
%                 end
%             end
%             if max(RH_peaks) < doublePk.mu(min_ind)
%                 fprintf('%s\n', 'You have deleted data that is less than the smaller GMM peak, which may/may not be an issue.');
%                 %SOMETIMES THE LESSER GMM PEAK ISN'T QUITE RIGHT, ESP. WHEN THERE ARE 2 CLEAR PEAKS IN THE ACTUAL LSP, BUT THEY'RE PRETTY CLOSE TOGETHER.
%                 %So deleting values less than this lesser peak won't be a problem when the GMM peaks aren't quite right. 
%                 %However, it might be an issue if the GMM peaks are about right.
%             end
%             
%             if length(find(RH_peaks)) > length(RH_peaks)/2 
%                 % if more nonzeros than zeros, we're still looking at more than one lsp peak per day, which isn't good
%                 % keep_refining = true
%                 max_good = max_good - 0.5*(sqrt(doublePk.Sigma(:,:,max_ind))); %subtract another half standard deviation 
%                 
%             else
%                 keep_refining = false;
%                 %we have zeroed out all lsp peaks that correspond to bad RH values, which may be more than half of the data in some cases.
%                 max_RH_peaks = RH_peaks(find(RH_peaks));
%             end
%             
%             
%         end
%         max_RH_peaks
%         avg_RH = mean(max_RH_peaks)
%         std_dev = std(RH_peaks);
%         %std_dev =sqrt(doublePk.Sigma(:,:,min_ind));
%    %NOTE THAT BECAUSE WE HAVE REMOVED VALUES, THIS ESTIMATE OF THE STD. DEV IS NOT THE BEST
%    %THIS ERROR IS NOT UNBIASED, SINCE WE'VE SYSTEMATICALLY REMOVED POINTS!
% 
% end


                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  