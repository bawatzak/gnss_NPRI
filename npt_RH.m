% This is an adaptation of the tools that accompany the publication below.
% Roesler, C.J. and K. M. Larson, Software Tools for GNSS Interferometric Reflectometry, 
% GPS Solutions Vol 22:80, doi:10.1007/s10291-018-0744-8, 2018 
%
% Adapted for NPRI by Ben Watzak
% June 2019
% adapted from sample_gnss_ir.m

function [ all_peaks, lsp_amps] = npt_RH(year,month, dom,freqtype, default )
% Want to functionalize the Newport_Main file to return the mean RH for each day
% Then, in the calling matlab file (npt_monthly_RH.m), it will create a matrix of [dom; mean_RH], then it can plot the full month
% 
% The function takes as input the following:
% year, month, and day of month (dom)
% freqtype: the type of GPS frequency being analyzed (either 1,2, or 5)
% default: a boolean that holds the following significance
%     True: use the default values for QC that I have set for newport
%     False: (not a good idea to use) user must input the QC criteria by hand, for every iteration
%  
% The function returns a vector of the Reflector height values at the LSP peak of every satellite track that passes QC


% the track_maxRH variable will store  a crude median reflector height for a single day/site
track_maxRH = [];

all_peaks = [];
lsp_amps = [];

 % pick up the newport files
    [filename, outputfile, nofile] = auto_files(year, month, dom, freqtype );
    if nofile
      return
    end

%close all
satlist = [1:32 100:150]; % use all GPS satellites and all GLONASS satellites as well
%Roesler & Larson's code defines Satellites as follows:
% GPS 1-99
% Glonass 101-199
% Galileo 201-299
% Beidou 301-399
%There probably aren't 50 glonass satellites 
%but it's a good upper limit so you can be sure that you get every glonass track in the file

LW = 2 ; % linewidth for quadrant plots
comL = '%'; 

if default == false
    disp('This is not set up to enter varying manual QC standards. Try using Newport_Main.m to experiment');
    return
end


pvf = 3; % polynomial order used to remove the direct signal.
% this can be smaller, especially for short elevation angle ranges.


% you can turn off all the plots by changing this to false.
plot2screen = false; 
% THIS IS WHERE FILE PICKUP WAS IN NEWPORT_MAIN... I MOVED IT TO THE CALLING FILE FOR THIS FUNCTION
 

% load the SNR data into variable x, open output file for LSP results
% nofile is a boolean that knows whether data have been found
[fid,x,nofile] = open_filesX(filename,outputfile,freqtype);

if nofile
  return
end


minRH   = 1.2; % meters
maxArcTime = 1; % one hour 

minPoints = 15; 


if default
     
     emin = 3.5; emax = 16; %this is based off analyzing the map of the fresnel zones for the newport site 

       ediff = 6 ;
    % This is the minimum required difference between observed min and max elevation angles, i.e. (maxObsE-minObsE) must be > ediff
 
    maxHeight = 6.3; % (i.e. exclude reflector heights beyond this value)
    desiredPrecision = 0.005; % 5mm 
    frange = [0 6]; % noise range is calculated between 0 and 6 meters.
    az_min = 15;    az_max = 145.1;    %default azimuth limits are the whole 360 degrees of possible azimuths
end


% get additional QC levels and print them to the screen
% You may want to supress the output in quicky_QC once you are sure that
% your data is good
[minAmp,pknoiseCrit,frange ] = quicky_QC(freqtype, maxHeight, desiredPrecision, ediff,frange);

% returns wavelength factor (lambda/2) and column number (using snr format)
[cf,ic] = get_waveL(freqtype);
% make header for the output txt file
output_header( fid );

% you can check azimuths in degree bins, though not necessary because
% you're not plotting anything in this function

az_diff = az_max - az_min;
azrange = 15; 
naz = round(az_diff/azrange); %number of bins

%for a=1:naz
      
      % window by these azimuths
      %azim1 = (a-1)*azrange + az_min;
      %azim2 = azim1 + azrange;

      % window by satellite
      for sat = satlist
          % x contains the SNR data
          % Column 1 satellite number
          % Column 2 elevation angle
          % Column 3 azimuth angle
          % Column 4 second of the day, GPS time (i.e. no leap seconds)
          % Column 7 S1
          % Column 8 S2
          % Column 9 S5
          i=find(x(:,2) < emax & x(:,2) > emin & x(:,1) == sat & x(:,3) > az_min & x(:,3) < az_max); 
          % i creates a single column vector of the row indices where x meets the QC criteria


          if length(i) > minPoints
             w = x(i,:); %creates matrix of only the SNR data in x that met QC criteria
             elevAngles = w(:,2); % vector of elevation angles in degrees
             
             
             data = 10.^(w(:,ic)/20);   % change SNR data from dB-Hz to linear units
             time = w(:,4)/3600;    %3600 sec/hour => time is in now hours instead of sec
             % these are UTC hours  %Coordinated Universal Time (UTC)
             meanUTC = mean(time);
             dt = time(end) - time(1);    % time span of track in hours
             azm = mean(w(:,3));    %average azimuth for a track, in degrees
             
             % remove direct signal.
             pf=polyfit(elevAngles, data, pvf);  
             pv = polyval(pf, elevAngles);
             sineE = sind(elevAngles);  % sin(elevation angles) in degrees
             saveSNR = data-pv; %remove the direct signal with a polynomial
             [sortedX,j] = sort(sineE);

             sortedY = saveSNR(j);
             % ofac (oversampling factor) and hifac (high frequency factor) define the LSP frequency grid spacing and maximum frequency.
             [ofac,hifac] = get_ofac_hifac( elevAngles,cf, maxHeight, desiredPrecision);

             % call the lomb scargle code.  Input data have been scaled so that f comes out in units of reflector height (meters)
             [f,p,dd,dd2]=lomb(sortedX/cf, sortedY, ofac,hifac);
             % returned values are arrays of frequencies considered (f), the associated spectral amplitude(p), ...
                % estimated noise significance of the power values (dd), and the 95% confident level amplitude(dd2).
             [ maxRH, maxRHAmp, pknoise ] = peak2noise(f,p,frange);

             maxObsE = max(elevAngles);
             minObsE = min(elevAngles);   

                if maxRHAmp > minAmp && maxRH > minRH  && dt < maxArcTime && pknoise > pknoiseCrit && (maxObsE-minObsE) > ediff
                    fprintf(fid,'%4.0f %2.0f %2.0f %6.2f %6.2f %6.1f %6.0f %3.0f %6.2f %6.2f %6.2f %4.0f %5.2f\n', ...
                       year,month, dom, maxRH,maxRHAmp,azm, sat, dt*60, minObsE, maxObsE, pknoise,freqtype,meanUTC);
                    %track_maxRH = [track_maxRH; maxRH];
                    
                    %now, look for psd peaks in the satellite track on the LSP
                    if freqtype ==2
                        num_peaks = 1; %used to be 2 
                        %because of potential double peaks, look at the two highest-valued PSD peaks per track
                    else 
                        num_peaks = 1; %because double peaks are not expected
                    end
                    [peaks, peak_RH] = findpeaks(p,f,'NPeaks',num_peaks,'SortStr','descend'); 
                    %returns up to num_peaks peak values of LSP and the RH where they occur
                    %record these values in the daily values vectors
                    all_peaks = [all_peaks; peak_RH ]; 
                    lsp_amps = [lsp_amps; peaks];

                else 
%                   fprintf(1,'%s RH %6.2f Amp %6.2f Azm %6.1f Sat %2.0f Tdiff %4.0f Emin %6.2f Emax %6.2f Peak2Noise %6.2f \n', ...
%                       'Fail QC',maxRH,maxRHAmp,azm, sat, dt*60, minObsE, maxObsE, pknoise);
%                   don't want to print all of this when automating the process
                end % did you pass QC test loop
              
           end %do you have enough points loop
       

      end % satellite loop

     
%end % azimuth loop

%med_RH = median(track_maxRH);



fclose(fid);
