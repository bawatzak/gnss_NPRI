% This is an adaptation of the tools that accompany the publication below.
% Roesler, C.J. and K. M. Larson, Software Tools for GNSS Interferometric Reflectometry, 
% GPS Solutions Vol 22:80, doi:10.1007/s10291-018-0744-8, 2018 
%
% Adapted for NPRI by Ben Watzak
% June 2019
% adapted from sample_gnss_ir.m

% For cryosphere applications, please also cite: Larson, K.M., J. Wahr, and P. Kuipers Munneke, Constraints on Snow Accumulation and Firn Density in Greenland 
% Using GPS Receivers, J. Glaciology, Vol. 61, No. 225, doi:10.3189/2015JoG14J130, 2015

constellations = input('0 for GPS only; 1 for multi-GNSS ');

% close all
if constellations == 0
    satlist = 1:32; %GPS only
elseif constellations == 1
    satlist = [1:32 100:150]; % use all GPS satellites, then GLONASS as well
else
    error('You did not enter a proper constellation type');
end

LW = 2 ; % linewidth for quadrant plots
comL = '%'; 


pvf = 3; % polynomial order used to remove the direct signal.
% this can be smaller, especially for short elevation angle ranges.

% the avg_maxRH variable will store reflector height for a single day/site
avg_maxRH = [];

all_peaks = [];
lsp_amps = [];

% you can turn off all the plots by changing this to false.
plot2screen = true; 

% pick up files
[filename, year, month, dom, station, outputfile,freqtype,nofile] = npt_files(1);
if nofile
  return
end
 

% ask if you want l2c satellites  
% note - this just selects L2C satellites. It does not mean your snr file has L2C data in it.
% We want to include L2C for Newport data. 
% However this option excludes all satellites not using L2C

if freqtype == 2
  ok = input('would you like only L2C ? y/n','s');
  if strcmp('y',ok)
    satlist=[1 3 5:10  12 15 17 24:29 30:32 ];
  end
end

% load the SNR data into variable x, open output file for LSP results
% nofile is a boolean that knows whether data have been found
[fid,x,nofile] = open_filesX(filename,outputfile,freqtype);

if nofile
  return
end

% pick the kind of plots you want.
disp('0. Plots for all azimuths together')
disp('1. Separate plots by azimuth bin');
plt_type = input('Plot choice: ');
 

% Rule of thumb: you should not think you are correctly estimating
% RH when it is less than 2*lambda, which is 40-50 cm, depending on the wavelength (l1 vs l2).

minRH   = 1.2; % meters... The water level should never, at least within foreseeable future, get up to this. 
% This minRH represents the distance from GNSS phase center to levelling collar (MWWL).
% They put the MWWL there because water level shouldn't ever get that high

maxArcTime = 1; % one hour 

% Mininum number of points. This value could be used as QC
% it is totally arbitrary for now.  
minPoints = 15; 

% allow users to analyze their own files or vary defaults
% NOT utilized for Newport
idef = input('do you want to use default values for elevation cutoffs, QC, (y/n) ','s');

if strcmp('y', idef)
    
     emin = 3.5; emax = 16; %this is based off the mapping of the fresnel zones for the newport site 

    ediff = 6 ;
    % This is the minimum required difference between observed min and max elevation angles, 
    % i.e. (maxObsE-minObsE) must be > ediff
    % This variable is a QC metric because you don't want a bunch of tiny arcs as these periodograms can be very unreliable.
 
    maxHeight = 6.3; % (i.e. exclude reflector heights beyond this value) 
    desiredPrecision = 0.005; % 5mm 
    frange = [0 6]; % noise range is calculated between 0 and 6 meters.
 
    az_min = 15;    az_max = 145.1;    
    
else %this is for manual entry of the QC parameters
  maxHeight = input('What is largest H_R value (m) you want to estimate? ');
  desiredPrecision=input('Precision of retrieval (m): ');
  emin =input('Minimum elevation angle (degrees): ');
  emax =input('Maximum elevation angle (degrees): ');
  frange =input('Noise calculated over this range (m) e.g. [0 6]: ');
  
  ediff = input('What is the minimum elevation angle difference you require? ');
  if emax < emin
    disp('Emax has to be bigger than Emin.');
    return
  end
  if maxHeight < 0
    disp('Max Reflector Height has to be positive.')
    return
  end
  
  az_min = input('Minimum azimuth angle (degrees): ');
  az_max = input('Maximum azimuth angle (degrees): ');
  if az_max < az_min
    disp('az_max has to be bigger than az_min.');
    return
  end
  if az_max < 0 
      az_max = az_max + 360;
  end
  if az_min < 0 
      az_min = az_min + 360;
  end
      
  
end
figure
% get additional QC levels and print them to the screen
[minAmp,pknoiseCrit,frange ] = quicky_QC(freqtype, maxHeight, desiredPrecision, ediff,frange);

% get wavelength factor (lambda/2) and column where data are stored (ic)
[cf,ic] = get_waveL(freqtype); % returns wavelength factor (lambda/2) and column number (using snr format)
% make header for the output txt file
wavelength = cf*2; %this will be used to catch double peaks

output_header( fid );

% checking azimuths in 15 degree bins
az_diff = az_max - az_min;
azrange = 15; 
naz = round(az_diff/azrange); %number of bins

    for a=1:naz
          if plt_type == 1 && plot2screen 
              % one plot per quadrant in azimuth range    
             figure
          end
          % window by these azimuths
          azim1 = (a-1)*azrange + az_min;
          azim2 = azim1 + azrange;

          % run thru every satellite track
          for sat = satlist
              % x contains the SNR data
              % Column 1 satellite number
              % Column 2 elevation angle
              % Column 3 azimuth angle
              % Column 4 second of the day, GPS time (i.e. no leap seconds)
              % Column 7 S1
              % Column 8 S2
              % Column 9 S5
          
              i=find(x(:,2) < emax & x(:,2) > emin & x(:,1) == sat & x(:,3) > azim1 & x(:,3) < azim2); 
              % i creates a single column vector of the row indices where x meets the QC criteria

              if length(i) > minPoints
                 w = x(i,:); %creates matrix of only the SNR data in x that met QC criteria for the current satellite
                 elevAngles = w(:,2); % vector of elevation angles in degrees


                 data = 10.^(w(:,ic)/20);   % change SNR data from dB-Hz to linear units
                 time = w(:,4)/3600;    %3600 sec/hour => time is in hours instead of sec
                 % these are UTC hours  %Coordinated Universal Time (UTC)
                 meanUTC = mean(time);
                 dt = time(end) - time(1);    % time span of track in hours
                 azm = mean(w(:,3));    %average azimuth for a track, in degrees

                 % remove direct signal. polyfit value does not need to be as large as this for some arcs.
                 pf=polyfit(elevAngles, data, pvf);  
                 
                 pv = polyval(pf, elevAngles);
 
                 sineE = sind(elevAngles);  % sin(elevation angles) in degrees
                 saveSNR = data-pv; %remove the direct signal with a polynomial
                 [sortedX,j] = sort(sineE);
                 % sort the data (in ascending order) so all tracks are rising
                 
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

                    %finally, we actually plot the LSP's
                    if maxRHAmp > minAmp && maxRH > minRH  && dt < maxArcTime && pknoise > pknoiseCrit && (maxObsE-minObsE) > ediff
                         fprintf(fid,'%4.0f %2.0f %2.0f %6.2f %6.2f %6.1f %6.0f %3.0f %6.2f %6.2f %6.2f %4.0f %5.2f\n', ...
                           year,month,dom,maxRH,maxRHAmp,azm, sat, dt*60, minObsE, maxObsE, pknoise,freqtype,meanUTC);
                         if plot2screen
                           if plt_type == 1     %(separate plots)
                              subplot(2,1,1) % raw SNR data
                              plot(asind(sineE), saveSNR, '-','linewidth',LW);  hold on;  %asind gives arcsin(X) in degrees      
                              subplot(2,1,2) % periodogram
                              plot(f,p,'linewidth',LW) ; hold on;
                           else
                               % plot all the periodograms on top of each other in gray
                               subplot(2,1,1)
                               plot(f,p,'color',[0.5 0.5 0.5]); hold on; 
                               % freq units on x-axis will be reflector heights (meters) 
                               % amplitude units on y-axis will be volts/volts  
                          end  
                        end
                        avg_maxRH = [avg_maxRH; maxRH];

                        % now, look for psd peaks in the satellite track on the LSP
                        if freqtype ==2
                           num_peaks = 2;
                            % because of potential double peaks, due to cross correlation, 
                            % look at the two highest-valued PSD peaks per track                       
                        else 
                            num_peaks = 1; %because double peaks are not expected
                        end
                        [peaks, peak_RH] = findpeaks(p,f,'NPeaks',num_peaks,'SortStr','descend'); 
                            %returns up to num_peaks peak values of LSP and the RH where they occur
                            %record these values in the daily values vectors
                        all_peaks = [all_peaks; peak_RH ];
                        lsp_amps = [lsp_amps; peaks];


                    else 
                      fprintf(1,'%s RH %6.2f Amp %6.2f Azm %6.1f Sat %2.0f Tdiff %4.0f Emin %6.2f Emax %6.2f Peak2Noise %6.2f \n', ...
                          'Fail QC',maxRH,maxRHAmp,azm, sat, dt*60, minObsE, maxObsE, pknoise);
                    end % did you pass QC test loop


               end %do you have enough points loop


          end % satellite loop


          if plot2screen 
            plot_labels( plt_type, station, [azim1 azim2],freqtype)
          end
    end % azimuth loop

if plot2screen && plt_type == 0
  plot_labels( plt_type, station, [az_min az_max], freqtype)
  
end

if max(avg_maxRH) > (.2445/.1905)* min(avg_maxRH)
    fprintf('%s\n', 'There is likely a double peak and this is likely L2 data')
end



median_RH( plt_type, avg_maxRH, plot2screen )
hold on;
if plt_type ==0
    
    if length(all_peaks) < 3
        %If we're dealing with a lack of data
        if freqtype ==2
            avg_RH = min(all_peaks);
            %means there's only one satellite track for L2... definitely can't expect to take an average of one thing
            %it's ok for the purposes of this script. We'll just output it to the screen
            fprintf('%s\n','There was only one satellite track. This is NOT an average then');
        else
            if length(all_peaks) < 2
                avg_RH = all_peaks(1);
                %means there's only one satellite track for L1 or L5... definitely can't expect to take an average of one thing
                %it's ok for the purposes of this script. We'll just output it to the screen
                fprintf('%s\n','There was only one satellite track. This is NOT an average then');
            end
        end
        %        
    else
        [gmm, avg_RH, ~] = peakcheck(all_peaks, lsp_amps, maxHeight, freqtype, plot2screen);
    end
    
     
    subplot(2,1,1)
    hold on;
    tx = ['Calculated Average H_R: ' sprintf('%4.2f ',avg_RH) '(m)'];
    ylim = get(gca,'Ylim');
    text(avg_RH + 0.5, ylim(2)-0.5*diff(ylim), tx,'Color','r')
    plot( [avg_RH avg_RH], ylim, 'r--','linewidth',.75)


end
fclose(fid);
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          õ÷MŽ<rDCIÿÿ         +?  % This is an adaptation of the tools that accompany the publication below.
                                                                                                                                                                  