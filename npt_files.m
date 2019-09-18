% This is an adaptation of the tools that accompany the publication below.
% Roesler, C.J. and K. M. Larson, Software Tools for GNSS Interferometric Reflectometry, 
% GPS Solutions Vol 22:80, doi:10.1007/s10291-018-0744-8, 2018 
%
% Adapted for NPRI by Ben Watzak
% June 2019
% adapted from sample_files.m

function [ filename, year, month, dom, station, outputfile,freqtype,nofile] = npt_files( constellations )
% function [ filename, year, dom, station, outputfile,freqtype,nofile ] = sample_files( )
% user for the tutorial either uses a provided sample file or inputs
% information for their own files.
% outputs are filename (SNR)
% year 
% dom
% outputfile is a txt file to save the LSP outputs
% freqtype is 1,2,or 5
% 

nofile = false;
% dummy values:
outputfile = ''; freqtype = 0; year = 0; dom = 0; station = ''; cyr = ''; 
filename = '';

  % input your own site
  station = input('station: ','s');   %should be 'NPT' for now, may change to lower case 'npt' later (06/14/19)
  year = input('4-digit year: ');
  if (year - 2000 < 0)
      cyr = sprintf('%04d', year + 2000); %make sure it's a 4-digit year
  else
      cyr = num2str(year);
  end
  month = input('2-digit month: ');
  cmo = sprintf('%02d', month ); %make sure it's a 2-digit month
  dom = input('2-digit day of month: ');
  cday = sprintf('%02d', dom ); %make sure it's a 2-digit day
  
% You can add a directory structure if you prefer
% This assumes snr99 files - but you can modify to allow different ones 
filename = [station cyr '_' cmo cday '.snr99']; 

if constellations == 0
   %This is for GPS-only files
   %For our project, this also meant SNR files that were created using nav
   %files (which came before the adoption of sp3)
   %filename = ['/path/to/yourfile/' filename]; 

else
    %filename = ['/path/to/yourfile/' filename]; %this is the directory
    %path for the snr file created using sp3 files and multi-GNSS data
end

if ~exist(filename)
  disp('SNR file does not exist. Exiting.')
   outputfile = ''; freqtype = 0; % return dummy values
   nofile = true;
  return
else
  fprintf(1,'SNR data should be in : %s \n', filename);

  freqtype = input('frequency (1, 2, and 5 are allowed): ');
  if ~ismember(freqtype, [1 2 5])
    disp(['Illegal frequency type: ' num2str(freqtype) ' Exiting'])
    nofile = true;
    return
  end 
  outputfile = [station cyr '_' cmo cday '_L' num2str(freqtype) '.txt'];
  fprintf(1,'LSP Output will go to : %s \n', outputfile);

end
