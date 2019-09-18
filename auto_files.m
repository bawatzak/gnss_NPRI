% This is an adaptation of the tools that accompany the publication below.
% Roesler, C.J. and K. M. Larson, Software Tools for GNSS Interferometric Reflectometry, 
% GPS Solutions Vol 22:80, doi:10.1007/s10291-018-0744-8, 2018 
%
% Adapted for NPRI by Ben Watzak
% June 2019
% This was adapted from sample_files.m

function [ filename, outputfile, nofile] = auto_files(year, month, dom, freqtype )
% function [ filename,  outputfile, nofile ] = auto_files( year, month, dom, freqtype)
% input is the year, month, day of the month (dom)
% and the frequency type (freqtype), which can be 1, 2, or 5 
% outputs are filename (SNR)
% outputfile is a txt file to save the LSP outputs


nofile = false;
% dummy values:
outputfile = '';  cyr = ''; 
filename = '';


  station = 'npt';  % Newport
  
  
  if (year - 2000 < 0)
      cyr = sprintf('%04d', year + 2000); %make sure it's a 4-digit year
  else
      cyr = num2str(year);
  end
  
  cmo = sprintf('%02d', month ); %make sure it's a 2-digit month

  cday = sprintf('%02d', dom ); %make sure it's a 2-digit day

% This assumes snr99 files - but you can modify to allow different ones

filename = [station cyr '_' cmo cday '.snr99']; 

% THIS IS HARD-CODED AND WILL NEED TO BE ENTERED ACCORDING TO YOUR FILE
% NAMING CONVENTION AND STORAGE LOCATION
% filename = ['/path/to/yourfile/' filename]; %this is the directory path for the snr file created using sp3 files




if ~exist(filename)
  fprintf(1,'SNR file does not exist for %s', filename)
  fprintf(1,'%s\n', '. Exiting.' )
   outputfile = ''; freqtype = 0; % return dummy values
   nofile = true;
  return
  
else
  % fprintf(1,'SNR data should be in : %s \n', filename);
  outputfile = [station cyr '_' cmo cday '_L' num2str(freqtype) '.txt'];
  % fprintf(1,'LSP Output will go to : %s \n', outputfile);

end
