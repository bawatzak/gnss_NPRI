Howdy!

The files in this folder are to accompany my undergraduate research project titled:
Using Multi-GNSS Interferometric Reflectometry to Monitor Sea Level in Newport, Rhode Island
The research was completed at the University of Rhode Island, Graduate School of Oceanography, through the SURFO program.
The paper is yet to be published.

This folder contains matlab code, most of which has been adapted from the work of Kristine Larson and Carolyn J. Roesler. 
Their work can and should first be studied through the following publication:
Roesler, C.J. and K. M. Larson, Software Tools for GNSS Interferometric Reflectometry, GPS Solutions, 
Vol 22:80, doi:10.1007/s10291-018-0744-8, 2018 

You are recommended to fist begin with their code to familiarize yourself with GNSS-IR, which is available at:
https://github.com/kristinemlarson 
In fact, in order to use any of the code that is published here, you must first have the rest of the toolkit found in the
above link. This gnss_NPRI folder contains only the files that I have created or adapted to go along with the method.

-------------------------------------------------------------------------------------------------------------------

You are urged to begin with Newport_Main.m, which will look at a single SNR file at a time. Once you are familiar with the basics
of the methods used in this project, you should then work your way into npt_monthly_RH.m and the other files which automate the
process to look at RH over time. You should also note that npt_monthly_RH.m and the functions called within are the ones that were
used for the final results of my work. The other files in this folder should only be used for preliminary and cursory analysis.
For example, peakcheck.m is not the final version that was used, but peakcheck_auto.m is.

The particular matlab files found here are adapted so as to implement a new method of taking an average over the course of a day.
This method makes use of a Gaussian Mixture Model and weighted averages. This code is being published with the intent that it is 
implemented in other locations, and hopefully improved upon. File locations are hard-coded, so this will need to change based on 
your file locations and naming conventions. 

I am well aware that these algorithms are not the most efficient, but they do work. Math is my area of expertise, not computer 
science. If somebody is interested in making these scripts more time and computationally efficient, especially for larger data 
sets, that is highly encouraged and would be appreciated.


Should you have any questions, you should first contact my advisor, Dr. Meng "Matt" Wei. He will be continuing with the project.
matt-wei@uri.edu


Ben Watzak
Texas A&M University - College Station
September 17, 2019