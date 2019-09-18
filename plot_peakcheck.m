% Author: Ben Watzak
% July 2019

% This code is meant to accompany the tools adapted from the publication below.
% Roesler, C.J. and K. M. Larson, Software Tools for GNSS Interferometric Reflectometry, 
% GPS Solutions Vol 22:80, doi:10.1007/s10291-018-0744-8, 2018 

% The following plotting commands are meant as a check to make sure peakcheck.m is working
% They are to be used with Newport_Main and plot_type == 0


function plot_peakcheck(maxHeight, GMM)
    xgrid = linspace(0,maxHeight,500)';
    hold on;
    plot(xgrid,pdf(GMM,xgrid));

    hold on;
    %at each respective peak, mark the mixing proportion
    %Mixing proportion j determines the proportion of the population composed by component j,
    pk1 = min(GMM.mu(:)); 
    pk1_ind = find(GMM.mu == pk1);
    mx1 = GMM.ComponentProportion(pk1_ind);

    % now just show it on the plot
    cpk1=sprintf('%4.2f ',pk1);  
    tx = ['Peak 1 H_R: ' cpk1 '(m)'];
    ylim = get(gca,'Ylim');
    text(0.25, ylim(2)-0.25*diff(ylim), tx,'Color','b')
    plot( [pk1 pk1], ylim, 'b--','linewidth',.75)


    cmx1 = sprintf('%4.2f ',mx1);
    tx = ['Peak 1 prop: ' cmx1 ]; 
    text(0.25, ylim(1)+0.25*diff(ylim), tx,'Color','b')
    hold on;
  
    pk2 = max(GMM.mu(:));
    pk2_ind = find(GMM.mu == pk2);
    mx2 = GMM.ComponentProportion(pk2_ind);
    cpk2=sprintf('%4.2f ',pk2);  
    tx = ['Peak 2 H_R: ' cpk2 '(m)'];
    text(maxHeight/2 + 1.5, ylim(2)-0.25*diff(ylim), tx,'Color','b')
    plot( [pk2 pk2], ylim, 'b--','linewidth',.75)

    cmx2 = sprintf('%4.2f ',mx2);
    tx = ['Peak 2 prop: ' cmx2 ]; 
    text(maxHeight/2 + 1.5, ylim(1)+0.25*diff(ylim), tx,'Color','b')
    xlabel('H_R, Reflector Ht. (m)');
    ylabel('pdf value');
  