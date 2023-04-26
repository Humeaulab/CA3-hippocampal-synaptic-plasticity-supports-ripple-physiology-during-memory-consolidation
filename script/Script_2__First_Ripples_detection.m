% function [ripples ripplesint rint freqprof] = discrimSWR(x,stratumR,stratumO,SR)
% close all
clear all
warning off

BaseDirectory='';
MatlabDirectory= '\\filer2-iins\HUMEAU\EVERYONE\DATA IN VIVO\Marilyn-Matlab'; % MATLAB directory containing the scripts for ripples analysis.

cd(BaseDirectory);
load('Info.mat')

RipplesDataDirectory = uigetdir(DataExtracted,'select the Ripples Folder');

cd(BaseDirectory)
if exist('Info.mat') == 0
    save('Info.mat', 'RipplesDataDirectory')
else
    save('Info.mat', 'RipplesDataDirectory','-append')
end



SR=1000;
SRi=1000;

rippleband=[100 250];

smoothF=round(SR/(rippleband(1)));
dSR=50;

i = 1;
SWsign=0;
counter = 1;

%% Ripples detection: pic, start, end on referenced traces.

cd(BaseDirectory)


for j = 1:length (g)
    cd(LFPsDirectory)
    load('HPC_Dataref_nex.mat')
    
    cd(MatlabDirectory)
    
    x00 = eval(sprintf('CH%d_HPCrds',g(j)));
    timeStampcolumn=x00(:,2);
    t00=(0:1/SR:length(x00)/SR-1/SR)'; 
    
    xD00 = BandPass((x00(:,1)),SR,rippleband(1),rippleband(2),3); 
    xD00(:,2)=timeStampcolumn;
    assignin('base',sprintf('CH%d_HPCRipplesR',g(j)),xD00);
    
    xrms00=smooth(abs(xD00(:,1)),smoothF,'moving'); 
    xrms00(:,2)=timeStampcolumn;
    assignin('base',sprintf('CH%d_HPCRipplesAbsR',g(j)),xrms00);
    
    Std = std(xrms00(:,1));
    MeanR=mean(xrms00(:,1));
    Threshold= MeanR+(5*Std)
    assignin('base',sprintf('CH%d_HPCThreshold',g(j)),Threshold);
    
    [pksr00,locR00,widthr00,promr00] = findpeaks(xrms00(:,1),'minpeakheight',(MeanR+5*Std),'MinPeakWidth',0.015*SR,...
        'MinPeakDistance',0.045*SR);
    LocR00=x00(locR00,2);
    
    assignin('base',sprintf('CH%d_HPCLocRipples',g(j)),LocR00);
    assignin('base',sprintf('CH%d_HPCHalfProminence',g(j)),widthr00);
    
    xrmsONOFF=xrms00-(2*Std);
    xrmsONOFF(:,2)=timeStampcolumn;
    assignin('base',sprintf('CH%d_HPCAbsDownR',g(j)),xrmsONOFF);
    matriceRipples=zeros(length(locR00),3); 
    matriceRipples(:,3)=locR00; 
    
    tic
    for i=1:size(matriceRipples,1)
        a=find(xrmsONOFF(1:locR00(i),1)<0);
        if isempty(a)==1
            matriceRipples(i,1)=0;
        else
            matriceRipples(i,1)=(a(end));
        end
        
        a=find(xrmsONOFF(locR00(i):end,1)<0);
        if isempty(a)==1
            matriceRipples(i,2)=x00(end);
        else
            matriceRipples(i,2)=(a(1))+matriceRipples(i,3);
        end
        if i<length(matriceRipples)
            if matriceRipples((i)+1,1)<=matriceRipples((i),2)
                matriceRipples((i)+1,1)=matriceRipples((i),2)+1;
            end
        end
    end
    
    toc
    
    conflictSW=find(diff(matriceRipples(:,1))<=0);
    
    for i=1:length(conflictSW)
        localminSW=find(xrmsONOFF(matriceRipples(conflictSW(i),3):matriceRipples(conflictSW(i)+1,3),1)==...
            min(xrmsONOFF(matriceRipples(conflictSW(i),3):matriceRipples(conflictSW(i)+1,3),1)));
        matriceRipples(conflictSW(i),2)=localminSW+matriceRipples(conflictSW(i),3)-1;
        matriceRipples(conflictSW(i)+1,1)=localminSW+matriceRipples(conflictSW(i),3)+1;
    end
    
    for i=1:length(matriceRipples)-1
        if matriceRipples(i,2)>= matriceRipples(i+1,1)
            matriceRipples(i+1,1)=matriceRipples(i,2)+1;
        end
    end
    
    for i=1:length(matriceRipples)
        on=round(matriceRipples(i,1));
        off=round(matriceRipples(i,2));
        if on==0Val
            on=1;
            
        end
    end
    
    LocR00D=x00(matriceRipples(:,1),2);
    LocR00F=x00(matriceRipples(:,2),2);
    
    matriceRipplesSec=zeros(length(locR00),3);
    matriceRipplesSec(:,1)=LocR00D;
    matriceRipplesSec(:,2)=LocR00F;
    matriceRipplesSec(:,3)=LocR00;
    
    assignin('base',sprintf('CH%d_HPCMatriceRipplesRSec',g(j)),matriceRipplesSec);
    assignin('base',sprintf('CH%d_HPCMatriceRipplesR',g(j)),matriceRipples);
    
    % Saving Ripples predetection
    cd(RipplesDataDirectory)
    
    if exist('HPC_Data_Ripples_loc.mat') == 0
        save('HPC_Data_Ripples_loc.mat', sprintf('CH%d_HPCRipplesR',g(j)),sprintf('CH%d_HPCRipplesAbsR',g(j)),...
            sprintf('CH%d_HPCLocRipples',g(j)),sprintf('CH%d_HPCAbsDownR',g(j)),sprintf('CH%d_HPCMatriceRipplesR',g(j)),...
            sprintf('CH%d_HPCMatriceRipplesRSec',g(j)),sprintf('CH%d_HPCHalfProminence',g(j)),...
            sprintf('CH%d_HPCThreshold',g(j)));
    else
        save('HPC_Data_Ripples_loc.mat', sprintf('CH%d_HPCRipplesR',g(j)),sprintf('CH%d_HPCRipplesAbsR',g(j)),...
            sprintf('CH%d_HPCLocRipples',g(j)),sprintf('CH%d_HPCAbsDownR',g(j)),...
            sprintf('CH%d_HPCMatriceRipplesR',g(j)),sprintf('CH%d_HPCMatriceRipplesRSec',g(j)),...
            sprintf('CH%d_HPCHalfProminence',g(j)),sprintf('CH%d_HPCThreshold',g(j)), '-append')
    end
    
    clearvars -except counter MatlabDirectory RipplesDataDirectory BaseDirectory BrutDirectory DataExtracted...
        LFPsDirectory DataBrutDirectory np g counter initialVars SR SRi rippleband...
        smoothF dSR i SWsign
    
end

load ( 'HPC_Data_Ripples_loc.mat');