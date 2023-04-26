clear all
warning off

BaseDirectory='\\filer2-iins\HUMEAU\EVERYONE\DATA IN VIVO\CC-202302-2994\E_Day 3';

cd(BaseDirectory);
load('Info.mat')
i = 1;

RipplesAnalysisDirectory = uigetdir(BaseDirectory,'select the Analysis Ripples Folder');

cd(BaseDirectory)
if exist('Info.mat') == 0
    save('Info.mat', 'RipplesAnalysisDirectory')
else
    save('Info.mat', 'RipplesAnalysisDirectory', '-append')
end

SR=1000;
SRi=1000;

rippleband=[100 250];
smoothF=round(SR/(rippleband(1)));

dSR=50;

SWsign=0;
counter = 1;

for k=1:length (PPosition)
    cd(RipplesDataDirectory);
    load('HPC_Data_Ripples_loc.mat',sprintf('CH%d_HPCMatriceRipplesRSec',g(PPosition(k))));
    load('HPC_MatricesRipples.mat',sprintf('CH%d_HPCLoc',g(PPosition(k))));
    
    Matrice = eval(sprintf('CH%d_HPCMatriceRipplesRSec',g(PPosition(k))));
    NewMatrice = eval(sprintf('CH%d_HPCMatriceRipplesRSec',g(PPosition(k))));
    SelectedLoc= eval(sprintf('CH%d_HPCLoc',g((PPosition(k)))));
    
    for m = 1:length(Matrice)
        Lia = ismember(Matrice(m,3),SelectedLoc);
        if Lia == 1
        else
            NewMatrice(m,1)=0;
            NewMatrice(m,2)=0;
            NewMatrice(m,3)=0;
        end
        clearvars Lia
    end
    
    [ii,jj]=find(~NewMatrice);
    P=ii;
    RepFauxP = unique(P(:,1));
    
    load('HPC_Data_Ripples_loc.mat',sprintf('CH%d_HPCMatriceRipplesR',g(PPosition(k))));
    MatriceLoc = eval(sprintf('CH%d_HPCMatriceRipplesR',g(PPosition(k))));
    for p = 1:length(RepFauxP)
        q= RepFauxP(p);
        MatriceLoc (q,:) = 0;
    end
    
    NewMatriceIndex (:,1) = nonzeros (MatriceLoc(:,1));
    NewMatriceIndex (:,2) = nonzeros (MatriceLoc(:,2));
    NewMatriceIndex (:,3) = nonzeros (MatriceLoc(:,3));
    
    NewMatriceb (:,1) = nonzeros (NewMatrice(:,1));
    NewMatriceb (:,2) = nonzeros (NewMatrice(:,2));
    NewMatriceb (:,3) = nonzeros (NewMatrice(:,3));
    
    assignin('base',sprintf('CH%d_HPCSelectedMatriceRipplesIndex',g(PPosition(k))),NewMatriceIndex);
    assignin('base',sprintf('CH%d_HPCSelectedMatriceRipples',g(PPosition(k))),NewMatriceb);
    
    
    if exist('HPC_Data_Ripplesrb.mat') == 0
        save('HPC_Data_Ripplesrb.mat', sprintf('CH%d_HPCSelectedMatriceRipples',g(PPosition(k))), sprintf('CH%d_HPCSelectedMatriceRipplesIndex',g(PPosition(k))));
    else
        save('HPC_Data_Ripplesrb.mat',sprintf('CH%d_HPCSelectedMatriceRipples',g(PPosition(k))), sprintf('CH%d_HPCSelectedMatriceRipplesIndex',g(PPosition(k))), '-append')
    end
    
    clearvars Matrice NewMatrice Selectedloc NewMatriceb NewMatriceIndex MatriceLoc RepFauxP [ii,jj] P
end

for j = 1:length (PPosition)
    cd(BaseDirectory);
    cd(LFPsDirectory);
    load('HPC_Dataref_nex.mat');
    
    x00 = eval(sprintf('CH%d_HPCrds',g(PPosition(j))));
    timeStampcolumn=x00(:,2);
    t00=(0:1/SR:length(x00)/SR-1/SR)'; 
    
    xD00 = BandPass((x00(:,1)),SR,rippleband(1),rippleband(2),3); 
    xD00(:,2)=timeStampcolumn;
    assignin('base',sprintf('CH%d_HPCrbRipples',g(PPosition(j))),xD00);
    xD00_inverted(:,1)=-(xD00(:,1));
    xD00_inverted(:,2)=xD00(:,2);
    assignin('base',sprintf('CH%d_HPCrbRipples_inverted',g(PPosition(j))),xD00_inverted);
    
    xrms00=smooth(abs(xD00(:,1)),smoothF,'moving'); 
    xrms00(:,2)=timeStampcolumn;
    assignin('base',sprintf('CH%d_HPCrbRipplesAbs',g(PPosition(j))),xrms00);
    
    
    cd(BaseDirectory);
    cd(RipplesDataDirectory);
    load (sprintf('HPC_Data_Ripplesrb.mat', sprintf('CH%d_HPCSelectedMatriceRipplesIndex',g(PPosition(j)))))
    matriceRipples = eval(sprintf('CH%d_HPCSelectedMatriceRipplesIndex',g(PPosition(j))));
   
    for i=1:size(matriceRipples)
        IntTemp=xD00(matriceRipples(i,1):matriceRipples(i,2),1);
        IntTemp_inverted=xD00_inverted(matriceRipples(i,1):matriceRipples(i,2),1);
        pks = findpeaks(IntTemp);
        pks_inverted = findpeaks(IntTemp_inverted);
        matriceRipples(i,11)=length(pks)+length(pks_inverted);
        temp=(matriceRipples(i,2)-matriceRipples(i,1))/SR;
        matriceRipples(i,10)=matriceRipples(i,11)/(2*temp);
        clear IntTemp IntTemp_inverted pks pks_inverted temp
        matriceRipples(i,4)=max(xrms00(round(matriceRipples(i,1)):round(matriceRipples(i,2)),1));
        matrice_amplitude = zeros(length(matriceRipples),4);
        [matrice_amplitude(i,1),matrice_amplitude(i,2)]=max(xD00(round(matriceRipples(i,1)):round(matriceRipples(i,2)),1));
        LocFiltree(i,1)= ((matriceRipples(i,1))+(matrice_amplitude(i,2))-1);
        [matrice_amplitude(i,3),matrice_amplitude(i,4)]=max(xD00_inverted(round(matriceRipples(i,1)):round(matriceRipples(i,2)),1));
        LocFiltree(i,3)= ((matriceRipples(i,1))+(matrice_amplitude(i,4))-1);
        matriceRipples(i,5)=matrice_amplitude(i,1)+(abs(matrice_amplitude(i,3)));
        matriceRipples(i,6)=sum(xrms00(round(matriceRipples(i,1)):round(matriceRipples(i,2)),1));
    end
    
    
    LocFiltree(:,2)=xD00(LocFiltree(:,1),2);
    LocFiltree(:,4)=xD00(LocFiltree(:,3),2);
    assignin('base',sprintf('CH%d_HPCrbLocFR',g(PPosition(j))),LocFiltree);
    assignin('base',sprintf('CH%d_HPCrbLocF',g(PPosition(j))),LocFiltree(:,2));
    assignin('base',sprintf('CH%d_HPCrbLocF_inverted',g(PPosition(j))),LocFiltree(:,4));
    
    matriceRipples(:,12)=LocFiltree(:,2);
    
    load('HPC_MatricesRipples.mat', sprintf('CH%d_HPCHP',g(PPosition(j))));
    tempou=eval(sprintf('CH%d_HPCHP',g(PPosition(j))));
    matriceRipples(:,13)= tempou;
    clearvars tempou
    
    load('HPC_Data_Ripplesrb.mat', sprintf('CH%d_HPCSelectedMatriceRipples',g(PPosition(j))));
    temp = eval(sprintf('CH%d_HPCSelectedMatriceRipples',g(PPosition(j))));
    matriceRipples(:,1)= temp(:,1);
    matriceRipples(:,2)= temp(:,2);
    matriceRipples(:,3)= temp(:,3);
    
    matriceRipples(:,7)=matriceRipples(:,2)-matriceRipples(:,1);
    matriceRipples(:,8)=matriceRipples(:,3)-matriceRipples(:,1);
    matriceRipples(:,9)=matriceRipples(:,2)-matriceRipples(:,3);
    
    RipplesInterval=matriceRipples(:,1:2);
    assignin('base',sprintf('CH%d_HPCrbRipplesInterval',g(PPosition(j))),RipplesInterval);
    
    interval=matriceRipples(:,1:2);
    assignin('base',sprintf('CH%d_HPCrbRipplesIntervals',g(PPosition(j))),interval);
    assignin('base',sprintf('CH%d_HPCrbMatriceRipples',g(PPosition(j))),matriceRipples);
    
    if exist('HPC_Data_Ripplesrb.mat') == 0
        save('HPC_Data_Ripplesrb.mat', sprintf('CH%d_HPCrbRipples',g(PPosition(j))),sprintf('CH%d_HPCrbRipplesAbs',g(PPosition(j))),...
            sprintf('CH%d_HPCrbMatriceRipples',g(PPosition(j))),sprintf('CH%d_HPCrbRipplesIntervals',g(PPosition(j))),...
            sprintf('CH%d_HPCrbLocF',g(PPosition(j))),sprintf('CH%d_HPCrbLocF_inverted',g(PPosition(j))),sprintf('CH%d_HPCrbLocFR',g(PPosition(j))))
    else
        save('HPC_Data_Ripplesrb.mat', sprintf('CH%d_HPCrbRipples',g(PPosition(j))),sprintf('CH%d_HPCrbRipplesAbs',g(PPosition(j))),...
            sprintf('CH%d_HPCrbMatriceRipples',g(PPosition(j))),sprintf('CH%d_HPCrbRipplesIntervals',g(PPosition(j))),...
            sprintf('CH%d_HPCrbLocF',g(PPosition(j))),sprintf('CH%d_HPCrbLocF_inverted',g(PPosition(j))),sprintf('CH%d_HPCrbLocFR',g(PPosition(j))),'-append')
    end
    
    clearvars -except RipplesDataDirectory BaseDirectory BrutDirectory DataExtracted LFPsDirectory...
        DataBrutDirectory nDC DPosition g CPosition nNC NPosition nPC PPosition np counter SR SRi rippleband smoothF dSR i SWsign
    
end


load('HPC_Data_Ripplesrb.mat');
for j = 1:length (PPosition)
    temp{j,1} = eval(sprintf('CH%d_HPCrbLocF',g(PPosition(j))));
end

temp2 = vertcat(temp{:});
HPC_AllRipplesF = temp2;
clearvars temp temp2

if exist('HPC_Data_Ripplesrb.mat') == 0
    save('HPC_Data_Ripplesrb.mat','HPC_AllRipplesF');
else
    save('HPC_Data_Ripplesrb.mat','HPC_AllRipplesF','-append');
end

