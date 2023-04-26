clear all
warning off

BaseDirectory='\\filer2-iins\HUMEAU\EVERYONE\DATA IN VIVO\CC-202302-2994\E_Day 3';

cd(BaseDirectory);
load('Info.mat')

prompt = {'Enter the number of Ripple Negative Channels'};
dlgtitle = 'Number of channels';
definput = {'2'};
dims = [1 50];
answer = inputdlg(prompt,dlgtitle,dims,definput);
nNC = str2num(answer{:});

prompt = {'Enter the Ripples Negative Channel number :Carefull, this is the order of appearance of the channels in matlab : Enter space-separated numbers:'};
dlgtitle = 'Number :';
definput = {'0'};
dims = [1 50];
answer = inputdlg(prompt,dlgtitle,dims,definput);
N = str2num(answer{:});
NPosition=N(:);


prompt = {'Enter the number of Ripple Positive Channels'};
dlgtitle = 'Number of channels';
definput = {'4'};
dims = [1 50];
answer = inputdlg(prompt,dlgtitle,dims,definput);
nPC = str2num(answer{:});

prompt = {'Enter the Ripples Positive Channel number :Carefull, this is the order of appearance of the channels in matlab : Enter space-separated numbers:'};
dlgtitle = 'Number :';
definput = {'0'};
dims = [1 50];
answer = inputdlg(prompt,dlgtitle,dims,definput);
P = str2num(answer{:});
PPosition=P(:);


cd(BaseDirectory)
if exist('Info.mat') == 0
    save('Info.mat', 'nNC','NPosition','nPC','PPosition')
else
    save('Info.mat', 'nNC','NPosition','nPC','PPosition', '-append')
end

%% Noise intervals detection

cd(LFPsDirectory)
noise_threshold=2;
noise_intervals=[];
for i=1:size(DPosition,1)
    temp=[];
    if ismember(g_all(DPosition(i)),g)==1
        load('HPC_Data_nex',sprintf('CH%d_HPCds',g_all(DPosition(i))))
        temp=eval(sprintf('CH%d_HPCds',g_all(DPosition(i))));
    end
    if ismember(g_all(DPosition(i)),gPFC)==1
        load('PFC_Data_nex',sprintf('CH%d_PFCds',g_all(DPosition(i))))
        temp=eval(sprintf('CH%d_PFCds',g_all(DPosition(i))));
    end
    if ismember(g_all(DPosition(i)),gPPC)==1
        load('PPC_Data_nex',sprintf('CH%d_PPCds',g_all(DPosition(i))))
        temp=eval(sprintf('CH%d_PPCds',g_all(DPosition(i))));
    end
    if isempty(temp)==0
        a=find(temp(:,1)>=noise_threshold|temp(:,1)<=-noise_threshold);
        matrix=[temp(a(1),2),0];
        for l=2:size(a,1)
            if a(l)>a(l-1)+1000
                matrix(end,2)=temp(a(l-1),2);
                matrix(end+1,:)=temp(a(l),2);
            end
        end
        f=find(matrix(:,1)==matrix(:,2));
        matrix(f,2)=matrix(f,2)+0.0001;
        noise_intervals=[noise_intervals;matrix];
    end
end

noise_intervals=sortrows(noise_intervals);

cd(RipplesDataDirectory)
save('HPC_MatricesRipples','noise_intervals','noise_threshold');

clearvars -except BaseDirectory BrutDirectory DataExtracted LFPsDirectory DataBrutDirectory np nDC DPosition...
        g CPosition gPFC CPositionPFC g_all nNC NPosition nPC PPosition RipplesDataDirectory

%% False positive ripples elimination

load ('HPC_MatricesRipples.mat')

for j = 1:length (PPosition) %Positive channels
    load('HPC_Data_Ripples_loc.mat', sprintf('CH%d_HPCLocRipples',g(PPosition(j))));
    loctemp = eval(sprintf('CH%d_HPCLocRipples',g(PPosition(j))));
    load('HPC_Data_Ripples_loc.mat',sprintf('CH%d_HPCHalfProminence',g(PPosition(j))));
    HPTemp=eval(sprintf('CH%d_HPCHalfProminence',g(PPosition(j))));
    
    temp = noise_intervals;
        
    for n= 1:length (loctemp)
        for m = 1:length(temp)
            a=find(temp(m,1)<=loctemp(n) && temp(m,2)>=loctemp(n)) ;
            if isempty(a)==1
            else stocka(n,1)=loctemp(n);
                stockHP(n,j)=HPTemp (n);
            end
        end
    end    
    
    Zstocka =  unique (stocka,'sorted');
    Zzstocka = nonzeros (stocka);
   
    for k= 1:length (Zzstocka)
        for l=1:length (loctemp)
            if Zzstocka (k) == loctemp(l)
                loctemp(l) = 0;
            end
        end
    end
    
    newloctemp = nonzeros (loctemp);
    
    ZzstockHP = nonzeros (stockHP(:,j));
    for k= 1:length (ZzstockHP)
        for l=1:length (HPTemp)
            if ZzstockHP (k) == HPTemp(l)
                HPTemp(l) = 0;
            end
        end
    end
    
    newHPTemp = nonzeros (HPTemp);
    
    assignin('base',sprintf('CH%d_HPCLoc',g(PPosition(j))),newloctemp);
    assignin('base',sprintf('CH%d_HPCHP',g(PPosition(j))),newHPTemp);
    if exist('HPC_MatricesRipples.mat') == 0;
        save('HPC_MatricesRipples.mat',sprintf('CH%d_HPCLoc',g(PPosition(j))),sprintf('CH%d_HPCHP',g(PPosition(j))));
    else
        save('HPC_MatricesRipples.mat', sprintf('CH%d_HPCLoc',g(PPosition(j))),sprintf('CH%d_HPCHP',g(PPosition(j))),'-append')
    end
end

clearvars -except BaseDirectory BrutDirectory DataExtracted LFPsDirectory DataBrutDirectory np nDC DPosition...
        g CPosition gPFC CPositionPFC g_all nNC NPosition nPC PPosition

