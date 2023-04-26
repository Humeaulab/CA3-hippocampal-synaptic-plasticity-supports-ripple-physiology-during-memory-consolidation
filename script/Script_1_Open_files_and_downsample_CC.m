
BaseDirectory=''; % Working folder
BrutDirectory=''; % Directory for folder containing raw electrophysiological recordings from openephys.
DataExtracted=''; % Analysis directory

cd(BaseDirectory)
LFPsDirectory=uigetdir(DataExtracted,'select the LFPs Folder'); % Select .mat files containing the downsampled and referenced LFPs from the hippocampus and mPFC.
DataBrutDirectory = {''}; % Raw electrophysiological files from openephys.

brain_areas=["HPC","PFC"]; 
gPFC = [17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32]'; % PFC channel numbers.
g = [11,12,13,14]'; % HPC channel numbers.
gHPC = [11,12,13,14]'; % HPC channel numbers.
g_all=sort([gPFC;g]);

% Give the number of parts in the exp.
prompt = {'Enter the number of electrophysiology parts'}; % If multiple files have been recorded during a single recording session.
dlgtitle = 'Number of part';
definput = {'4'};
dims = [1 50];
answer = inputdlg(prompt,dlgtitle,dims,definput);
np = str2num(answer{1});


cd(BaseDirectory)
if exist('Info.mat') == 0 % Info.mat contains experiment details.
    save('Info.mat', 'BaseDirectory','BrutDirectory','DataExtracted','LFPsDirectory','DataBrutDirectory',...
        'gPFC','g','gHPC','g_all','brain_areas','np')
else
    save('Info.mat', 'BaseDirectory','BrutDirectory','DataExtracted','LFPsDirectory','DataBrutDirectory',...
        'gPFC','g','gHPC','g_all','brain_areas','np', '-append')
end

clearvars -except BaseDirectory BrutDirectory DataExtracted LFPsDirectory DataBrutDirectory selected_channels mPFC_chan dHPC_chan brain_areas

%% import of electrophysiology data from OpenEphys files

cd(BaseDirectory)
load('Info.mat')

for i = 1:size(g_all,1)
    Rawchannel_ds=[];
    for k=1:np
        cd(DataBrutDirectory{k})
        filename = fullfile(cd, sprintf('100_CH%d.continuous', g_all(i,1)));
        [data, timestamps, info] = load_open_ephys_data_faster(filename);
        Rawdata(:,1) = data/1000;
        Rawdata(:,2) = timestamps; 
        Rawdata_ds=downsample(Rawdata,30);
        if exist(sprintf('Record_start_%d',k))==0
            Record_start = Rawdata(1,2);
            assignin('base',sprintf('Record_start_%d',k),Record_start);
        end
        if exist(sprintf('Rawdata_info_%d',k))==0
            Rawdata_info = info;
            assignin('base',sprintf('Rawdata_info_%d',k),Rawdata_info);
        end
        Rawchannel_ds=[Rawchannel_ds;Rawdata_ds];
        
        clear data timestamps Rawdata Rawdata_ds filename
    end
    for l=1:size(brain_areas,2)
        region=eval(sprintf('g%s',brain_areas(l)));
        if ismember(g_all(i),region)==1
            assignin('base',sprintf('%s_%sds',info.header.channel,brain_areas(l)),Rawchannel_ds);
            cd(BaseDirectory)
            cd(LFPsDirectory)
            if exist(sprintf('%s_Data_nex.mat',brain_areas(l))) == 0
                save(sprintf('%s_Data_nex.mat', brain_areas(l)),sprintf('%s_%sds',info.header.channel,brain_areas(l)))
            else
                save(sprintf('%s_Data_nex.mat',brain_areas(l)), sprintf('%s_%sds',info.header.channel,brain_areas(l)),'-append')
            end
        end
    end
    
    clear Rawchannel_ds info
end

cd(BaseDirectory)
save('Info.mat','Record_start','Rawdata_info', '-append')

clearvars -except BaseDirectory

%% After check on NEX, last preprocessing step
cd(BaseDirectory)
load('Info.mat')

prompt = {'Enter the number of Dead Channels'};
dlgtitle = 'Number of channels';
definput = {'2'};
dims = [1 50];
answer = inputdlg(prompt,dlgtitle,dims,definput);
nDC = str2num(answer{:});

prompt = {'Enter the Dead Channel number :Carefull, this is the order of appearance of the channels in matlab : Enter space-separated numbers:'};
dlgtitle = 'Number :';
definput = {'0'};
dims = [1 50];
answer = inputdlg(prompt,dlgtitle,dims,definput);
D = str2num(answer{:});
DPosition=D(:);
DHPC=[];
DPFC=[];
for i=1:size(D,2)
   if ismember(g_all(D(i)),g)==1
       a=find(g==g_all(D(i)));
       DHPC=[DHPC;a];
   end
   if ismember(g_all(D(i)),gPFC)==1
       a=find(gPFC==g_all(D(i)));
       DPFC=[DPFC;a];
   end
end


%Returns channels that are not Dead
[C,CPosition]=setdiff(g,g(DHPC),'sorted');
[CPFC,CPositionPFC]=setdiff(gPFC,gPFC(DPFC),'sorted');


[C_all,CPosition_all]=setdiff(g_all,g_all(D),'sorted');

save('Info.mat','DPosition','CPosition','CPositionPFC','CPosition_all','nDC','-append')

%% Data referencing
deadch = [];
counter = 1;

for i = 1:size(brain_areas,2)
    cd(LFPsDirectory)
    load(sprintf('%s_Data_nex.mat',brain_areas(i)))
    region=eval(sprintf('g%s',brain_areas(i)));
    for l =1:length(CPosition_all)
        if find(deadch == g_all(CPosition_all(l)))
        else
            if ismember(g_all(CPosition_all(l)),region)==1
                tempo= eval(sprintf('CH%d_%sds',g_all(CPosition_all(l)),brain_areas(i)));
                tempref(:,counter) = tempo(:,1);
            
                counter = counter + 1;
                clear tempo
            end
        end    
    end
    ref = median(tempref,2); % median, can be mean. add option for substracting one channel only
    assignin('base', sprintf('%s_ref',brain_areas(i)), ref);
    
    cd(BaseDirectory)
    cd(LFPsDirectory)
    
    if exist(sprintf('%sRef_Data.mat',brain_areas(i))) == 0
        save(sprintf('%sRef_Data.mat',brain_areas(i)), sprintf('%s_ref',brain_areas(i)))
    else
        save(sprintf('%sRef_Data.mat',brain_areas(i)), sprintf('%s_ref',brain_areas(i)), '-append')
    end
    counter = 1;
    clearvars -except BaseDirectory BrutDirectory DataExtracted LFPsDirectory DataBrutDirectory nDC DPosition...
        gPFC g gHPC g_all CPosition_all deadch counter i brain_areas
end


clearvars -except BaseDirectory BrutDirectory DataExtracted LFPsDirectory DataBrutDirectory nDC DPosition g_all CPosition_all gPFC g gHPC brain_areas
counter = 1;

cd(BaseDirectory)

for i = 1:size(brain_areas,2)
    cd(LFPsDirectory)
    region=eval(sprintf('g%s',brain_areas(i)));
    load((sprintf('%s_Data_nex.mat',brain_areas(i))))
    load(sprintf('%sRef_Data.mat',brain_areas(i)), sprintf('%s_ref',brain_areas(i)))
    
    for k = 1:size(region,1)
        temp = eval(sprintf('CH%d_%sds',region(k),brain_areas(i)));
        tempref= eval(sprintf('%s_ref',brain_areas(i)));
        tempi = temp(:,1) - tempref;
        fs = 1000; % sampling frequency
        d = designfilt('bandstopiir', 'FilterOrder', 2, 'HalfPowerFrequency1', 49, 'HalfPowerFrequency2', 51,...
            'DesignMethod', 'butter', 'SampleRate', fs); % Notch filter 50Hz.
        tempif(:,1) = filtfilt(d,tempi(:,1));
        temp2=[tempif,temp(:,2)];
        assignin('base', sprintf('CH%d_%srds',region(k),brain_areas(i)), temp2);
        counter = counter + 1;
        if exist(sprintf('%s_Dataref_nex.mat',brain_areas(i))) == 0
            save(sprintf('%s_Dataref_nex.mat',brain_areas(i)), sprintf('CH%d_%srds',region(k),brain_areas(i)))
        else
            save(sprintf('%s_Dataref_nex.mat',brain_areas(i)), sprintf('CH%d_%srds',region(k),brain_areas(i)), '-append')
        end
        clear temp temp2 tempref tempi tempif area
    end
end

clearvars -except BaseDirectory BrutDirectory DataExtracted LFPsDirectory DataBrutDirectory nDC DPosition g_all gPFC g gHPC CPosition_all brain_areas
counter = 1;
cd(BaseDirectory)

