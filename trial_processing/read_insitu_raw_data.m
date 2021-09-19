function [insituname] = read_insitu_raw_data(day,trial)
% read, reformat, and save in situ data
% this code will use the processed data from the HWRL insitu instruments
% and ideally save it in a 'more' user friendly format

%% STEP 1: Create paths

% general path in computer to raw and processed data
genpath = 'E:\data\';
% create path for raw data
rawpath = ['G:\PRJ-1873\inter\Random',day,'\Trial',trial,'\'];

%% STEP 2: Data Types

data_type       = {'press','u','v','w','wg'};
% inst_no.press   = {'1','2','4','5','6','7','8','9'};
% inst_no.u       = {'1','3','4','5','7','8'};
% inst_no.v       = {'1','2','3','4','6','8','9','10','11','12'};
% inst_no.w       = {'1','2','4','5','9','10','12'};


%% STEP 3: Record General Information

% First processes the pressure and wavegage data
data_type = {'wg';'press';'u';'v';'w'};
for i = 1:length(data_type)
    datafiles           = dir([rawpath,data_type{i},'*']);
    
    for j = 1:length(datafiles(:))
        file        = datafiles(j).name(1:end-4);
        instnum     = file(regexp(file,'\d'));
        
        if strcmp(data_type{i},'w') && isempty(str2num(file(2:end)))
            continue
        end
        
        % Load text file
        fileID      = fopen([rawpath,file,'.txt'],'r');
        filetemp    = textscan(fileID,'%s','Delimiter','\n');
        ifind       = cell2table(filetemp{1,1});
        for k = 1:length(ifind{:,1})
            linetemp = str2double(ifind{k,1});
            if ~isnan(linetemp)
                break
            end
        end
        
        % grab introduction information
        intro = filetemp{1,1}(1:k-1,1);
        
        % find location of instruments
        [xyzd,stillwat] = find_insitu_xyzd(intro);
        
        % pick data
        if contains(file,'press') || contains(file,'wg')
            eval([data_type{i},'.',file,' = str2double(filetemp{1,1}(k:end,1));'])
            eval([data_type{i},'.xyzd.',file,' = xyzd;'])
        elseif contains(file,'u') 
            eval(['vel.',file,' = str2double(filetemp{1,1}(k:end,1));'])
            eval(['vel.xyzd.adv',instnum,' = xyzd;'])
        elseif contains(file,'v') || contains(file,'w')
            eval(['vel.',file,' = str2double(filetemp{1,1}(k:end,1));'])
        end  
        clear *temp x y z xyz
    end
    
    % General information
    % unit
    untxttemp   = find(contains(intro,'% CalibrationUnits'));
    unittemp    = intro{untxttemp}(22:end);
    
    if strcmp(data_type{i},'wg') || strcmp(data_type{i},'press')
        eval([data_type{i},'.unit   = unittemp;'])
        % rate
        Hztxttemp   = find(contains(intro,'% SampleRate'));
        Hztemp      = str2num(intro{Hztxttemp}(15:end-3)); 
        eval([data_type{i},'.Hz     = Hztemp;'])
    elseif strcmp(data_type{i},'u') || strcmp(data_type{i},'v') || strcmp(data_type{i},'w')
        eval(['vel.unit             = unittemp;'])
    end

    clear temp
end


%% STEP 4: Copy and store general trial information

% time in date format
Tdatetemp = find(contains(intro,'% CalibratedStartDateTimeUTC'));
tdate = intro{Tdatetemp};
tdate = tdate(36:end);

% time in num format
Tnumtemp = find(contains(intro,'% CalibratedStartDateNumUTC'));
tnum = intro{Tnumtemp};
tnum = str2double(tnum(31:end));

% trial conditions
condtemp = find(contains(intro,'% TrialConditions'));
cond = intro{condtemp};
cond = cond(39:end);

clear *temp Intro

Tinfo.time = tdate;
Tinfo.random = day;
Tinfo.trial = trial;
Tinfo.tnum = tnum;
Tinfo.cond = cond;
Tinfo.stillwat = stillwat;

%% STEP 5: Convert Pressure to watlev

% % find the period
% id      = strfind(cond,'T=');
% T       = str2num(cond(id+2));
% 
% % depth of instruments
% % NOTE STILL NEED TO FIND DEPTH AT THE XYZ LOCATIONS OF PRESSURE GAGES
% % This can be based on the x-location of the instrument d.
% % JK FOR NOW...
% % Just giving an offset of the instrument below
% offsetinst = 0.03;
% 
% [eta] = convert_press2eta(press,T,offsetinst);

%% STEP 6: Save Information 

clear random trial tnum cond watlev

display('Make sure all trial details are the same')

% % create path for processed data
procpath = [genpath,'processed\insitu\'];%random',day,'\Trial',trial,'\'];
% eval(['!mkdir ',genpath,'processed\insitu\Trials\random',day])
trialfolder = [tdate(6:7),'-',tdate(9:10),'-',tdate(1:4),'-',tdate(12:13),tdate(15:16),'UTC'];
eval(['!mkdir ',procpath,trialfolder])

sname = [tdate(6:7),'-',tdate(9:10),'-',tdate(1:4),'-',tdate(12:13),tdate(15:16),'UTC-insitu'];
eval(['save -v7.3 ',procpath,trialfolder,'\',sname,' vel',' press',' wg',' Tinfo'])%' eta'

insituname = trialfolder;
end

