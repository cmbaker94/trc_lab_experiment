function Tinfo = trial_files(Tinfo)
% This function will grab information about the data based on the wave
% conditions
% INPUT: 
% Tinfo - wave conditions
% OUTPUT: 
% sz, is, lidarfile, Tcam - trial times/names for different instruments

Tinfo.datapath    = 'E:\';

if Tinfo.spread == 0
    if Tinfo.Hs == 0.3
        if Tinfo.Tp == 2
            % In situ
            Tinfo.sz.tdate          = '09-01-2018-2213UTC';
            Tinfo.is.tdate          = '09-06-2018-1559UTC';
            % Lidar (LI)
            Tinfo.lidar.tdate       = '2018-09-01-22-13-48';
            % Stereo Reconstructions (SR)
            Tinfo.cam.tstart        = '09-01-2018-2214UTC';      	% time starting collection based on spreadsheet
            Tinfo.cam.tdate         = '09-01-2018-2155UTC';         % trial date and time - format ex: 09-01-2018-2155UTC
            Tinfo.offsets           = [-485;-69]; %-487 % index offset relative to camera [in situ, lidar]
        end
    elseif Tinfo.Hs == 0.35
        Tinfo.sz.tdate            = '';
        Tinfo.is.tdate            = '09-08-2018-1742UTC';
        % Lidar (LI)
        Tinfo.lidar.tdate        = '2018-09-08-17-42-36';
        % Stereo Reconstructions (SR)
        Tinfo.cam.tstart        = '09-08-2018-1742UTC';         % time starting collection based on spreadsheet
        Tinfo.cam.tdate         = '09-08-2018-1730UTC';         % trial date and time - format ex: 09-01-2018-2155UTC
        Tinfo.offsets           = [0;0]; % index offset relative to camera [in situ, lidar]
    elseif Tinfo.Hs == 0.25
        Tinfo.sz.tdate            = '09-01-2018-2001UTC';
        Tinfo.is.tdate            = '09-06-2018-1952UTC';
        % Lidar (LI)
        Tinfo.lidar.tdate        = '2018-09-01-20-02-03';
        % Stereo Reconstructions (SR)
        Tinfo.cam.tstart        = '09-01-2018-2002UTC';         % time starting collection based on spreadsheet
        Tinfo.cam.tdate         = '09-01-2018-1950UTC';         % trial date and time - format ex: 09-01-2018-2155UTC
        Tinfo.offsets           = [0;0]; % index offset relative to camera [in situ, lidar]
    else
        Tinfo.sz.tdate            = [];
        Tinfo.is.tdate            = [];
        % Lidar (LI)
        Tinfo.lidar.tdate        = [];
        % Stereo Reconstructions (SR)
        Tinfo.cam.tstart        = [];         % time starting collection based on spreadsheet
        Tinfo.cam.tdate         = [];         % trial date and time - format ex: 09-01-2018-2155UTC
        Tinfo.offsets           = [0;0]; % index offset relative to camera [in situ, lidar]
    end
elseif Tinfo.spread == 10
    if Tinfo.Hs == 0.25
        if Tinfo.Tp == 2
            if Tinfo.tide == 1.07
                % In situ
                Tinfo.sz.tdate          = '08-30-2018-2026UTC';
                Tinfo.is.tdate          = '09-06-2018-2342UTC';
                % Lidar (LI)
                Tinfo.lidar.tdate       = '2018-09-06-23-42-19';
                % Stereo Reconstructions (SR)
                Tinfo.cam.tstart        = '09-06-2018-2342UTC';     	% time starting collection based on spreadsheet
                Tinfo.cam.tdate         = '09-06-2018-2334UTC';         % trial date and time - format ex: 09-01-2018-2155UTC
                Tinfo.offsets           = [0;0]; % index offset relative to camera [in situ, lidar]
            end
        end
    end
elseif Tinfo.spread == 20
    if Tinfo.Hs == 0.3
        if Tinfo.Tp == 2
             % In situ
            Tinfo.sz.tdate          = '08-30-2018-2222UTC';
            Tinfo.is.tdate          = '09-06-2018-1655UTC';
            % Lidar (LI)
            Tinfo.lidar.tdate       = '2018-09-06-16-55-27';
            % Stereo Reconstructions (SR)
            Tinfo.cam.tstart        = '09-06-2018-1655UTC';     	% time starting collection based on spreadsheet
            Tinfo.cam.tdate         = '09-06-2018-1646UTC';         % trial date and time - format ex: 09-01-2018-2155UTC
            Tinfo.offsets           = [0;0]; % index offset relative to camera [in situ, lidar]
%             display('DROPPIING!!!!')
%             % In situ
%             Tinfo.sz.tdate          = '08-30-2018-2222UTC';
%             Tinfo.is.tdate          = '09-06-2018-1655UTC';
%             % Lidar (LI)
%             Tinfo.lidar.tdate       = '2018-08-30-22-22-39';
%             % Stereo Reconstructions (SR)
%             Tinfo.cam.tstart        = '08-30-2018-2222UTC';     	% time starting collection based on spreadsheet
%             Tinfo.cam.tdate         = '08-30-2018-2216UTC';         % trial date and time - format ex: 09-01-2018-2155UTC
%             Tinfo.offsets           = [0;0]; % index offset relative to camera [in situ, lidar]
        elseif Tinfo.Tp == 3
            % In situ
            Tinfo.sz.tdate          = '09-01-2018-1626UTC';%'08-30-2018-2222UTC';
            Tinfo.is.tdate          = '09-06-2018-1655UTC';
            % Lidar (LI)
            Tinfo.lidar.tdate       = '2018-09-01-16-26-49';%'2018-08-30-22-22-39';
            % Stereo Reconstructions (SR)
            Tinfo.cam.tstart        = '09-01-2018-1627UTC';     	% time starting collection based on spreadsheet
            Tinfo.cam.tdate         = '09-01-2018-1530UTC';         % trial date and time - format ex: 09-01-2018-2155UTC
            Tinfo.offsets           = [0;0]; % index offset relative to camera [in situ, lidar]
        end
    elseif Tinfo.Hs == 0.25
        if Tinfo.Tp == 2
            if Tinfo.tide == 1.07
                % In situ
                Tinfo.sz.tdate          = '08-30-2018-1904UTC';
                Tinfo.is.tdate          = '09-06-2018-2048UTC';
                % Lidar (LI)
                Tinfo.lidar.tdate       = '2018-09-06-20-48-56';
                % Stereo Reconstructions (SR)
                Tinfo.cam.tstart        = '09-06-2018-2049UTC';
                Tinfo.cam.tdate         = '09-06-2018-2043UTC';
                Tinfo.offsets           = [0;0]; % index offset relative to camera [in situ, lidar]
            end
        end
    else
        Tinfo.sz.tdate            = [];
        Tinfo.is.tdate            = [];
        % Lidar (LI)
        Tinfo.lidar.tdate        = [];
        % Stereo Reconstructions (SR)
        Tinfo.cam.tstart        = [];         % time starting collection based on spreadsheet
        Tinfo.cam.tdate         = [];         % trial date and time - format ex: 09-01-2018-2155UTC
        Tinfo.offsets           = [0;0]; % index offset relative to camera [in situ, lidar]
    end
elseif Tinfo.spread == 30
    if Tinfo.Hs == 0.3
        % In situ
        Tinfo.sz.tdate          = '08-29-2018-2255UTC';         % 2119+20 min
        Tinfo.is.tdate          = '09-06-2018-1748UTC';
        % Lidar (LI)
        Tinfo.lidar.tdate       = '2018-08-29-22-55-28';
        % Stereo Reconstructions (SR)
        Tinfo.cam.tstart        = '08-29-2018-2255UTC';         % time starting collection based on spreadsheet
        Tinfo.cam.tdate         = '08-29-2018-2236UTC';         % trial date and time - format ex: 09-01-2018-2155UTC
        Tinfo.offsets           = [3655; 341];
        Tinfo.offsets           = [0;0]; % index offset relative to camera [in situ, lidar]
    elseif Tinfo.Hs == 0.25
        if Tinfo.Tp == 2
            if Tinfo.tide == 1.07
                % In situ
                Tinfo.sz.tdate          = '09-01-2018-2101UTC';
                Tinfo.is.tdate          = '09-07-2018-0035UTC';
                % Lidar (LI)
                Tinfo.lidar.tdate       = '2018-09-01-21-01-57';
                % Stereo Reconstructions (SR)
                Tinfo.cam.tstart        = '09-01-2018-2102UTC';
                Tinfo.cam.tdate         = '09-01-2018-2052UTC';
                Tinfo.offsets           = [-182; -36]; % index offset relative to camera [in situ, lidar]
            end
        end
    else
        Tinfo.sz.tdate            = [];
        Tinfo.is.tdate            = [];
        % Lidar (LI)
        Tinfo.lidar.tdate        = [];
        % Stereo Reconstructions (SR)
        Tinfo.cam.tstart        = [];         % time starting collection based on spreadsheet
        Tinfo.cam.tdate         = [];         % trial date and time - format ex: 09-01-2018-2155UTC
        Tinfo.offsets           = [0;0]; % index offset relative to camera [in situ, lidar]
    end
elseif Tinfo.spread == 40
    if Tinfo.Hs == 0.3
        if Tinfo.tide == 1.07
            if Tinfo.Tp == 2
                % In situ
                Tinfo.sz.tdate            = '08-29-2018-2358UTC'; % 2119+20 min
                Tinfo.is.tdate            = '09-06-2018-1841UTC';
                % Lidar (LI)
                Tinfo.lidar.tdate = '2018-08-29-23-58-16';
                % Stereo Reconstructions (SR)
                Tinfo.cam.tstart        = '08-29-2018-2359UTC';                   % time starting collection based on spreadsheet
                Tinfo.cam.tdate         = '08-29-2018-2349UTC';      % trial date and time - format ex: 09-01-2018-2155UTC
                Tinfo.offsets           = [-(372/8)*100+417,-(371/8)*10]; % [-(372/8)*100+417,-(372/8)*10]; index offset relative to camera [in situ, lidar]
%                             display('This case is dropping frames!!!!')
%                             % In situ
%                             Tinfo.sz.tdate            = '08-30-2018-2129UTC'; % 2119+20 min
%                             Tinfo.is.tdate            = '09-06-2018-1841UTC';
%                             % Lidar (LI)
%                             Tinfo.lidar.tdate = '2018-08-30-21-29-26_Velodyne-HDL-32-Data_gridded';
%                             % Stereo Reconstructions (SR)
%                             Tinfo.cam.tstart        = '08-30-2018-2129UTC';                   % time starting collection based on spreadsheet
%                             Tinfo.cam.tdate         = '08-30-2018-2119UTC';      % trial date and time - format ex: 09-01-2018-2155UTC
%                             Tinfo.offsets           = [0,0]; %2925 index offset relative to camera [in situ, lidar]
            elseif Tinfo.Tp == 3
                % In situ
                Tinfo.sz.tdate          = '09-01-2018-1757UTC';
                Tinfo.is.tdate          = '09-07-2018-1737UTC';
                % Lidar (LI)
                Tinfo.lidar.tdate       = '2018-09-01-17-57-52';
                % Stereo Reconstructions (SR)
                Tinfo.cam.tstart        = '09-01-2018-1758UTC';      	% time starting collection based on spreadsheet
                Tinfo.cam.tdate         = '09-01-2018-1740UTC';         % trial date and time - format ex: 09-01-2018-2155UTC
                Tinfo.offsets           = [0;0]; %-487 % index offset relative to camera [in situ, lidar]
            end
        elseif Tinfo.tide == 1.00
            % In situ
            Tinfo.sz.tdate            = '08-31-2018-2232UTC';%'08-30-2018-2129UTC'; % 2119+20 min
            Tinfo.is.tdate            = '09-06-2018-1841UTC';
            % Lidar (LI)
            Tinfo.lidar.tdate = '2018-08-31-22-32-25';%'2018-08-30-21-29-26_Velodyne-HDL-32-Data_gridded';
            % Stereo Reconstructions (SR)
            Tinfo.cam.tstart        = '08-31-2018-2232UTC';                   % time starting collection based on spreadsheet
            Tinfo.cam.tdate         = '08-31-2018-2225UTC';      % trial date and time - format ex: 09-01-2018-2155UTC
            Tinfo.offsets           = [3130,293]; %2925 index offset relative to camera [in situ, lidar]
        end
    elseif Tinfo.Hs == 0.25
         % In situ
        Tinfo.sz.tdate            = '08-30-2018-1655UTC'; % 2119+20 min
        Tinfo.is.tdate            = '09-06-2018-2150UTC'; %day 9, trila 7, 2018-09-06-21-50
        % Lidar (LI)
        Tinfo.lidar.tdate         = '2018-08-30-16-55-17';
        % Stereo Reconstructions (SR)
        Tinfo.cam.tstart        = '08-30-2018-1655UTC';                   % time starting collection based on spreadsheet
        Tinfo.cam.tdate         = '08-30-2018-1634UTC';      % trial date and time - format ex: 09-01-2018-2155UTC
        Tinfo.offsets           = [((0-2.17)/8)*100; 0]; % difference between lidar and cam unknown, index offset relative to camera [in situ, lidar]
    elseif Tinfo.Hs == 0.2
        Tinfo.sz.tdate          = '08-30-2018-1534UTC';         % 2119+20 min
        Tinfo.is.tdate          = '09-06-2018-2248UTC';
        % Lidar (LI)
        Tinfo.lidar.tdate       = '2018-08-30-15-34-16';
        % Stereo Reconstructions (SR)
        Tinfo.cam.tstart        = '08-30-2018-1534UTC';                   % time starting collection based on spreadsheet
        Tinfo.cam.tdate         = '08-30-2018-1518UTC';      % trial date and time - format ex: 09-01-2018-2155UTC
        Tinfo.offsets           = [0;0]; % index offset relative to camera [in situ, lidar]
    else
        Tinfo.sz.tdate            = [];
        Tinfo.is.tdate            = [];
        % Lidar (LI)
        Tinfo.lidar.tdate        = [];
        % Stereo Reconstructions (SR)
        Tinfo.cam.tstart        = [];         % time starting collection based on spreadsheet
        Tinfo.cam.tdate         = [];         % trial date and time - format ex: 09-01-2018-2155UTC
        Tinfo.offsets           = [0;0]; % index offset relative to camera [in situ, lidar]
    end
end