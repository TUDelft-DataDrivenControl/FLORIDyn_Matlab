function T_phi = CLC_read_probe(path_2_file, average_time, default_phi)
% Function to get closed loop behaviour also with the wind direction
%   Meant to be called with a probe sampling at the location of the
%   turbines. 
% INPUTS
%   path_2_file     Path to probes file
%   average_time    Time window over which the function should average
%   default_phi     If not enough / no past values are written, what should 
%                   be returned.


probes.nProbes        = 0;
probes.probesLocation = [];
probes.t              = [];
probes.u              = [];
probes.v              = [];
probes.w              = [];
it                    = 1;
% Read file
fID = fopen(path_2_file);

tline = fgetl(fID);
while ischar(tline)
    % Read probe locations
    if startsWith(tline,'# Probe')
        probes.nProbes = probes.nProbes + 1;
        iNum = find(tline=='(',1);
        probes.probesLocation(probes.nProbes,:) = ...
            str2num(tline(iNum+1:end-1));
    end

    % Read probe data
    if ~startsWith(tline,'#')
        iNum = find(tline~=' ',1);
        if isempty(iNum); continue; end 
        probes.t(it) = str2double(tline(iNum:iNum + find(tline(iNum:end)==' ')-2));
        entries = split(tline,'(');
        for iE = 2:length(entries)
            entries{iE} = replace(entries{iE},')','');
            i_spaces = find(entries{iE}==' ',3);
            if length(i_spaces==2);i_spaces(3) = length(entries{iE})+1;end
            probes.u(it,iE-1) = str2double(entries{iE}(1:i_spaces(1)-1));
            probes.v(it,iE-1) = str2double(entries{iE}(i_spaces(1)+1:i_spaces(2)-1));
            probes.w(it,iE-1) = str2double(entries{iE}(i_spaces(2)+1:i_spaces(3)-1));
        end
        it = it + 1;
    end

    % Get next line
    tline = fgetl(fID);
end
probes.phi = 270 - rad2deg(atan2(probes.v, probes.u));


%% Calculate output
% Get entries that fall within average time window
if ~isempty(probes.t)
    average_measurements = probes.t>= probes.t(end)-average_time;
else
    average_measurements = [];
end

if or(isempty(average_measurements), sum(average_measurements)==0)
    % No data avaiable -> Return default value
    T_phi = ones(1,probes.nProbes)*default_phi;
else
    % Return mean
    T_phi = mean(probes.phi(average_measurements,:));
end

end