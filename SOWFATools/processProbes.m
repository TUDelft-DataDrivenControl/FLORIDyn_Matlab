function processProbes(pathToProbes)
%pathToProbes = './Simulations/2021_9T_Data/SOWFA_Data/Probes';
probesDir = dir(pathToProbes);

nDir = length(probesDir);

for iD = 3:nDir
    if ~probesDir(iD).isdir
        continue;
    elseif ~startsWith(probesDir(iD).name,'probe','IgnoreCase',true)
        continue;
    end
    % Open probe file to read and create goal file to write
    fID_source  = fopen(...
        [pathToProbes '/' probesDir(iD).name '/20000/U'],'r');
    fID_goal    = fopen(...
        [pathToProbes '/' probesDir(iD).name '_U.csv'],'w');
    
    % Start reading lines from the probe file and use the data
    tline = fgetl(fID_source);
    nP = 0;
    while tline ~= -1
        if strcmp(tline(1:7),'# Probe')
            % Probe line
            % Can extract the number and location of the probes
            nP = nP + 1;
        elseif ~strcmp(tline(1),'#')
            % Data line
            t = str2double(erase(tline(1:15),' '));
            val = tline(16:end);
            val = replace(val,')               (',';');
            
            iS = 1;
            while strcmp(val(iS),' ')
                iS = iS+1;
            end
            val = val(iS+1:end-1);
            valNum = str2num(val);
            u_mean = sqrt(sum(mean(valNum).^2));
            fprintf(fID_goal,'%6.2f %12.8f\r\n',t,u_mean);
        end
        % Read next line, will equal -1 when file is over
        tline = fgetl(fID_source);
    end
    % Save changes and go to the next file
    fclose(fID_source);
    fclose(fID_goal);
end

