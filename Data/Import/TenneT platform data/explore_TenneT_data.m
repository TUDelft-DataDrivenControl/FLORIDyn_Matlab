dirs = dir("/Users/marcusbecker/surfdrive/PhD_Surf/01_Research/01_FLORIDyn/20_Data/06_TenneT_platform");

%% TODO 
% 1. Create reduced files
% 2. Replace 9999 values
% 3. Smooth data
% 4. 


%%
figure
hold on
for iD = 3:length(dirs)
    filename = [dirs(iD).folder filesep dirs(iD).name];
    time_dir_uv_199m = importfile_ZephIR(filename);
    plot(time_dir_uv_199m.TimeAndDate,time_dir_uv_199m.WindDirectiondegAt119m)
end

%% Figure




%%
figure()
plot(cos(deg2rad(time_dir_uv_199m.WindDirectiondegAt119m)),...
    sin(deg2rad(time_dir_uv_199m.WindDirectiondegAt119m)))