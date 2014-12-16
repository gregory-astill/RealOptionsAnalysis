
%% Save current directory
% Get the current directory
curDir = pwd;
% Get the path of startup
path = fileparts(which('startup'));
%Save the current directory to a file
curDirFile = [path filesep 'curDir.mat'];
save(curDirFile,'curDir');

