% Make MATLAB set a breakpoint when there is an error
dbstop error;

% Set up figures with bigger text and
% bigger lines for better printing and slides
% use png() function to make picutres out of 
% your images
set(0,'DefaultAxesFontSize',14)
set(0,'DefaultAxesLineWidth',2)
set(0,'DefaultLineLineWidth',2)
set(0,'DefaultTextFontSize',14)
set(0,'DefaultAxesXGrid','On')
set(0,'DefaultAxesYGrid','On')

% Every engineer wants mor digits
format long

%% Get the previous path
% get the path of this script
path = fileparts(mfilename('fullpath'));
% look for curDir.mat next to this script
curDirFile = [path filesep 'curDir.mat'];
if exist(curDirFile,'file')
    load(curDirFile)
    if exist(curDir,'dir')
        cd (curDir);
    else
        fprintf('Could not cd to %s\n',curDir);
    end
end

clear path curDirFile curDir
