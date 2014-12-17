function display4532(varargin)

%To find what line is printing to the command window, place a break point
%on the builtin line and run the code. Then click step, to see where in the
%code you were before.

%Code to insert breakpoint from command window
%dbstop in double/display4532 at 9

%This code caused errors when the function was named display
%https://www.mathworks.com/matlabcentral/newsreader/view_thread/93759
builtin('display', varargin{:});