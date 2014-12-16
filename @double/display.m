function display(varargin)

%To find what line is printing to the command window, place a break point
%on the builtin line and run the code. Then click step, to see where in the
%code you were before.

%Code to insert breakpoint from command window
%dbstop in double/display at 9
builtin('display', varargin{:});