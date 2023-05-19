function res = split(varargin)
% SPLIT - overload of both the classic MATLAB string split function and the 
% CORA 2018 split function. Based on the input arguments, this function
% chooses which of the two functions to call.
%
% The paths to the MATLAB and CORA directories need to be set individually
% to your machine.
%
% Author:  Lucas Lymburner
% Created: June 1, 2021

% Store current directory
orig_dir = cd;

% Check the input argument:
%   CORA 2018 split function only takes 2 arguments, neither of which are 
%   strings
arg0 = varargin{1};
if nargin ~= 2 || isstring(arg0) || ischar(arg0) || iscellstr(arg0)
    cd('/usr/local/MATLAB/R2023a/toolbox/matlab/strfun')
else
    cd('/home/roahm/Documents/MyToolBox/CORA_2018/global functions/globOptimization')
end

% Call the relevant split function after chaning to the correct directory
res = split(varargin{:});

% Change back to starting directory
cd(orig_dir);
end