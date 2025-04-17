function p = xPath(varargin)
% xPath  Convert a single absolute path or build one from pieces.
%        – If you pass ONE string containing slashes, it converts it.
%        – If you pass ≥2 args, it builds a native path with fullfile.

if nargin==1                                 % absolute → convert
    s = varargin{1};

    if ispc                                   % need Windows form
        s = regexprep(s,'^[\/]','');          % drop leading “/” (if any)
        p = strrep(s,'/','\');                % flip slashes

    else                                      % need Unix/mac form
        s = regexprep(s,'^[A-Za-z]:\\','');   % drop drive (C:\ …)
        p = strrep(s,'\','/');                % flip slashes
        if ~startsWith(p,'/'), p = ['/' p]; end
    end

else                                          % build from pieces
    p = fullfile(varargin{:});                % native separator
end
end
