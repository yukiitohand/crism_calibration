function [x] = get_quadrantPxl(quadrant_id,varargin)
% [x] = get_quadrantPxl(quadrant_id)
%  get an array of spatial pixels for the input quadrant ID. In case of
%  binning = 0 (spatial averaging width = 1)
%   1: 1:160
%   2: 161:320
%   3: 321:480
%   4: 481:640
%   note: start from 1 not 0.
%
%  Optional Parameters
%    'BINNING' : ID of binning mode {0,1,2,3}
%                (default) 0
%    'BINX'    : binning size (PIXEL_AVERAGING_WIDTH in LBL)
%                (default) 1 
%     *Note: The relationship between BINNING and BINX
%        BINNING   BINX
%              0      1
%              1      2
%              2      5
%              3     10
%    If you only specify 'BINNING', 'BINX' is obtained automatically. If
%    you specify 'BINX', whatever 'BINNING' parameter will not be
%    considered.

binning = 0;
binx = 1;
if (rem(length(varargin),2)==1)
    error('Optional parameters should always go by pairs');
else
    for i=1:2:(length(varargin)-1)
        switch upper(varargin{i})
            case 'BINNING'
                binning = varargin{i+1};
                binx = get_binning(binning);
            case 'BINX'
                binx = varargin{i+1};
            otherwise
                % Hmmm, something wrong with the parameter string
                error(['Unrecognized option: ''' varargin{i} '''']);   
        end
    end
end

S = 640/binx;

switch quadrant_id
    case 1
        x = 1:S/4;
    case 2
        x = (S/4+1):S/2;
    case 3
        x = (S/2+1):(3*S/4);
    case 4
        x = ((3*S/4)+1):S;
    otherwise
        error('Quadrant%d is not defined',quadrant_id);
end

    