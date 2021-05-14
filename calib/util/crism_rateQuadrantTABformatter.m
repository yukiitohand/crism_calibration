function [datafrmtd] = crism_rateQuadrantTABformatter(rate_id,tab,colname,varargin)
% [datafrmtd] = crism_rateQuadrantTABformatter(rate_id,tab,colname,varargin)
%   format a column (colname) of the CRISM TAB data with columns defined 
%   as 'RATE','QUADRANT', and something, like
%         0,1,31.0
%         0,2,31.0
%         0,3,31.0
%         0,4,31.0
%         1,1,31.0
%         1,2,31.0
%         1,3,31.0
%         1,4,31.0
%         3,1,31.0
%         3,2,31.0
%         3,3,31.0
%         3,4,31.0
%         4,1,31.0
%         4,2,31.0
%         4,3,31.0
%         4,4,31.0
%
%  Input Parameters
%     rate_id: frame rate id at each line, [L x 1] or scalar
%     tab: struct of table data having two fields 'data' and 'colinfo'
%          the 'data' field is the data of the table having 'RATE' and 
%          'QUADRANT'
%     colname: field name of the target columns to be formatted to an image
%  Output Parameters
%     datafrmtd: [L x 640 x 1]
%                 tab.data.(colname) data formatted as an image.
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
%    Please specify only one of 'BINNING', and 'BINX'

binx = 1;
if (rem(length(varargin),2)==1)
    error('Optional parameters should always go by pairs');
else
    for i=1:2:(length(varargin)-1)
        switch upper(varargin{i})
            case 'BINNING'
                binning = varargin{i+1};
                binx = crism_get_binning(binning);
            case 'BINX'
                binx = varargin{i+1};
            otherwise
                error('Unrecognized option: %s', varargin{i});   
        end
    end
end

S = 640 / binx;

idxs = 4.*[tab.data.RATE;] + [tab.data.QUADRANT;];
tbl_srtd(idxs) = tab.data;

L = length(rate_id);
rate_id = rate_id(:);
datafrmtd = zeros([L,S,1]);
for i=1:4 % quadrant ID
    x = crism_get_quadrantPxl(i,'BINX',binx);
    idx = 4.*rate_id + i;
    datatmp = [tbl_srtd(idx).(colname)];
    datatmp = datatmp(:);
    datafrmtd(:,x,1) = repmat(datatmp,[1,length(x),1]);
end

end
