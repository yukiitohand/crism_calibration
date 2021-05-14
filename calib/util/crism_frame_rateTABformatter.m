function [datafrmtd] = frame_rateTABformatter(rate,tab,colname,varargin)
% [datafrmtd] = frame_rateTABformatter(rate,tab,colname)
%   format a column (colname) of the CRISM TAB data with columns defined 
%   as 'FRAME_RATE', and something, like
%     0,  0.00000, -0.0004, -0.00045,  100.0,   0.00
%     1,  0.00000, -0.0004, -0.00045,  100.0,   0.00
%     3,  0.00000, -0.0004, -0.00045,  100.0,   0.00
%     4,  0.00000, -0.0004, -0.00045,  100.0,   0.00
%
%  Input Parameters
%     rate: frame rate at each line, [L x 1] or scalar
%     tab: struct of table data having two fields 'data' and 'colinfo'
%          the 'data' field is the data of the table having 'FRAME_RATE'
%     colname: field name of the target columns to be formatted to an image
%  Optional Parameters
%    'COLUMNS' : number of columns to be generated 
%                (default) 640
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
%    Please specify only one of 'COLUMNS','BINNING', and 'BINX'
%  Output Parameters
%     datafrmtd: [L x 'COLUMNS' x 1]
%                 tab.data.(colname) data formatted as an image.

columns = 640;
if (rem(length(varargin),2)==1)
    error('Optional parameters should always go by pairs');
else
    for i=1:2:(length(varargin)-1)
        switch upper(varargin{i})
            case 'BINNING'
                binning = varargin{i+1};
                binx = get_binning(binning);
                columns = columns/binx;
            case 'BINX'
                binx = varargin{i+1};
                columns = columns/binx;
            case 'COLUMNS'
                columns = varargin{i+1};
            otherwise
                % Hmmm, something wrong with the parameter string
                error(['Unrecognized option: ''' varargin{i} '''']);   
        end
    end
end
idxs = [tab.data.FRAME_RATE;]+1;
tbl_srtd(idxs) = tab.data;
datatmp = [tbl_srtd(rate+1).(colname)];
datatmp = datatmp(:);
datafrmtd = repmat(datatmp,[1,columns,1]);

end