function [RSPl] = binning_RSP(RSPj,DMdata,rowNumTableRSPj,varargin)
% [RSPl] = binning_RSP(RSPj,DMdata)
%  apply binning to spectralradiometric responsitivity using the detector
%  mask defined in (CDR DM data).
%   Input
%     RSPj   : spectral radiometric response
%     DMdata : CRISM obj, CDR DM data
%   Output Parameters
%     RSPl: binned spectral radiometric responsitivity
%   Optional Parameters
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
% [~,S,Band] = size(RSPj);

if isempty(DMdata.img), DMdata.readimg(); end
if isempty(DMdata.ROWNUM_TABLE), DMdata.read_ROWNUM_TABLE(); end

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
                error(['Unrecognized option: %s', varargin{i});   
        end
    end
end


[~,S_b,Bands_b] = size(DMdata.img);

S_rspl = crism_getSampleSize(binx);
S_rspj = size(RSPj,2);

if S_rspl==S_rspj
    RSPl = RSPj;
else
    if all(rowNumTableRSPj==DMdata.ROWNUM_TABLE)
        RSPk = squeeze(RSPj); %[S x L]

        h = fspecial('average',[binx,1]); % convolution vector
        RSPk_b = conv2(RSPk,h,'valid'); % take convolution and 
        RSPk_b = RSPk_b(1:binx:end,:); % downsample

        RSPl = reshape(RSPk_b,[1,S_b,Bands_b]);

    else
        error('not implemented for wavelengthfilter > 0');
    end
end