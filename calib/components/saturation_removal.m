function [DN14c,mask_saturation,mask_dead] = saturation_removal(DN14b,VLdata,mask4095,varargin)
% [DN14c,mask_saturation] = saturation_removal(DN14b,VLdata,mask4095)
%   mask and replace "digital" and "analogue" saturated pixels with nans
%  INPUTS
%    DN14b: DN14b image [LxSxB]
%    VLdata: CRISMdata obj, CDR VL data
%    mask4095: boolean image [LxSxB] of "digital" saturated pixels
%  Optional Parameters
%    'RATE' (default) 1, frame rate id
%   rate  rateHz[Hz]
%      0       1
%      1    3.75
%      2       5
%      3      15
%      4      30
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
%    'IS_SPHERE' boolean, (default) false
%  OUTPUTS
%    DN14c: image processed, saturated pixels are replaced with nans.
%    mask_saturation: [LxSxB] boolean image indicating the saturated pixels
%                     with 1s.
%    mask_dead: [LxSxB] boolean image indicating dead pixels. sphere
%               measurement only.

if isempty(VLdata.tab), VLdata.readTAB(); end

rate_id = 1;
binx = 1;
is_sphere = false;
if (rem(length(varargin),2)==1)
    error('Optional parameters should always go by pairs');
else
    for i=1:2:(length(varargin)-1)
        switch upper(varargin{i})
            case 'RATE_ID'
                rate_id = varargin{i+1};
            case 'BINNING'
                binning = varargin{i+1};
                binx = get_binning(binning);
            case 'BINX'
                binx = varargin{i+1};
            case 'IS_SPHERE'
                is_sphere = varargin{i+1};
            otherwise
                % Hmmm, something wrong with the parameter string
                error(['Unrecognized option: ''' varargin{i} '''']);   
        end
    end
end

[L,S,B] = size(DN14b);
DN14c = DN14b;

% first "digital saturation" is removed.
if ~isempty(mask4095)
    DN14c(mask4095) = nan;
end


IR_14_BIT_LIMIT = rateQuadrantTABformatter(rate_id,VLdata.tab,'IR_14_BIT_LIMIT','BINX',binx);
IR_LINEARITY_LIMIT = rateQuadrantTABformatter(rate_id,VLdata.tab,'IR_LINEARITY_LIMIT','BINX',binx);
IR_NOISE_LIMIT = rateQuadrantTABformatter(rate_id,VLdata.tab,'IR_NOISE_LIMIT','BINX',binx);
IR_SENSITIVITY_LIMIT = rateQuadrantTABformatter(rate_id,VLdata.tab,'IR_SENSITIVITY_LIMIT','BINX',binx);

% Next "analogue saturation" is removed.
flg_asat = DN14>=(IR_14_BIT_LIMIT.*IR_LINEARITY_LIMIT);
DN14c(flg_asat) = nan;

% next detect dead pixels
if is_sphere
    DN14c_ext = reshape(DN14c,[L*S,B]);
    DN14c_med = reshape(nanmedian(DN14c_ext,1),[1,1,B]);
    mask_dead = DN14c < (IR_SENSITIVITY_LIMIT.*DN14c_med);
    DN14c(mask_dead) = nan;
else
    mask_dead = '';
end

if ~isempty(mask4095)
    mask_saturation = or(mask4095,flg_asat);
else
    mask_saturation = flg_asat;
end

end