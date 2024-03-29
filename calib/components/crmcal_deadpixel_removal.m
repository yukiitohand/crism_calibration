function [RT14jj,mask_dead] = crmcal_deadpixel_removal(RT14j,VLdata,DMdata,varargin)
% [RT14jj,mask_dead] = crmcal_deadpixel_removal(RT14j,VLdata,DMdata,varargin)
%   Mask and replace dead pixels with nans. Sphere measurement only. Only
%   scene pixels are evaluated.
%  INPUTS
%    RT14j: RT14j image [LxSxB]
%    VLdata: CRISMdata obj, CDR VL data
%    DMdata: CDR DM data
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
%  OUTPUTS
%    RT14jj: image processed, dead pixels are replaced with nans. non-scene
%    pixels are untouched.
%    mask_dead: [LxSxB] boolean image indicating dead pixels. 

if isempty(VLdata.tab), VLdata.readTAB(); end

rate_id = 1;
binx = 1;
if (rem(length(varargin),2)==1)
    error('Optional parameters should always go by pairs');
else
    for i=1:2:(length(varargin)-1)
        switch upper(varargin{i})
            case 'RATE_ID'
                rate_id = varargin{i+1};
            case 'BINNING'
                binning = varargin{i+1};
                binx = crism_get_binx(binning);
            case 'BINX'
                binx = varargin{i+1};
            otherwise
                error('Unrecognized option: %s', varargin{i});   
        end
    end
end

[L,S,B] = size(RT14j);
RT14jj = RT14j;

% IR_NOISE_LIMIT = rateQuadrantTABformatter(rate_id,VLdata.tab,'IR_NOISE_LIMIT','BINX',binx);
IR_SENSITIVITY_LIMIT = crism_rateQuadrantTABformatter(rate_id,VLdata.tab,'IR_SENSITIVITY_LIMIT','BINX',binx);


% next detect dead pixels
RT14j_dm     = crmcal_apply_DM(RT14j,DMdata);
% RT14j_dm_ext = reshape(RT14j_dm,[L*S,B]);
% RT14j_dm_med = reshape(nanmedian(RT14j_dm_ext,1),[1,1,B]);
RT14j_dm_med = median(RT14j_dm,[1,2],'omitnan');
mask_dead = RT14j_dm < (IR_SENSITIVITY_LIMIT.*RT14j_dm_med);
RT14jj(mask_dead) = nan;



end