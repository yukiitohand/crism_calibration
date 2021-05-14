function [imgBI] = minipipeline_calibration_IR_BI_yuki(EDRBIdataList,PPdata,BSdata,HDdata,HKdata,DMdata,BPdata,varargin)
% [imgBI] = minipipeline_calibration_IR_BI_yuki(EDRBIdataList_s,PPdata,BSdata,HDdata,HKdata,DMdata,BPdata,varargin)
%  re-calculate CDR BI data from the collection of EDR BI data
%  INPUTS
%   EDRBIdataList: CRISMdata array of EDR BI data. Make sure frame rate is
%                  same for all the data in the list
%   PPdata: CDR PPdata
%   BSdata: CDR BSdata
%   HDdata: CDR HDdata
%   HKdata: CDR HKdata
%   DMdata: CDR DMdata (if 'BPRMVL' is set false, can be an empty)
%   BPdata: CDR BPdata (if 'BPRMVL' is set false, can be an empty)
%
%  OUTPUTS
%   imgBI: bias image [1 x S x B] (S: samples and B: bands)
%
%  OPTIONAL PARAMETERS
%   'DN4095_RMVL': binary, whether or not to perform replacement of saturated
%                  pixels or not.
%                  (default) false
%   'MEAN_ROBUST': integer {0,1}, mode for how mean operation is performed.
%        0: DN14a_i = nanmean(DN14(:,:,:),1);
%        1: DN14a_i = robust_v2('mean',DN14,1,'NOutliers',2);
%      (default) 1
%   'BPRMVL'     : binary, whether or not to perform bad pixel removal. 
%                  (default) false

dn4095_rmvl = false;
mean_robust = 1;
bprmvl = false;
if (rem(length(varargin),2)==1)
    error('Optional parameters should always go by pairs');
else
    for i=1:2:(length(varargin)-1)
        switch upper(varargin{i})
            case 'DN4095_RMVL'
                dn4095_rmvl = varargin{i+1};
            case 'BPRMVL'
                bprmvl = varargin{i+1};
            case 'MEAN_ROBUST'
                mean_robust = varargin{i+1};
            otherwise
                error('Undefined option: %s', varargin{i});
        end
    end
end

%-------------------------------------------------------------------------%
% check frame rate
% propEDRBIdataList = [EDRBIdataList.prop];
frame_rate_list = zeros(1,length(EDRBIdataList));
for i=1:length(EDRBIdataList)
    frame_rate_list(i) = EDRBIdataList(i).lbl.MRO_FRAME_RATE{1};
end
if length(unique(frame_rate_list)) > 1
    error('frame rate of the data does not seem to be unique.');
end

%-------------------------------------------------------------------------%
% get exposuretime
expo_timeList = zeros(length(EDRBIdataList),1);
integ_list = zeros(length(EDRBIdataList),1);
for i=1:length(EDRBIdataList)
    frame_rate = EDRBIdataList(i).lbl.MRO_FRAME_RATE{1};
    integ = EDRBIdataList(i).lbl.MRO_EXPOSURE_PARAMETER;
    integ_list(i) = integ;
    [t] = crism_get_integrationTime(integ,frame_rate,'Hz');
    expo_timeList(i) = t;
end

%-------------------------------------------------------------------------%
% covert raw 12bit biases to 14bits, then mean
DN14a = []; 
for i=1:length(EDRBIdataList)
    % read raw 12bit image -----------------------------------------------%
    DN12 = EDRBIdataList(i).readimg();
    % convert 12bits to 14bits -------------------------------------------%
    EDRBIdataList(i).read_ROWNUM_TABLE();
    [ DN14 ] = DN12toDN14( DN12,PPdata,EDRBIdataList(i).ROWNUM_TABLE );
    % flag saturated pixels ----------------------------------------------%
    if dn4095_rmvl
        flg_saturation = (DN12==4095);
        DN14(flg_saturation) = nan;
    end
    % take mean ----------------------------------------------------------%
    switch mean_robust
        case 0
            DN14a_i = nanmean(DN14(:,:,:),1);
        case 1
            DN14a_i = robust_v2('mean',DN14,1,'NOutliers',2);
        otherwise
            error('Undefined mean_robust=%d',mean_robust);
    end
    DN14a = cat(1,DN14a,DN14a_i);
end

%-------------------------------------------------------------------------%
% compute the term with a0I
[L,S,B] = size(DN14a);
if isempty(BSdata.tab),BSdata.readTAB(); end
term_a0IList = [];
for i=1:length(EDRBIdataList)
    binx = EDRBIdataList(i).lbl.PIXEL_AVERAGING_WIDTH;
    hkt = EDRBIdataList(i).readHKT();
    hkt = crism_correctHKTwithHD(hkt,HDdata);
    hkt = crism_correctHKTwithHK(hkt,HKdata);
    rate = cat(1,hkt.data.RATE);
    rate = unique(rate);
    if length(rate)>1
        error('something wrong with rate information in HKT of %s',EDRBIdataList(i).basename);
    end
    a0I = rateQuadrantTABformatter(rate,BSdata.tab,'A0','BINX',binx);
    row_lambdaList = reshape(EDRBIdataList(i).ROWNUM_TABLE,[1 1 B])+1; % +1 is already performed.
    integ_t = integ_list(i);
    term_integ_t = (502/480)*(480-integ_t);
    
    term_a0I = repmat(a0I,[1,1,B]).*heaviside(repmat(row_lambdaList,[1,S,1])-repmat(term_integ_t,[1,S,B]));
    term_a0IList = cat(1,term_a0IList,term_a0I);
end

%-------------------------------------------------------------------------%
% perform linear regression??
Y3 = DN14a+term_a0IList;
c0 = nan(1,S,B);
c1 = nan(1,S,B);
T = [ones(L,1) expo_timeList];
for b=1:B
    Y2 = Y3(:,:,b);
    chat = T\Y2;
    c0(1,:,b) = chat(1,:);
    c1(1,:,b) = chat(2,:);
end

imgBI = c0;

%-------------------------------------------------------------------------%
% bad pixel removal
if bprmvl
    [ imgBI,BP ] = apriori_badpixel_removal( imgBI,BPdata,BPdata,DMdata,'InterpOpt',1 );
end


end