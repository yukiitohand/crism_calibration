function [BIdata_o,imgBI] = minipipeline_calibration_IR_BI_wCDRBI_yuki(CDRBIdata,varargin)
% [BIdata_o,imgBI] = minipipeline_calibration_IR_BI_wCDRBI_yuki(CDRBIdata,varargin)
%  re-calculate CDR BI data without performing bad pixel interpolation
%  INPUTS
%    CDRBIdata: CRISMdata obj, reference CDR BI data.
%  OUTPUTS
%    BIdata_o : CRISMdata obj, basname is pointed to the original CDRBIdata
%               storing imgBI at the attribute img.
%    imgBI    : bias image [1 x S x B] (S: samples and B: bands)
%  Optional Parameters
%   'DN4095_RMVL': binary, whether or not to perform replacement of saturated
%                  pixels or not.
%                  (default) false
%   'MEAN_ROBUST': integer {0,1}, mode for how mean operation is performed.
%        0: DN14a_i = nanmean(DN14(:,:,:),1);
%        1: DN14a_i = robust_v2('mean',DN14,1,'NOutliers',2);
%      (default) 1
%   'BPRMVL'     : binary, whether or not to perform bad pixel removal. 
%                  (default) false
%   'DWLD','DOWNLOAD' : if download the data or not, 2: download, 1:
%                       access an only show the path, 0: nothing
%                       (default) 0
%   'OUT_FILE'       : path to the output file
%                       (default) ''
%   'Force'          : binary, whether or not to force performing
%                      pds_downloader. (default) false

dn4095_rmvl = false;
mean_robust = 1;
bprmvl = false;
dwld = 0;
force = false;
outfile = '';

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
            case {'DWLD','DOWNLOAD'}
                dwld = varargin{i+1};
            case 'FORCE'
                force = varargin{i+1};
            case 'OUT_FILE'
                outfile = varargin{i+1};
        end
    end
end

% if isempty(CDRBIdata.basenamesCDR), CDRBIdata.load_basenamesCDR(); end
% if isempty(CDRBIdata.basenames_SOURCE_OBS)
%     CDRBIdata.load_basenames_SOURCE_OBS(); 
% end
CDRBIdata.load_basenamesCDR('Download',dwld,'Force',force,'OUT_file',outfile);
CDRBIdata.load_basenames_SOURCE_OBS('Download',dwld,'Force',force,'OUT_file',outfile); 

% get source BIdata download if not exists.
EDRBIdataList = CDRBIdata.read_SOURCE_OBS('BI');
% for i=1:length(CDRBIdata.basenames_SOURCE_OBS.BI)
%     basename_edrbi = CDRBIdata.basenames_SOURCE_OBS.BI{i};
%     [dirfullpath_local_obsbi,~,~,~,~]...
%         = get_dirpath_observation(basename_edrbi,'Download',dwld,...
%         'Force',force,'OUT_file',outfile);
%     bidata = CRISMdata(basename_edrbi,dirfullpath_local_obsbi);
%     EDRBIdataList = [EDRBIdataList bidata];
% end

% read HKP
for i=1:length(EDRBIdataList)
    EDRBIdataList(i).download(dwld);
end


% select ones with the same frame rate.
frame_rate = CDRBIdata.lbl.MRO_FRAME_RATE{1};
EDRBIdataList_s = [];
for i=1:length(EDRBIdataList)
    if EDRBIdataList(i).lbl.MRO_FRAME_RATE{1} == frame_rate
        EDRBIdataList_s = [EDRBIdataList_s EDRBIdataList(i)];
    end
end

PPdata = CDRBIdata.readCDR('PP');
BSdata = CDRBIdata.readCDR('BS');
HDdata = CDRBIdata.readCDR('HD');
HKdata = CDRBIdata.readCDR('HK');

% BPdata and DMdata are only necessary when bprmvl is on.
    function [BPdata] = get_BPdata_fromEDRBI(EDRBIdataList_s,dwld)
        sclk_EDRBI_list = zeros(1,length(EDRBIdataList_s));
        for ii=1:length(EDRBIdataList_s)
            [sclk_stop_i] = EDRBIdataList_s(ii).get_sclk_stop();
            sclk_EDRBI_list(ii) = sclk_stop_i;
        end

        sclk_stop = max(sclk_EDRBI_list);
        propBP_search = crism_create_propCDR4basename;
        propBP_search.acro_calibration_type = 'BP';
        propBP_search.binning = crism_get_binning_id(EDRBIdataList_s(1).lbl.PIXEL_AVERAGING_WIDTH);
        propBP_search.wavelength_filter = EDRBIdataList_s(1).lbl.MRO_WAVELENGTH_FILTER;
        propBP_search.frame_rate = crism_get_frame_rate_id(EDRBIdataList_s(1).lbl.MRO_FRAME_RATE{1});
        % for the BP products, exposure is almost always 0. Do not set it.
        % propBP_search.exposure = DFdata.lbl.MRO_EXPOSURE_PARAMETER;
        propBP_search.sensor_id = EDRBIdataList_s(1).lbl.MRO_SENSOR_ID;
        propBP_search.sclk = floor(sclk_stop);
        [basenameBPmrb,propBPmrb] = crism_searchCDRmrb(propBP_search,'dwld',1,'force',1);
        crism_get_dirpath_cdr(basenameBPmrb,'dwld',dwld);
        BPdata = CRISMdata(basenameBPmrb,'');
        for ii=1:length(EDRBIdataList_s)
            basenameEDRBI_i = EDRBIdataList_s(ii).basename;
            if isempty(extractMatchedBasename_v2(basenameEDRBI_i,BPdata.lbl.SOURCE_PRODUCT_ID))
                error('It seems no BPdata associated with %s.',basenameEDRBI_i);
            end
        end
    end
if bprmvl
    % find BPdata
    BPdata = get_BPdata_fromEDRBI(EDRBIdataList_s,dwld);
    BPdata.load_basenamesCDR();
    DMdata = BPdata.readCDR('DM');
else
    DMdata = []; BPdata = [];
end

%--------------------------------------------------------------------------
% main calculation
[imgBI] = minipipeline_calibration_IR_BI_yuki(...
    EDRBIdataList_s,PPdata,BSdata,HDdata,HKdata,DMdata,BPdata,...
    'DN4095_RMVL',dn4095_rmvl,'BPRMVL',bprmvl,'MEAN_ROBUST',mean_robust);
BIdata_o = CRISMdata(CDRBIdata.basename,'');
BIdata_o.img = imgBI;

% % get exposuretime
% expo_timeList = zeros(length(EDRBIdataList_s),1);
% for i=1:length(EDRBIdataList_s)
%     integ = EDRBIdataList_s(i).lbl.MRO_EXPOSURE_PARAMETER;
%     [t] = crism_get_integrationTime(integ,frame_rate,'Hz');
%     expo_timeList(i) = t;
% end
% 
% % covert raw 12bit biases to 14bits, then mean
% DN14a = [];
% PPdata = CDRBIdata(i).readCDR('PP');
% for i=1:length(EDRBIdataList_s)
%     EDRBIdataList_s.read_ROWNUM_TABLE();
%     DN12 = EDRBIdataList_s.readimg();
%     [ DN14 ] = crmcal_DN12toDN14( DN12,PPdata,EDRBIdataList_s.rownum_table );
%     DN14a_i = robust_v2('mean',DN14,1,'NOutliers',2);
%     DN14a = cat(3,DN14a,DN14a_i);
% end
% 
% % compute the term with a0I
% [L,S,B] = size(DN14a);
% BSdata = CDRBIdata.readCDR('BS');
% if isempty(BSdata.tab),BSdata.readTAB(); end
% binx = CDRBIdata.lbl.PIXEL_AVERAGING_WIDTH;
% a0I = crism_rateQuadrantTABformatter(rate,BSdata.tab,'A0','BINX',binx);
% row_lambdaList = reshape(rownum_table,[1 1 B])+1; % +1 is already performed.
% term_a0IList = [];
% for i=1:length(EDRBIdataList_s)
%     integ_t = expo_timeList(i);
%     term_integ_t = (502/480)*(480-integ_t);
%     term_a0I = repmat(a0I,[1,1,B]) ...
%     .*heaviside(repmat(row_lambdaList,[1,S,1])-repmat(term_integ_t,[1,S,B]));
%     term_a0IList = cat(3,term_a0IList,term_a0I);
% end
% 
% % perform linear regression??
% Y3 = DN14a+term_a0IList;
% c0 = nan(1,S,B);
% c1 = nan(1,S,B);
% T = [ones(L,1) expo_timeList];
% for b=1:B
%     Y2 = Y3(:,:,B);
%     chat = T\Y2;
%     c0(1,:,b) = chat(1,:);
%     c1(1,:,b) = chat(2,:);
% end
% 
% imgBI = c0;



end
