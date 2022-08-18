function [SPdata_o,RT14j_woc,RT14j,RT14h2_bk1_o,RT14h2_bk2_o,BPdata1,BPdata2] = crmcal_minipipeline_IR_SP_wCDRSP_yuki(...
    SPdata,bkoption,varargin)
% [SPdata_o,RT14j_woc,RT14j,RT14h2_bk1_o,RT14h2_bk2_o,BPdata1,BPdata2] = crmcal_minipipeline_IR_SP_wCDRSP_yuki(...
%     SPdata,bkoption,varargin)
%   Mini pipeline for the calibration of the CRISM images using SP CDR
%   data. This is a wrapper function for 
%    "minipipeline_calibration_IR_SP_yuki"
%
%  Input Parameters
%   SPdata: CDR SPdata, CRISMdata obj
%   bkoption: option for the use of background image. {1,2}
%             1: linear background estimation using prior and post dark
%                measurements.
%             2: flat background estimation using only prior dark
%                measurements.
%   OUTPUTS
%    SPdata_o: SPdata that stores processed image at img. Band inverse is
%               not performed. RT14j_woc is stored.
%    RT14j_woc: non interpolated version of SPdata
%    RT14j: interpolated version of SPdata
%    RT14h2_bk1_o: non interpolated version of processed dark frame
%   Optional Parameters
%    'SAVE_MEMORY'
%       saving memory or not. true or false
%       (default) true
%    'DWLD','DOWNLOAD' : if download the data or not, 2: download, 1:
%                       access an only show the path, 0: nothing
%                       (default) 0
%    'Force_dwld'     : binary, whether or not to force performing
%                      pds_downloader. (default) false
%   'DWLD_INDEX_CACHE_UPDATE' : boolean, whether or not to update index.html 
%        (default) false
%   'VERBOSE_DWLD'   : boolean, whether or not to show the downloading
%                      operations.
%                      (default) true
%   'DWLD_OVERWRITE' : if overwrite the file if exists
%                      (default) 0
%    'APBPRMVL'
%       option for a priori bad pixel removal {'HighOrd', 'None'}
%       'HighOrd': the image where a priori bad pixel removal is performed
%                  is used for the estimation of the higher order leaked 
%                  light
%       'None'   : the image where a priori bad pixel removal is NOT performed
%                  is used for the estimation of the higher order leaked 
%                  light
%       (default) 'HighOrd'
%    'MEAN_ROBUST' : integer {0,1}, mode for how mean operation is performed.
%        0: DN14e_df = nanmean(DN14d_df(:,:,:),1);
%        1: DN14e_df = robust_v2('mean',DN14d_df,1,'NOutliers',2);
%      (default) 1
%    'SATURATiON_RMVL': integer, how to perform replacement of saturated
%           pixles {0,1,2}
%           0: no removal
%           1: digital saturation is removed
%           2: analogue saturation is also removed
%           (default) 2
%    'MEAN_DN14'  : binary,when mean operation is performed
%                  1: before non-linearity correction
%                  0: last (after divided by integration time
%                  (default) 1
%   ****** Parameters for manual BK production ****************************
%   'BK_SATURATiON_RMVL': integer, how to perform replacement of saturated
%           pixles {0,1,2}
%           0: no removal
%           1: digital saturation is removed
%           2: analogue saturation is also removed
%           (default) 2
%   'BK_MEAN_ROBUST': integer {0,1}, mode for how mean operation is performed.
%        0: DN14e_df = nanmean(DN14d_df(:,:,:),1);
%        1: DN14e_df = robust_v2('mean',DN14d_df,1,'NOutliers',2);
%      (default) 1
%   'BK_BPRMVL'     : binary, whether or not to perform bad pixel removal. 
%                  (default) false     
%   'BK_MEAN_DN14'  : binary,when mean operation is performed
%                  1: before non-linearity correction
%                  0: last (after divided by integration time
%                  (default) 1
save_mem = false;
apbprmvl = 'HighOrd';
saturation_rmvl = 2;
mean_DN14 = true;
mean_robust = true;
bk_saturation_rmvl = 2;
bk_mean_robust = 1;
bk_bprmvl = false;
bk_mean_DN14 = true;
dwld = 0;
force_dwld = false;
dwld_index_cache_update = 0;
verbose_dwld = 1;
dwld_overwrite = 0;
BIdata = [];

if (rem(length(varargin),2)==1)
    error('Optional parameters should always go by pairs');
else
    for i=1:2:(length(varargin)-1)
        switch upper(varargin{i})
            case 'SAVE_MEMORY'
                save_mem = varargin{i+1};
            case 'APBPRMVL'
                apbprmvl = varargin{i+1};
                if ~any(strcmpi(apbprmvl,{'HighOrd','None'}))
                    error('apbprmvl (%s) should be either {"HighOrd","None"}',apbprmvl);
                end
            case 'MEAN_DN14'
                mean_DN14 = varargin{i+1};
            case 'MEAN_ROBUST'
                mean_robust = varargin{i+1};
            case 'SATURATION_RMVL'
                saturation_rmvl = varargin{i+1};
            case 'BK_SATURATION_RMVL'
                bk_saturation_rmvl = varargin{i+1};
            case 'BK_BPRMVL'
                bk_bprmvl = varargin{i+1};
            case 'BK_MEAN_ROBUST'
                bk_mean_robust = varargin{i+1};
            case 'BK_MEAN_DN14'
                bk_mean_DN14 = varargin{i+1};
            case {'DWLD','DOWNLOAD'}
                dwld = varargin{i+1};
            case 'FORCE_DWLD'
                force_dwld = varargin{i+1};
            case 'DWLD_INDEX_CACHE_UPDATE'
                dwld_index_cache_update = varargin{i+1};
            case 'DWLD_OVERWRITE'
                dwld_overwrite = varargin{i+1};
            case 'VERBOSE_DWLD'
                verbose_dwld = varargin{i+1};
            case 'BIDATA'
                BIdata = varargin{i+1};
            otherwise
                error('Unrecognized option: %s', varargin{i});
        end
    end
end

binx_sp = SPdata.lbl.PIXEL_AVERAGING_WIDTH;
binning_id_sp = crism_get_binning_id(binx_sp);

% if isempty(SPdata.basenamesCDR)
    SPdata.load_basenamesCDR('Download',dwld, ...
        'OVERWRITE',dwld_overwrite,'INDEX_CACHE_UPDATE',dwld_index_cache_update);
% end
% if isempty(SPdata.basenames_SOURCE_OBS)
    SPdata.load_basenames_SOURCE_OBS('Download',dwld, ...
        'OVERWRITE',dwld_overwrite,'INDEX_CACHE_UPDATE',dwld_index_cache_update);
% end

%-------------------------------------------------------------------------%
% get EDRSPdata from SPdata
basenameEDRSP = SPdata.basenames_SOURCE_OBS.SP;
if iscell(basenameEDRSP)
    error('Mulitple EDR SP is used. not supported.');
end
EDRSPdata = CRISMdata(basenameEDRSP,'');
EDRSPdata.download(dwld);

%-------------------------------------------------------------------------%
% get DFdata from SPdata
DFdata1 = get_DFdata(EDRSPdata,sprintf('%02s',dec2hex(hex2dec(EDRSPdata.prop.obs_counter)-1)),dwld);
DFdata1.download(dwld);
DFdata2 = get_DFdata(EDRSPdata,sprintf('%02s',dec2hex(hex2dec(EDRSPdata.prop.obs_counter)+1)),dwld);
DFdata2.download(dwld);

%-------------------------------------------------------------------------%
% get BPdata using sclk of the DFdata
BPdata1 = get_BPdata_fromDF(DFdata1,binning_id_sp,dwld);
BPdata2 = get_BPdata_fromDF(DFdata2,binning_id_sp,dwld);

%-------------------------------------------------------------------------%
% get BIdata using source_obs of SPdata
EDRBIdataList = SPdata.read_SOURCE_OBS('BI');

% select EDR BI data with same frame rate.
frame_rate = SPdata.lbl.MRO_FRAME_RATE.value;
EDRBIdataList_s = CRISMdata.empty(0,0);
for i=1:length(EDRBIdataList)
    if ~isempty(EDRBIdataList(i).lbl) && EDRBIdataList(i).lbl.MRO_FRAME_RATE.value == frame_rate
        EDRBIdataList_s = [EDRBIdataList_s EDRBIdataList(i)];
    end
end

if isempty(EDRBIdataList_s)
    error('There is no EDR BI data with the same frame-rate as SP data.');
end

% sometimes two or more different sets of BI EDR is selected.
edrbi_dirnames_s = unique({EDRBIdataList_s.dirname});
edrbi_groups_s = cellfun(@(x) find(strcmpi(x,edrbi_dirnames_s)), {EDRBIdataList_s.dirname});
sclk_edrbi_groups = zeros(1,length(edrbi_dirnames_s));
for gi=1:length(edrbi_dirnames_s)
    idxes = find(gi==edrbi_groups_s);
    sclk_mean_list = nan(1,length(idxes));
    for iii = 1:length(idxes)
        [sclk_stop_i] = EDRBIdataList_s(idxes(iii)).get_sclk_stop();
        [sclk_start_i] = EDRBIdataList_s(idxes(iii)).get_sclk_start();
        sclk_mean_list(iii) = (sclk_start_i + sclk_stop_i) / 2;
    end
    sclk_edrbi_groups(gi) = nanmean(sclk_mean_list);
end
[~,edrbi_closest] = min(abs(SPdata.prop.sclk-sclk_edrbi_groups));
EDRBIdataList_ss = EDRBIdataList_s(edrbi_groups_s==edrbi_closest);

if isempty(BIdata)
    [BIdata] = get_BIdata_fromEDRBI(EDRBIdataList_ss,binning_id_sp,dwld);
end

%-------------------------------------------------------------------------%
% get BKdata using sclk of the DFdata
BKdata1 = get_BKdata_fromDF(DFdata1,binning_id_sp,dwld);
BKdata2 = get_BKdata_fromDF(DFdata2,binning_id_sp,dwld);

%-------------------------------------------------------------------------%
% Read othe CDRs
PPdata = SPdata.readCDR('PP');
BSdata = SPdata.readCDR('BS'); DBdata = SPdata.readCDR('DB');
EBdata = SPdata.readCDR('EB'); HDdata = SPdata.readCDR('HD');
HKdata = SPdata.readCDR('HK');
GHdata = SPdata.readCDR('GH');
LCdata = SPdata.readCDR('LC');
DMdata = SPdata.readCDR('DM');
LLdata = SPdata.readCDR('LL');
VLdata = SPdata.readCDR('VL');

%-------------------------------------------------------------------------%
% main pipeline
[SPdata_o,RT14j_woc,RT14j,RT14h2_bk1_o,RT14h2_bk2_o] = crmcal_minipipeline_IR_SP_yuki(...
    EDRSPdata,DFdata1,DFdata2,BKdata1,BKdata2,BPdata1,BPdata2,BIdata,...
    PPdata,BSdata,DBdata,EBdata,HDdata,HKdata,GHdata,VLdata,DMdata,LCdata,LLdata,...
    bkoption,'SAVE_MEMORY',save_mem,'APBPRMVL',apbprmvl,'MEAN_DN14',mean_DN14,...
    'SATURATION_RMVL',saturation_rmvl,'Mean_Robust',mean_robust,...
    'BK_SATURATION_RMVL',bk_saturation_rmvl,'BK_BPRMVL',bk_bprmvl,...
    'BK_MEAN_ROBUST',bk_mean_robust,'BK_MEAN_DN14',bk_mean_DN14,'SPdata_ref',SPdata,...
    'BINNING_SP',binning_id_sp);

% [BP1nan] = crism_formatBP1nan(BPdata1);
% BPpri1nan = crism_formatBPpri1nan(BPdata1,BPdata2);
%[BIdata_o,imgBI] = crmcal_minipipeline_IR_BI_wCDRBI_yuki(BIdata,'DN4095_RMVL',0,'BPRMVL',0,'MEAN_ROBUST',1);
end

function [DFdata] = get_DFdata(EDRSPdata,obs_counter,dwld)
    propEDRDF = crism_create_propOBSbasename();
    propEDRDF.obs_counter = obs_counter;
    propEDRDF.obs_class_type = EDRSPdata.prop.obs_class_type;
    propEDRDF.obs_id = EDRSPdata.prop.obs_id;
    propEDRDF.product_type = EDRSPdata.prop.product_type;
    propEDRDF.sensor_id = EDRSPdata.prop.sensor_id;
    propEDRDF.activity_id = 'DF';

    [dir_info,basenameEDRDF] = crism_search_observation_fromProp(propEDRDF,'dwld',dwld);

    DFdata = CRISMdata(basenameEDRDF,'');
    DFdata.download(dwld);
end

function [BPdata] = get_BPdata_fromDF(DFdata,binning_id_sp,dwld)
    propBP_search = crism_create_propCDR4basename;
    propBP_search.acro_calibration_type = 'BP';
    propBP_search.binning = crism_get_binning_id(DFdata.lbl.PIXEL_AVERAGING_WIDTH);
    propBP_search.wavelength_filter = DFdata.lbl.MRO_WAVELENGTH_FILTER;
    propBP_search.frame_rate = crism_get_frame_rate_id(DFdata.lbl.MRO_FRAME_RATE.value);
    % for the BP products, exposure is almost always 0. Do not set it.
    % propBP_search.exposure = DFdata.lbl.MRO_EXPOSURE_PARAMETER;
    propBP_search.sensor_id = DFdata.lbl.MRO_SENSOR_ID;
    [sclk_df_stop,p_df_stop] = DFdata.get_sclk_stop();
    propBP_search.sclk = floor(sclk_df_stop);
    propBP_search.partition = p_df_stop;
    propBP_search.binning = binning_id_sp;
    [basenameBPmrb,propBPmrb] = crism_searchCDRmrb(propBP_search,'dwld',1);
    crism_get_dirpath_cdr(basenameBPmrb,'dwld',dwld);
    BPdata = CRISMdata(basenameBPmrb,'');
    if isempty(extractMatchedBasename_v2(DFdata.basename,BPdata.lbl.SOURCE_PRODUCT_ID))
        fprintf('It seems no BPdata associated with %s.\n in local database',DFdata.basename);
        fprintf('Connect to remote server...\n');
        [basenameBPmrb,propBPmrb] = crism_searchCDRmrb(propBP_search,'dwld',1);
        crism_get_dirpath_cdr(basenameBPmrb,'dwld',dwld);
        BPdata = CRISMdata(basenameBPmrb,'');
        if isempty(extractMatchedBasename_v2(DFdata.basename,BPdata.lbl.SOURCE_PRODUCT_ID))
            error('It seems no BPdata associated with %s.',DFdata.basename);
        end
    end
end

function [BIdata] = get_BIdata_fromEDRBI(EDRBIdataList_s,binning_id_sp,dwld)
    % get sclk
    sclk_EDRBI_list = zeros(1,length(EDRBIdataList_s));
    for ii=1:length(EDRBIdataList_s)
        [sclk_stop_i] = EDRBIdataList_s(ii).get_sclk_stop();
        sclk_EDRBI_list(ii) = sclk_stop_i;
    end

    sclk_stop = max(sclk_EDRBI_list);

    propBI_search = crism_create_propCDR4basename;
    propBI_search.acro_calibration_type = 'BI';
    propBI_search.binning = crism_get_binning_id(EDRBIdataList_s(1).lbl.PIXEL_AVERAGING_WIDTH);
    propBI_search.wavelength_filter = EDRBIdataList_s(1).lbl.MRO_WAVELENGTH_FILTER;
    propBI_search.frame_rate = crism_get_frame_rate_id(EDRBIdataList_s(1).lbl.MRO_FRAME_RATE.value);
    propBI_search.sensor_id = EDRBIdataList_s(1).lbl.MRO_SENSOR_ID;
    propBI_search.sclk = floor(sclk_stop);
    propBI_search.binning = binning_id_sp;
    [basenameBImrb,propBImrb] = crism_searchCDRmrb(propBI_search,'dwld',1);
    crism_get_dirpath_cdr(basenameBImrb,'dwld',dwld);
    BIdata = CRISMdata(basenameBImrb,'');

    is_in_local_database = 1;
    for ii=1:length(EDRBIdataList_s)
        basenameEDRBI_i = EDRBIdataList_s(ii).basename;
        if isempty(extractMatchedBasename_v2(basenameEDRBI_i,BIdata.lbl.SOURCE_PRODUCT_ID))
            fprintf('It seems no BIdata associated with %s.\n',basenameEDRBI_i);
            is_in_local_database = 0;
        end
    end
    if ~is_in_local_database
        [basenameBImrb,propBImrb] = crism_searchCDRmrb(propBI_search,'dwld',1);
        crism_get_dirpath_cdr(basenameBImrb,'dwld',dwld);
        BIdata = CRISMdata(basenameBImrb,'');
        for ii=1:length(EDRBIdataList_s)
            basenameEDRBI_i = EDRBIdataList_s(ii).basename;
            if isempty(extractMatchedBasename_v2(basenameEDRBI_i,BIdata.lbl.SOURCE_PRODUCT_ID))
                error('It seems no BIdata associated with %s.',basenameEDRBI_i);
            end
        end
    end

end

function [BKdata] = get_BKdata_fromDF(DFdata,binning_id_sp,dwld)
    propBK_search = crism_create_propCDR4basename;
    propBK_search.acro_calibration_type = 'BK';
    propBK_search.binning = crism_get_binning_id(DFdata.lbl.PIXEL_AVERAGING_WIDTH);
    propBK_search.wavelength_filter = DFdata.lbl.MRO_WAVELENGTH_FILTER;
    propBK_search.frame_rate = crism_get_frame_rate_id(DFdata.lbl.MRO_FRAME_RATE.value);
    propBK_search.exposure = DFdata.lbl.MRO_EXPOSURE_PARAMETER;
    propBK_search.sensor_id = DFdata.lbl.MRO_SENSOR_ID;
    [sclk_df_stop,p_df_stop] = DFdata.get_sclk_stop();
    propBK_search.sclk = floor(sclk_df_stop);
    propBK_search.partition = p_df_stop;
    propBK_search.binning = binning_id_sp;
    [basenameBKmrb,propBKmrb] = crism_searchCDRmrb(propBK_search,'dwld',1);
    crism_get_dirpath_cdr(basenameBKmrb,'dwld',dwld);
    BKdata = CRISMdata(basenameBKmrb,'');
    % check if the most recent before product is actually the one we are
    % looking for.
    if isempty(extractMatchedBasename_v2(DFdata.basename,BKdata.lbl.SOURCE_PRODUCT_ID))
        fprintf('It seems no BKdata associated with %s.\n in local database',DFdata.basename);
        fprintf('Connect to remote server...\n');
        [basenameBKmrb,propBKmrb] = crism_searchCDRmrb(propBK_search,'dwld',1);
        crism_get_dirpath_cdr(basenameBKmrb,'dwld',dwld);
        BKdata = CRISMdata(basenameBKmrb,'');
        if isempty(extractMatchedBasename_v2(DFdata.basename,BKdata.lbl.SOURCE_PRODUCT_ID))
            error('It seems no BKdata associated with %s.',DFdata.basename);
        end
    end

end