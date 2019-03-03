function [SPdata_o,RT14j_woc,RT14j,RT14h2_bk1_o,RT14h2_bk2_o] = minipipeline_calibration_IR_SP_wCDRSP_yuki(...
    SPdata,bkoption,varargin)
% [RT14j_woc,RT14j,RT14h2_bk1_o,RT14h2_bk2_o] = minipipeline_calibration_IR_SP_wCDRSP_yuki(...
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
%     'APBPRMVL'
%       option for a priori bad pixel removal {'HighOrd', 'None'}
%       'HighOrd': the image where a priori bad pixel removal is performed
%                  is used for the estimation of the higher order leaked 
%                  light
%       'None'   : the image where a priori bad pixel removal is NOT performed
%                  is used for the estimation of the higher order leaked 
%                  light
%       (default) 'HighOrd'
%   'DWLD','DOWNLOAD' : if download the data or not, 2: download, 1:
%                       access an only show the path, 0: nothing
%                       (default) 0
%   'OUT_FILE'       : path to the output file
%                       (default) ''
%   'Force'          : binary, whether or not to force performing
%                      pds_downloader. (default) false
%   ****** Parameters for manual BK production ****************************
%   'BK_DN4095_RMVL': binary, whether or not to perform replacement of saturated
%                  pixels or not.
%                  (default) true
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
bk_dn4095_rmvl = true;
bk_mean_robust = 1;
bk_bprmvl = false;
bk_mean_DN14 = true;
dwld = 0;
force = false;
outfile = '';

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
            case 'BK_DN4095_RMVL'
                bk_dn4095_rmvl = varargin{i+1};
            case 'BK_BPRMVL'
                bk_bprmvl = varargin{i+1};
            case 'BK_MEAN_ROBUST'
                bk_mean_robust = varargin{i+1};
            case 'BK_MEAN_DN14'
                bk_mean_DN14 = varargin{i+1};
            case {'DWLD','DOWNLOAD'}
                dwld = varargin{i+1};
            case 'FORCE'
                force = varargin{i+1};
            case 'OUT_FILE'
                outfile = varargin{i+1};
            otherwise
                % Hmmm, something wrong with the parameter string
                error(['Unrecognized option: ''' varargin{i} '''']);
        end
    end
end

if isempty(SPdata.basenamesCDR)
    SPdata.load_basenamesCDR('Download',dwld,'Force',force,'OUT_file',outfile); 
end
if isempty(SPdata.basenames_SOURCE_OBS)
    SPdata.load_basenames_SOURCE_OBS('Download',dwld,'Force',force,'OUT_file',outfile); 
end

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
    function [DFdata] = get_DFdata(EDRSPdata,obs_counter,dwld)
        propEDRDF = create_propOBSbasename();
        propEDRDF.obs_counter = obs_counter;
        propEDRDF.obs_class_type = EDRSPdata.prop.obs_class_type;
        propEDRDF.obs_id = EDRSPdata.prop.obs_id;
        propEDRDF.product_type = EDRSPdata.prop.product_type;
        propEDRDF.sensor_id = EDRSPdata.prop.sensor_id;
        propEDRDF.activity_id = 'DF';

        [~,~,~,~,~,basenameEDRDF] = get_dirpath_observation_fromProp(propEDRDF,'dwld',dwld);

        DFdata = CRISMdata(basenameEDRDF,'');
        DFdata.download(2);
    end
DFdata1 = get_DFdata(EDRSPdata,sprintf('%02s',dec2hex(hex2dec(EDRSPdata.prop.obs_counter)-1)),dwld);
DFdata1.download(dwld);
DFdata2 = get_DFdata(EDRSPdata,sprintf('%02s',dec2hex(hex2dec(EDRSPdata.prop.obs_counter)+1)),dwld);
DFdata2.download(dwld);

%-------------------------------------------------------------------------%
% get BPdata using sclk of the DFdata
    function [BPdata] = get_BPdata_fromDF(DFdata,dwld)
        propBP_search = create_propCDR4basename;
        propBP_search.acro_calibration_type = 'BP';
        propBP_search.binning = get_binning_id(DFdata.lbl.PIXEL_AVERAGING_WIDTH);
        propBP_search.wavelength_filter = DFdata.lbl.MRO_WAVELENGTH_FILTER;
        propBP_search.frame_rate = get_frame_rate_id(DFdata.lbl.MRO_FRAME_RATE{1});
        % for the BP products, exposure is almost always 0. Do not set it.
        % propBP_search.exposure = DFdata.lbl.MRO_EXPOSURE_PARAMETER;
        propBP_search.sensor_id = DFdata.lbl.MRO_SENSOR_ID;
        [sclk_df_stop,p_df_stop] = DFdata.get_sclk_stop();
        propBP_search.sclk = floor(sclk_df_stop);
        propBP_search.partition = p_df_stop;
        [basenameBPmrb,propBPmrb] = crism_searchCDRmrb(propBP_search,'dwld',1,'force',1);
        get_dirpath_cdr(basenameBPmrb,'dwld',dwld);
        BPdata = CRISMdata(basenameBPmrb,'');
        if isempty(extractMatchedBasename_v2(DFdata.basename,BPdata.lbl.SOURCE_PRODUCT_ID))
            error('It seems no BPdata associated with %s.',DFdata.basename);
        end
    end
BPdata1 = get_BPdata_fromDF(DFdata1,dwld);
BPdata2 = get_BPdata_fromDF(DFdata2,dwld);

%-------------------------------------------------------------------------%
% get BIdata using source_obs of SPdata
EDRBIdataList = SPdata.read_SOURCE_OBS('BI');

% select EDR BI data with same frame rate.
frame_rate = SPdata.lbl.MRO_FRAME_RATE{1};
EDRBIdataList_s = CRISMdata.empty(0,0);
for i=1:length(EDRBIdataList)
    if EDRBIdataList(i).lbl.MRO_FRAME_RATE{1} == frame_rate
        EDRBIdataList_s = [EDRBIdataList_s EDRBIdataList(i)];
    end
end

    function [BIdata] = get_BIdata_fromEDRBI(EDRBIdataList_s,dwld)
        % get sclk
        sclk_EDRBI_list = zeros(1,length(EDRBIdataList_s));
        for ii=1:length(EDRBIdataList_s)
            [sclk_stop_i] = EDRBIdataList_s(ii).get_sclk_stop();
            sclk_EDRBI_list(ii) = sclk_stop_i;
        end

        sclk_stop = max(sclk_EDRBI_list);

        propBI_search = create_propCDR4basename;
        propBI_search.acro_calibration_type = 'BI';
        propBI_search.binning = get_binning_id(EDRBIdataList_s(1).lbl.PIXEL_AVERAGING_WIDTH);
        propBI_search.wavelength_filter = EDRBIdataList_s(1).lbl.MRO_WAVELENGTH_FILTER;
        propBI_search.frame_rate = get_frame_rate_id(EDRBIdataList_s(1).lbl.MRO_FRAME_RATE{1});
        propBI_search.sensor_id = EDRBIdataList_s(1).lbl.MRO_SENSOR_ID;
        propBI_search.sclk = floor(sclk_stop);
        [basenameBImrb,propBImrb] = crism_searchCDRmrb(propBI_search,'dwld',1,'force',1);
        get_dirpath_cdr(basenameBImrb,'dwld',dwld);
        BIdata = CRISMdata(basenameBImrb,'');
        for ii=1:length(EDRBIdataList_s)
            basenameEDRBI_i = EDRBIdataList_s(ii).basename;
            if isempty(extractMatchedBasename_v2(basenameEDRBI_i,BIdata.lbl.SOURCE_PRODUCT_ID))
                error('It seems no BIdata associated with %s.',basenameEDRBI_i);
            end
        end
    end

[BIdata] = get_BIdata_fromEDRBI(EDRBIdataList_s);

%-------------------------------------------------------------------------%
% get BKdata using sclk of the DFdata
    function [BKdata] = get_BKdata_fromDF(DFdata,dwld)
        propBK_search = create_propCDR4basename;
        propBK_search.acro_calibration_type = 'BK';
        propBK_search.binning = get_binning_id(DFdata.lbl.PIXEL_AVERAGING_WIDTH);
        propBK_search.wavelength_filter = DFdata.lbl.MRO_WAVELENGTH_FILTER;
        propBK_search.frame_rate = get_frame_rate_id(DFdata.lbl.MRO_FRAME_RATE{1});
        propBK_search.exposure = DFdata.lbl.MRO_EXPOSURE_PARAMETER;
        propBK_search.sensor_id = DFdata.lbl.MRO_SENSOR_ID;
        [sclk_df_stop,p_df_stop] = DFdata.get_sclk_stop();
        propBK_search.sclk = floor(sclk_df_stop);
        propBK_search.partition = p_df_stop;
        [basenameBKmrb,propBKmrb] = crism_searchCDRmrb(propBK_search,'dwld',1,'force',1);
        get_dirpath_cdr(basenameBKmrb,'dwld',dwld);
        BKdata = CRISMdata(basenameBKmrb,'');
        % check if the most recent before product is actually the one we are
        % looking for.
        if isempty(extractMatchedBasename_v2(DFdata.basename,BKdata.lbl.SOURCE_PRODUCT_ID))
            error('It seems no BKdata associated with %s.',DFdata.basename);
        end
    end
BKdata1 = get_BKdata_fromDF(DFdata1,dwld);
BKdata2 = get_BKdata_fromDF(DFdata2,dwld);

%-------------------------------------------------------------------------%
% main pipeline
[SPdata_o,RT14j_woc,RT14j,RT14h2_bk1_o,RT14h2_bk2_o] = minipipeline_calibration_IR_SP_yuki(...
    EDRSPdata,DFdata1,DFdata2,BKdata1,BKdata2,BPdata1,BPdata2,BIdata,...
    PPdata,BSdata,DBdata,EBdata,HDdata,HKdata,GHdata,DMdata,LCdata,LLdata,...
    bkoption,'SAVE_MEMORY',save_mem,'APBPRMVL',apbprmvl,...
    'BK_DN4095_RMVL',bk_dn4095_rmvl,'BK_BPRMVL',bk_bprmvl,...
    'BK_MEAN_ROBUST',bk_mean_robust,'BK_MEAN_DN14',bk_mean_DN14);

end