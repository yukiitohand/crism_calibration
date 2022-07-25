function [SPdata_o,RT14j_woc_bn_mod,RT14j_bn,RT14j_bn_mod,MP] = minipipeline_calibration_VNIR_SP_wCDRSP_yuki(...
    SPdata,TRRIFdata,varargin)
%   Main pipeline for the calibration of the CRISM images using CDR SPdata,
%   VNIR. This is a wrapper function for 
%    "crmcal_minipipeline_VNIR_SP_yuki"
%   Binning is also applied.
%  Input Parameters
%   SPdata: CDR SPdata, CRISMdata obj
%   TRRIFdata: TRRIF VNIR data corresponding to SPdata
%   bkoption: option for the use of background image. {1,2}
%             1: linear background estimation using prior and post dark
%                measurements.
%             2: flat background estimation using only prior dark
%                measurements.
%   OUTPUTS
%    SPdata_o: SPdata that stores processed image at img. RT14j_woc is stored.
%    RT14j_woc_bn_mod: non interpolated version of the image of SPdata applied
%                  MP parameters
%    RT14j_bn: interpolated version of SPdata (image)
%    RT14j_mod: MP corrected RT14j
%    MP: scalar, shutter mirror parameter
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
%   'MEAN_ROBUST': integer {0,1}, mode for how mean operation is performed.
%        0: DN14e_df = nanmean(DN14d_df(:,:,:),1);
%        1: DN14e_df = robust_v2('mean',DN14d_df,1,'NOutliers',2);
%      (default) 1  
%   'SATURATION_RMVL': integer, how to perform replacement of saturated
%           pixles {0,1,2}
%           0: no removal
%           1: digital saturation is removed
%           2: analogue saturation is also removed
%           (default) 2
%   'BK_MEAN_ROBUST': integer {0,1}, mode for how mean operation is performed.
%        0: DN14e_df = nanmean(DN14d_df(:,:,:),1);
%        1: DN14e_df = robust_v2('mean',DN14d_df,1,'NOutliers',2);
%      (default) 1  
%   ****** Parameters for manual BK production ****************************
%   'BK_MEAN_ROBUST': integer {0,1}, mode for how mean operation is performed.
%        0: DN14e_df = nanmean(DN14d_df(:,:,:),1);
%        1: DN14e_df = robust_v2('mean',DN14d_df,1,'NOutliers',2);
%      (default) 1  
%
% for future release
%     'APBPRMVL'
%       option for a priori bad pixel removal {'HighOrd', 'None'}
%       'HighOrd': the image where a priori bad pixel removal is performed
%                  is used for the estimation of the higher order leaked 
%                  light
%       'None'   : the image where a priori bad pixel removal is NOT performed
%                  is used for the estimation of the higher order leaked 
%                  light
%       (default) 'None'
%   'MEAN_DN14'  : binary,when mean operation is performed
%                  1: before non-linearity correction
%                  0: last (after divided by integration time
%      (default) 1
save_mem = false;
% apbprmvl = 'None';
saturation_rmvl = 2;
mean_robust = 1;
bk_mean_robust = 1;
dwld = 0;
force_dwld = false;
dwld_index_cache_update = 0;
verbose_dwld = 1;
dwld_overwrite = 0;
% BIdata = [];
% mean_DN14 = true;

if (rem(length(varargin),2)==1)
    error('Optional parameters should always go by pairs');
else
    for i=1:2:(length(varargin)-1)
        switch upper(varargin{i})
            case 'SAVE_MEMORY'
                save_mem = varargin{i+1};
            % case 'APBPRMVL'
            %     apbprmvl = varargin{i+1};
            %     if ~any(strcmpi(apbprmvl,{'HighOrd','None'}))
            %         error('apbprmvl (%s) should be either {"HighOrd","None"}',apbprmvl);
            %     end
            % case 'MEAN_DN14'
            %    mean_DN14 = varargin{i+1};
            case 'SATURATION_RMVL'
                saturation_rmvl = varargin{i+1};
            case 'BK_MEAN_ROBUST'
                bk_mean_robust = varargin{i+1};
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
            % case 'BIDATA'
            %    BIdata = varargin{i+1};
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
    SPdata.load_basenames_SOURCE_OBS('Download',dwld,...
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
SPdata.read_SOURCE_OBS('DF');
DFdata1 = SPdata.source_obs.DF(1);
if length(SPdata.source_obs.DF)==1
    DFdata2 = DFdata1;
else
    DFdata2 = SPdata.source_obs.DF(2);
end
DFdata1.download(dwld);
DFdata2.download(dwld);

%-------------------------------------------------------------------------%
% get BIdata using sclk of the DFdata
% BIdata1 = get_BIdata_fromDF(DFdata1,dwld);
% BIdata2 = get_BIdata_fromDF(DFdata2,dwld);

%-------------------------------------------------------------------------%
% get BPdata using sclk of the DFdata
% BPdata1 = get_BPdata_fromDF(DFdata1,dwld);
% BPdata2 = get_BPdata_fromDF(DFdata2,dwld);

%-------------------------------------------------------------------------%
% Read other CDRs
% PPdata = TRRIFdata.readCDR('PP');
% DBdata = TRRIFdata.readCDR('DB');
% EBdata = TRRIFdata.readCDR('EB');
% HDdata = TRRIFdata.readCDR('HD');
% HKdata = TRRIFdata.readCDR('HK');
% GHdata = TRRIFdata.readCDR('GH');
% LCdata = TRRIFdata.readCDR('LC');
% DMdata = TRRIFdata.readCDR('DM');
% VLdata = TRRIFdata.readCDR('VL');
PPdata = crism_searchCDR6mrb('PP',SPdata.prop.sclk,'sensor_id','S');
DBdata = crism_searchCDR6mrb('DB',SPdata.prop.sclk,'sensor_id','S');
EBdata = crism_searchCDR6mrb('EB',SPdata.prop.sclk,'sensor_id','S');
HDdata = crism_searchCDR6mrb('HD',SPdata.prop.sclk,'sensor_id','J');
HKdata = crism_searchCDR6mrb('HK',SPdata.prop.sclk,'sensor_id','J');
GHdata = crism_searchCDR6mrb('GH',SPdata.prop.sclk,'sensor_id','S');
LCdata = crism_searchCDR6mrb('LC',SPdata.prop.sclk,'sensor_id','S');
VLdata = crism_searchCDR6mrb('VL',SPdata.prop.sclk,'sensor_id','S');

DMdata = crism_searchCDR4mrb('DM',SPdata.prop.sclk,'sensor_id','S',...
    'wavelength_filter',EDRSPdata.lbl.MRO_WAVELENGTH_FILTER,'binning',crism_get_binning_id(EDRSPdata.lbl.PIXEL_AVERAGING_WIDTH));
DMdata_SP = crism_searchCDR4mrb('DM',SPdata.prop.sclk,'sensor_id','S',...
    'wavelength_filter',SPdata.prop.wavelength_filter,'binning',SPdata.prop.binning);
%-------------------------------------------------------------------------%
% main pipeline
[~,RT14j_woc,RT14j] = crmcal_minipipeline_VNIR_SP_yuki( EDRSPdata,DFdata1,DFdata2,...
    PPdata,DBdata,EBdata,HDdata,HKdata,GHdata,VLdata,DMdata_SP,LCdata,...
    'save_memory',save_mem,'saturation_rmvl',saturation_rmvl,...
    'bk_mean_robust',bk_mean_robust,'mean_robust',mean_robust,...
    'BINNING_SP',binning_id_sp);

% [BP1nan] = crism_formatBP1nan(BPdata1);
% BPpri1nan = crism_formatBPpri1nan(BPdata1,BPdata2);
%[BIdata_o,imgBI] = crmcal_minipipeline_IR_BI_wCDRBI_yuki(BIdata,'DN4095_RMVL',0,'BPRMVL',0,'MEAN_ROBUST',1);

%-------------------------------------------------------------------------%
% binning
RT14j_woc_bn = crism_bin_image_frames(RT14j_woc,'binning',SPdata.prop.binning);
RT14j_bn = crism_bin_image_frames(RT14j,'binning',SPdata.prop.binning);

RT14j_woc_bn = crmcal_apply_DM(RT14j_woc_bn,DMdata_SP);
RT14j_bn     = crmcal_apply_DM(RT14j_bn,DMdata_SP);

%-------------------------------------------------------------------------%
% post processing
SSdata = crism_searchCDR4mrb('SS',SPdata.prop.sclk,'sensor_id','S',...
    'wavelength_filter',SPdata.prop.wavelength_filter,'binning',SPdata.prop.binning,...
    'version',2);
SHdata = crism_searchCDR4mrb('SH',SPdata.prop.sclk,'sensor_id','S',...
    'wavelength_filter',SPdata.prop.wavelength_filter,'binning',SPdata.prop.binning);
% SSdata = TRRIFdata.readCDR('SS');
% SHdata = TRRIFdata.readCDR('SH');
SPdata_o = CRISMdata(SPdata.basename,'');
SPdata_o.img = RT14j_woc_bn;
[MP] = crmcal_calculate_MP(SPdata_o,SSdata,SHdata);

SHdata.readimg();
SC = SHdata.img(2,:,:);

RT14j_bn_mod = RT14j_bn ./ (1 + MP.* SC);
RT14j_woc_bn_mod = RT14j_woc_bn ./ (1 + MP.* SC);

SPdata_o.img = RT14j_woc_bn;

end

