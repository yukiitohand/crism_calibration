function [RDn,RDn_woc,RDn_bk1_o,RDn_bk2_o] = crmcal_pipeline_IR_IF_wTRRIF_yuki(...
    TRRIFdata,TRRIFdataVNIR,crism_obs,bkoption,varargin)
% [RDn,RDn_woc,RDn_bk1_o,RDn_bk2_o] = crmcal_pipeline_IR_IF_wTRRIF_yuki(...
%     TRRIFdata,TRRIFdataVNIR,crism_obs,bkoption,varargin)
%   Main pipeline for the calibration of the CRISM images using TRR IF
%   data. This is a wrapper function for 
%    "crmcal_pipeline_IR_IF_yuki"
%
%  Input Parameters
%   TRRIFdata: TRR IF, CRISMdata obj
%   TRRIFdataVNIR: TRR IF, CRISMdata obj, corresponding VNIR data
%   bkoption: option for the use of background image. {1,2}
%             1: linear background estimation using prior and post dark
%                measurements.
%             2: flat background estimation using only prior dark
%                measurements.
%   OUTPUTS
%    RDn:       [L x S x B] radiance image
%    RDn_woc: [L x S x B] radiance image, non-interpolated version
%    RDn_bk1_o: processed prior dark measurement
%    RDn_bk2_o: processed post dark measurement
%   Optional Parameters
%    'SAVE_MEMORY'
%       saving memory or not. true or false
%       (default) true
%    'DWLD','DOWNLOAD' : if download the data or not, 2: download, 1:
%                       access an only show the path, 0: nothing
%                       (default) 0
%   'Force_dwld'     : binary, whether or not to force performing
%                      pds_downloader. (default) false
%   'DWLD_INDEX_CACHE_UPDATE' : boolean, whether or not to update index.html 
%        (default) false
%   'VERBOSE_DWLD'   : boolean, whether or not to show the downloading
%                      operations.
%                      (default) true
%   'DWLD_OVERWRITE' : if overwrite the file if exists
%                      (default) 0
%    'MODE_SP'
%       {'MAN','SOC'}
%       'SOC': CDR SP data produced by science operation center (SOC) is
%              used, MP=0 in this case.
%       'MAN': SP and MP are manually produced in this toolbox.
%       (default) 'MAN'
%    'MEAN_ROBUST' : integer {0,1}, mode for how mean operation is performed.
%        0: DN14e_df = nanmean(DN14d_df(:,:,:),1);
%        1: DN14e_df = robust_v2('mean',DN14d_df,1,'NOutliers',2);
%      (default) 1
%    'APBPRMVL'
%       option for a priori bad pixel removal {'HighOrd', 'None'}
%       'HighOrd': the image where a priori bad pixel removal is performed
%                  is used for the estimation of the higher order leaked 
%                  light
%       'None'   : the image where a priori bad pixel removal is NOT performed
%                  is used for the estimation of the higher order leaked 
%                  light
%       (default) 'HighOrd'
%    'SATURATiON_RMVL': integer, how to perform replacement of saturated
%           pixles {0,1,2}
%           0: no removal
%           1: digital saturation is removed
%           2: analogue saturation is also removed
%           (default) 2
%    'FLAT_FIELD' : 
%      boolean, whether or not to perform flat field correction using 
%      NUdata or not.
%      (default) true
%    'MEAN_DN14'  : binary,when mean operation is performed
%                  1: before non-linearity correction
%                  0: last (after divided by integration time
%                  (default) 1
%    'SP_APBPRMVL'
%       option for a priori bad pixel removal {'HighOrd', 'None'}
%       'HighOrd': the image where a priori bad pixel removal is performed
%                  is used for the estimation of the higher order leaked 
%                  light
%       'None'   : the image where a priori bad pixel removal is NOT performed
%                  is used for the estimation of the higher order leaked 
%                  light
%       (default) 'HighOrd'
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
%
%   'SPdata'
%       CRISMdata obj, CDR SPdata
%       (default) []
%   'BIdata'
%       CRISMdata obj, CDR BIdata
%       (default) []

save_mem = false;
dwld = 0;
force_dwld = false;
dwld_index_cache_update = 0;
verbose_dwld = 1;
dwld_overwrite = 0;
mode_SP = 'MAN';
SPdata  = [];
SPdataVNIR = [];
%
apbprmvl = 'HighOrd';
saturation_rmvl = 2;
flat_field = true;
%
sp_apbprmvl = 'HighOrd';
sp_saturation_rmvl = 2;
sp_mean_DN14 = true;
sp_mean_robust = true;
%
bk_saturation_rmvl = 2;
bk_mean_robust = 1;
bk_bprmvl = false;
bk_mean_DN14 = true;
%
BIdata = [];


if (rem(length(varargin),2)==1)
    error('Optional parameters should always go by pairs');
else
    for i=1:2:(length(varargin)-1)
        switch upper(varargin{i})
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % for some top level option
            case 'SAVE_MEMORY'
                save_mem = varargin{i+1};
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
            case 'MODE_SP'
                mode_SP = varargin{i+1};
                if ~any(strcmpi(mode_SP,{'MAN','SOC'}))
                    error('apbprmvl (%s) should be either {"MAN","SOC"}',apbprmvl);
                end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % for main processing
            case 'APBPRMVL'
                apbprmvl = varargin{i+1};
                if ~any(strcmpi(apbprmvl,{'HighOrd','None'}))
                    error('apbprmvl (%s) should be either {"HighOrd","None"}',apbprmvl);
                end
            case 'SATURATION_RMVL'
                saturation_rmvl = varargin{i+1};
            case 'FLAT_FIELD'
                flat_field = varargin{i+1};
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % for SP processing 
            case 'SPDATA'
                SPdata = varargin{i+1};
            case 'SPDATAVNIR'
                SPdataVNIR = varargin{i+1};
            case 'SP_APBPRMVL'
                sp_apbprmvl = varargin{i+1};
                if ~any(strcmpi(apbprmvl,{'HighOrd','None'}))
                    error('sp apbprmvl (%s) should be either {"HighOrd","None"}',sp_apbprmvl);
                end
            case 'SP_MEAN_DN14'
                sp_mean_DN14 = varargin{i+1};
            case 'SP_MEAN_ROBUST'
                sp_mean_robust = varargin{i+1};
            case 'SP_SATURATION_RMVL'
                sp_saturation_rmvl = varargin{i+1};
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % for BK processing (also for BK of SP)
            case 'BK_SATURATION_RMVL'
                bk_saturation_rmvl = varargin{i+1};
            case 'BK_BPRMVL'
                bk_bprmvl = varargin{i+1};
            case 'BK_MEAN_ROBUST'
                bk_mean_robust = varargin{i+1};
            case 'BK_MEAN_DN14'
                bk_mean_DN14 = varargin{i+1};
            
                
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            case 'BIDATA'
                BIdata = varargin{i+1};
                
            otherwise
                error('Undefined keyword: %s',varargin{i});
        end
    end
end

if (isempty(SPdata) && ~isempty(SPdataVNIR)) || (~isempty(SPdata) && isempty(SPdataVNIR))
    error('You only give SPdata or SPdataVNIR. Enter both if you do custom SP');
end

TRRIFdata.load_basenamesCDR('Download',dwld,...
    'OVERWRITE',dwld_overwrite,'INDEX_CACHE_UPDATE',dwld_index_cache_update);
TRRIFdata.load_basenames_SOURCE_OBS('Download',dwld,...
    'OVERWRITE',dwld_overwrite,'INDEX_CACHE_UPDATE',dwld_index_cache_update);

%-------------------------------------------------------------------------%
% get EDRdata from TRRIFdata
basenameEDR = TRRIFdata.basenames_SOURCE_OBS.SC;
if iscell(basenameEDR)
    error('Mulitple EDR is used. not supported.');
end
EDRdata = CRISMdata(basenameEDR,'');
EDRdata.download(dwld);

%-------------------------------------------------------------------------%
% get DFdata
[DFdata1,DFdata2] = crism_get_DFdata4SC(EDRdata,crism_obs);
if isempty(DFdata2)
    DFdata2 = DFdata1;
end

is_same_df1_df2 = strcmpi(DFdata1.basename,DFdata2.basename);

%-------------------------------------------------------------------------%
% get BIdata 
if isempty(BIdata)
    TRRIFdata.readCDR('BI'); BIdata = [];
    for i=1:length(TRRIFdata.cdr.BI)
        bidata = TRRIFdata.cdr.BI(i);
        if bidata.lbl.MRO_FRAME_RATE.value == TRRIFdata.lbl.MRO_FRAME_RATE.value
            if bidata.lbl.PIXEL_AVERAGING_WIDTH == TRRIFdata.lbl.PIXEL_AVERAGING_WIDTH
                BIdata = [BIdata bidata];
            end
        end
    end
end
if length(BIdata)>1
    error('BIdata detected twice');
end

%-------------------------------------------------------------------------%
% get BKdata using DFdata
obs_id_scene = TRRIFdata.get_obsid;
binning_id_scene = crism_get_binning_id(TRRIFdata.lbl.PIXEL_AVERAGING_WIDTH);
wvfilter_id_scene = TRRIFdata.lbl.MRO_WAVELENGTH_FILTER;
TRRIFdata.readCDR('BK');
BKdata1 = []; BKdata2 = [];
for j=1:length(TRRIFdata.cdr.BK)
    bkdata = TRRIFdata.cdr.BK(j);
    if strcmpi(bkdata.get_obsid, obs_id_scene)
        if strcmpi(bkdata.get_obs_number,DFdata1.get_obs_number)
            BKdata1 = [BKdata1 bkdata];
        elseif ~isempty(DFdata2) && strcmpi(bkdata.get_obs_number,DFdata2.get_obs_number)
            BKdata2 = [BKdata2 bkdata];
        else % Exceptions: sometimes BKdata is produced by different DF measurements.
            if bkdata.prop.binning == binning_id_scene && bkdata.prop.wavelength_filter == wvfilter_id_scene
                if hex2dec(bkdata.get_obs_number()) < hex2dec(TRRIFdata.get_obs_number())
                    BKdata1 = [BKdata1 bkdata];
                    % load different DFdata1
                    bkdata.load_basenames_SOURCE_OBS('dwld',2);
                    DFdata1 = CRISMdata(bkdata.basenames_SOURCE_OBS.DF, '');
                else
                    BKdata2 = [BKdata2 bkdata];
                    bkdata.load_basenames_SOURCE_OBS('dwld',2);
                    DFdata2 = CRISMdata(bkdata.basenames_SOURCE_OBS.DF, '');
                end
            end

        end
    end
end
if is_same_df1_df2
    BKdata2 = BKdata1;
end

if length(BKdata1)>1 || length(BKdata2)>1
    error('BKdata detected twice');
end

%-------------------------------------------------------------------------%
% get BPdata
[BPdata1,BPdata2,BPdata_post] = crism_load_BPdataSC_fromDF(TRRIFdata,DFdata1.basename,DFdata2.basename);
if length(BPdata1)>1 || length(BPdata2)>1
    error('BPdata detected twice');
end

%-------------------------------------------------------------------------%
% Read other CDRs
PPdata = TRRIFdata.readCDR('PP');
BSdata = TRRIFdata.readCDR('BS'); DBdata = TRRIFdata.readCDR('DB');
EBdata = TRRIFdata.readCDR('EB'); HDdata = TRRIFdata.readCDR('HD');
HKdata = TRRIFdata.readCDR('HK');
GHdata = TRRIFdata.readCDR('GH');
LCdata = TRRIFdata.readCDR('LC');
DMdata = TRRIFdata.readCDR('DM');
LLdata = TRRIFdata.readCDR('LL');
VLdata = TRRIFdata.readCDR('VL');
%
TRRIFdata.readCDR('SP');
for i=1:length(TRRIFdata.cdr.SP)
    spdata = TRRIFdata.cdr.SP(i);
    spdata_prop = crism_getProp_basenameCDR4(spdata.basename);
    switch upper(spdata_prop.sensor_id)
        case 'L'
            if isempty(SPdata)
                SPdata = spdata;
            end
        case 'S'
            if isempty(SPdataVNIR)
                SPdataVNIR = spdata;
            end
        otherwise
            error('sensor_id %s is wrong',sensor_id);
    end
end
%
SSdata = TRRIFdata.readCDR('SS');
SHdata = TRRIFdata.readCDR('SH');
NUdata = TRRIFdata.readCDR('NU');

%-------------------------------------------------------------------------%
% Process Background
[RT14g_bkgd,BKdata1_o,RT14g_df1] = crmcal_minipipeline_IR_BK_wCDRBK_yuki(...
    BKdata1,BIdata,BPdata1,'DWLD',dwld,'FORCE_DWLD',force_dwld, ...
    'DWLD_OVERWRITE',dwld_overwrite,'VERBOSE_DWLD',verbose_dwld, ...
    'DWLD_INDEX_CACHE_UPDATE',dwld_index_cache_update,'VERBOSE_DWLD',verbose_dwld,...
    'SAVE_MEMORY',save_mem,...
    'SATURATION_RMVL',bk_saturation_rmvl,'BPRMVL',bk_bprmvl,...
    'MEAN_ROBUST',bk_mean_robust,'MEAN_DN14',bk_mean_DN14);
if ~isempty(BKdata2)
    [RT14g_bkgd,BKdata2_o,RT14g_df2] = crmcal_minipipeline_IR_BK_wCDRBK_yuki(...
        BKdata2,BIdata,BPdata2,'DWLD',dwld,'FORCE_DWLD',force_dwld,...
        'DWLD_OVERWRITE',dwld_overwrite,'VERBOSE_DWLD',verbose_dwld, ...
        'DWLD_INDEX_CACHE_UPDATE',dwld_index_cache_update,...
        'SAVE_MEMORY',save_mem,...
        'SATURATION_RMVL',bk_saturation_rmvl,'BPRMVL',bk_bprmvl,...
        'MEAN_ROBUST',bk_mean_robust,'MEAN_DN14',bk_mean_DN14);
else
    RT14g_bkgd = []; BKdata2_o = []; RT14g_df2 = [];
end
% [~,BKdata1_o,RT14g_df1] = crmcal_minipipeline_IR_BK_yuki(...
%     DFdata1,PPdata,BSdata,DBdata,EBdata,HDdata,HKdata,BIdata,DMdata,...
%     BPdata1,GHdata,VLdata,LCdata,'SATURATION_RMVL',bk_saturation_rmvl,'BPRMVL',bk_bprmvl,...
%     'MEAN_ROBUST',bk_mean_robust,'MEAN_DN14',bk_meanDN14);
% [~,BKdata2_o,RT14g_df2] = crmcal_minipipeline_IR_BK_yuki(...
%     DFdata2,PPdata,BSdata,DBdata,EBdata,HDdata,HKdata,BIdata,DMdata,...
%     BPdata2,GHdata,VLdata,LCdata,'SATURATION_RMVL',bk_saturation_rmvl,'BPRMVL',bk_bprmvl,...
%     'MEAN_ROBUST',bk_mean_robust,'MEAN_DN14',bk_meanDN14);

%-------------------------------------------------------------------------%
% main pipeline
switch upper(mode_SP)
    case 'MAN'
        %-----------------------------------------------------------------%
        % calculate MP
        [SPdataVNIR_o,RT14jVNIR_woc_mod,RT14j,RT14jVNIR_mod,MP] = crmcal_minipipeline_VNIR_SP_wCDRSP_yuki(...
            SPdataVNIR,TRRIFdataVNIR,'DWLD',dwld,'FORCE_DWLD',force_dwld, ...
            'DWLD_OVERWRITE',dwld_overwrite,'VERBOSE_DWLD',verbose_dwld, ...
            'DWLD_INDEX_CACHE_UPDATE',dwld_index_cache_update);

        % [SPdataMP,SSdataMP,SHdataMP] = crmcal_selectCDR4MP(SPdataVNIR);
        % [MP] = crmcal_calculate_MP(SPdataVNIR_o,SSdataMP,SHdataMP);

        %-----------------------------------------------------------------%
        % process SP (MP IS NOT APPLIED TO SPDATA_o)
        [SPdata_o,RT14j_woc,RT14j,RT14h2_bk1_o,RT14h2_bk2_o,BPdata1_sp,BPdata2_sp]...
                  = crmcal_minipipeline_IR_SP_wCDRSP_yuki(...
                    SPdata,bkoption,'DWLD',dwld,'FORCE_DWLD',force_dwld, ...
                    'DWLD_OVERWRITE',dwld_overwrite,'VERBOSE_DWLD',verbose_dwld, ...
                    'DWLD_INDEX_CACHE_UPDATE',dwld_index_cache_update,...
                    'SAVE_MEMORY',save_mem,'APBPRMVL',sp_apbprmvl,...
                    'MEAN_DN14',sp_mean_DN14,'SATURATION_RMVL',sp_saturation_rmvl,'MEAN_ROBUST',sp_mean_robust,...
                    'BK_SATURATION_RMVL',bk_saturation_rmvl,'BK_BPRMVL',bk_bprmvl,...
                    'BK_MEAN_ROBUST',bk_mean_robust,'BK_MEAN_DN14',bk_mean_DN14);
        % Note
        %  BK_* are assumed to be same as the ones for BK processing.
        [RDn,RDn_woc,RDn_bk1_o,RDn_bk2_o] = crmcal_pipeline_IR_IF_yuki(...
            EDRdata,DFdata1,DFdata2,BKdata1,BKdata2,BPdata1,BPdata2,BIdata,...
            PPdata,BSdata,DBdata,EBdata,HDdata,HKdata,GHdata,VLdata,DMdata,LCdata,...
            LLdata,SPdata,SSdata,SHdata,NUdata,...
            SPdata_o,BKdata1_o,BKdata2_o,RT14g_df1,RT14g_df2,bkoption,MP,...
            'Saturation_rmvl',saturation_rmvl,'apbprmvl',apbprmvl,'SAVE_MEMORY',save_mem,'FLAT_FIELD',flat_field);
    case 'SOC'
        MP = 0;
        SPdata_o = [];
        [RDn,RDn_woc,RDn_bk1_o,RDn_bk2_o] = crmcal_pipeline_IR_IF_yuki(...
            EDRdata,DFdata1,DFdata2,BKdata1,BKdata2,BPdata1,BPdata2,BIdata,...
            PPdata,BSdata,DBdata,EBdata,HDdata,HKdata,GHdata,VLdata,DMdata,LCdata,...
            LLdata,SPdata,SSdata,SHdata,NUdata,...
            SPdata_o,BKdata1_o,BKdata2_o,RT14g_df1,RT14g_df2,bkoption,MP,...
            'Saturation_rmvl',saturation_rmvl,'apbprmvl',apbprmvl,'SAVE_MEMORY',save_mem);
    otherwise
        error('Undefined MODE_SP=%s',mode_SP);
end
% [BP1nan] = crism_formatBP1nan(BPdata1);
% BPpri1nan = crism_formatBPpri1nan(BPdata1,BPdata2);
%[BIdata_o,imgBI] = crmcal_minipipeline_IR_BI_wCDRBI_yuki(BIdata,'DN4095_RMVL',0,'BPRMVL',0,'MEAN_ROBUST',1);
end
