function [SPdata_o,RT14j_woc_mod,RT14j_mod] = minipipeline_calibration_VNIR_SP_wCDRSP_yuki(...
    SPdata,TRRIFdata,varargin)
save_mem = false;
apbprmvl = 'HighOrd';
saturation_rmvl = 2;
bk_saturation_rmvl = 2;
bk_mean_robust = 1;
bk_bprmvl = false;
bk_mean_DN14 = true;
dwld = 0;
force = false;
outfile = '';
BIdata = [];
mean_DN14 = false;

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
            case 'FORCE'
                force = varargin{i+1};
            case 'OUT_FILE'
                outfile = varargin{i+1};
            case 'BIDATA'
                BIdata = varargin{i+1};
            otherwise
                % Hmmm, something wrong with the parameter string
                error(['Unrecognized option: ''' varargin{i} '''']);
        end
    end
end

% if isempty(SPdata.basenamesCDR)
    SPdata.load_basenamesCDR('Download',dwld,'Force',force,'OUT_file',outfile); 
% end
% if isempty(SPdata.basenames_SOURCE_OBS)
    SPdata.load_basenames_SOURCE_OBS('Download',dwld,'Force',force,'OUT_file',outfile); 
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
DFdata2 = SPdata.source_obs.DF(2);
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
PPdata = TRRIFdata.readCDR('PP');
DBdata = TRRIFdata.readCDR('DB');
EBdata = TRRIFdata.readCDR('EB');
HDdata = TRRIFdata.readCDR('HD');
HKdata = TRRIFdata.readCDR('HK');
GHdata = TRRIFdata.readCDR('GH');
LCdata = TRRIFdata.readCDR('LC');
DMdata = TRRIFdata.readCDR('DM');
VLdata = TRRIFdata.readCDR('VL');

%-------------------------------------------------------------------------%
% main pipeline
[SPdata_o,RT14j_woc,RT14j] = minipipeline_calibration_VNIR_SP_yuki( EDRSPdata,DFdata1,DFdata2,...
    PPdata,DBdata,EBdata,HDdata,HKdata,GHdata,VLdata,DMdata,LCdata);

% [BP1nan] = formatBP1nan(BPdata1);
% BPpri1nan = formatBPpri1nan(BPdata1,BPdata2);
%[BIdata_o,imgBI] = minipipeline_calibration_IR_BI_wCDRBI_yuki(BIdata,'DN4095_RMVL',0,'BPRMVL',0,'MEAN_ROBUST',1);

SSdata = TRRIFdata.readCDR('SS');
SHdata = TRRIFdata.readCDR('SH');
SPdata_o = CRISMdata(SPdata.basename,'');
SPdata_o.img = RT14j_woc;
[MP] = calculate_MP(SPdata_o,SSdata,SHdata);

SHdata.readimg();
SC = SHdata.img(2,:,:);

RT14j_mod = RT14j ./ (1 + MP.* SC);
RT14j_woc_mod = RT14j_woc ./ (1 + MP.* SC);

SPdata_o.img = RT14j_woc;

end
