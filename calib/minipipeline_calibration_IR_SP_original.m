function [SPdata_o,RT14j_woc_mod,RT14j_mod,RT14h2_bk1_o_mod,RT14h2_bk2_o_mod] = minipipeline_calibration_IR_SP_original(...
    SPdata,SPdataVNIR,bkoption,varargin)
% [RT14j_woc,RT14j,RT14h2_bk1_o,RT14h2_bk2_o] = minipipeline_calibration_IR_SP_original(...
%     SPdata,bkoption,varargin)
%   Mini pipeline for the calibration of the CRISM images using SP CDR
%   data, and applying a shutter mirror parameter.
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
%   'DWLD','DOWNLOAD' : if download the data or not, 2: download, 1:
%                       access an only show the path, 0: nothing
%                       (default) 0
%   'OUT_FILE'       : path to the output file
%                       (default) ''
%   'Force'          : binary, whether or not to force performing
%                      pds_downloader. (default) false
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
bk_saturation_rmvl = 2;
bk_mean_robust = 1;
bk_bprmvl = false;
bk_mean_DN14 = true;
dwld = 0;
force = false;
outfile = '';
BIdata = [];
mean_DN14 = true;
mean_robust = true;

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
            case 'MEAN_ROBUST'
                mean_robust = varargin{i+1};
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

if isempty(SPdata.basenamesCDR)
    SPdata.load_basenamesCDR('Download',dwld,'Force',force,'OUT_file',outfile); 
end
if isempty(SPdata.basenames_SOURCE_OBS)
    SPdata.load_basenames_SOURCE_OBS('Download',dwld,'Force',force,'OUT_file',outfile); 
end

% apply pipeline
[SPdata_o,RT14j_woc,RT14j,RT14h2_bk1_o,RT14h2_bk2_o,BPdata1,BPdata2]...
          = minipipeline_calibration_IR_SP_wCDRSP_yuki(...
            SPdata,bkoption,'DWLD',dwld,'SAVE_MEMORY',save_mem,'APBPRMVL',apbprmvl,...
            'MEAN_DN14',mean_DN14,'SATURATION_RMVL',saturation_rmvl,'MEAN_ROBUST',mean_robust,...
            'BK_SATURATION_RMVL',bk_saturation_rmvl,'BK_BPRMVL',bk_bprmvl,...
            'BK_MEAN_ROBUST',bk_mean_robust,'BK_MEAN_DN14',bk_mean_DN14);
%
%
% then apply 1/(1+MP*SC)
[SPdataMP,SSdataMP,SHdataMP] = selectCDR4MP(SPdataVNIR);
[MP] = calculate_MP(SPdataMP,SSdataMP,SHdataMP);
%SSdata = TRRIFdata.readCDR('SS');
SHdata = SPdata.readCDR('SH');

SHdata.readimg();
SC = SHdata.img(2,:,:);

coef = (1 + MP.* SC);

RT14j_woc_mod = RT14j_woc .* coef;
RT14j_mod = RT14j .* coef;
RT14h2_bk1_o_mod = RT14h2_bk1_o .* coef;
RT14h2_bk2_o_mod = RT14h2_bk2_o .* coef;

SPdata_o.img = RT14j_woc_mod;


end
