function [RT14g_bkgd,BKdata_o] = minipipeline_calibration_IR_BK_wCDRBK_yuki(CDRBKdata,CDRBIdata,CDRBPdata,varargin)
% [RT14g_bkgd,BKdata_o] = minipipeline_calibration_IR_BK_wCDRBK_yuki(CDRBKdata,CDRBIdata,CDRBPdata,varargin)
%  re-calculate Background CDR using SOURCE_PRODUCTS stored in lbl of CDR
%  BKdata.
%   INPUTS
%    CDRBKdata: CRISMdata object of the CDR BK data
%    BIdata: CRISMdata object of CDR BI data (estimation should be done
%            outside for now.
%    BPdata: CRISMdata object of CDR BP data (estimation should be done
%            outside for now.
%   OUTPUTS
%    RT14g_bkgd: produced bakground image [1,S,B] (S: samples,B: bands) 
%    BKdata_o: BKdata that stores processed image at img. Band inverse is
%               not performed.
%   OPTIONAL PARAMETERS
%   'DN4095_RMVL': binary, whether or not to perform replacement of saturated
%                  pixels or not.
%                  (default) true
%   'MEAN_ROBUST': integer {0,1}, mode for how mean operation is performed.
%        0: DN14e_df = nanmean(DN14d_df(:,:,:),1);
%        1: DN14e_df = robust_v2('mean',DN14d_df,1,'NOutliers',2);
%      (default) 1
%   'BPRMVL'     : binary, whether or not to perform bad pixel removal. 
%                  (default) false     
%   'MEAN_DN14'  : binary,when mean operation is performed
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
dn4095_rmvl = true;
mean_robust = 1;
bprmvl = false;
mean_DN14 = true;
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
            case 'MEAN_DN14'
                mean_DN14 = varargin{i+1};
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

if isempty(CDRBKdata.basenamesCDR)
    CDRBKdata.load_basenamesCDR('Download',dwld,'Force',force,'OUT_file',outfile);
end
if isempty(CDRBKdata.basenames_SOURCE_OBS)
    CDRBKdata.load_basenames_SOURCE_OBS('Download',dwld,'Force',force,'OUT_file',outfile); 
end

%-------------------------------------------------------------------------%
% Read associated EDR DF image and read all the relevant CDRs
DFdata = CDRBKdata.read_SOURCE_OBS('DF');
DFdata.download(dwld);
PPdata = CDRBKdata.readCDR('PP');
BSdata = CDRBKdata.readCDR('BS'); DBdata = CDRBKdata.readCDR('DB');
EBdata = CDRBKdata.readCDR('EB'); HDdata = CDRBKdata.readCDR('HD');
HKdata = CDRBKdata.readCDR('HK');
GHdata = CDRBKdata.readCDR('GH');
LCdata = CDRBKdata.readCDR('LC');
DMdata = CDRBKdata.readCDR('DM');

%-------------------------------------------------------------------------%
[RT14g_bkgd,BKdata_o] = minipipeline_calibration_IR_BK_yuki(...
    DFdata,PPdata,BSdata,DBdata,EBdata,HDdata,HKdata,CDRBIdata,DMdata,CDRBPdata,...
    GHdata,LCdata,'DN4095_RMVL',dn4095_rmvl,'BPRMVL',bprmvl,...
    'MEAN_ROBUST',mean_robust,'MEAN_DN14',mean_DN14);


end