function [RT14g_bkgd,BKdata_o] = calculate_IR_BK_wBK(BKdata,bkgd_robust)
% [RT14g_bkgd,BKdata_o] = calculate_IR_BK_wBK(BKdata,bkgd_robust)
%  re-calculate Background CDR using SOURCE_PRODUCTS stored in lbl of
%  BKdata.
%   INPUTS
%    BKdata: CRISMdata object of the CDR BK data
%    bkgd_robust: integer, currently {0,1} are supported.
%                 0: namean is used
%                 1: robust_v2('mean',RT14g_df1,1,'NOutliers',2)
%   OUTPUTS
%    RT14g_bkgd: produced bakground image [1,S,B] (S: samples,B: bands) 
%    BKdata_o: BKdata that stores processed image at img. Band inverse is
%               not performed.
%               

if isempty(BKdata.basenamesCDR), BKdata.load_basenamesCDR(); end
if isempty(BKdata.basenames_SOURCE_OBS)
    BKdata.load_basenames_SOURCE_OBS(); 
end

%-------------------------------------------------------------------------%
% Read associated EDR DF image and read all the relevant CDRs
basenameDF = BKdata.basenames_SOURCE_OBS.DF;
DFdata = CRISMdata(basenameDF,'');
PPdata = BKdata.readCDR('PP');
BSdata = BKdata.readCDR('BS'); DBdata = BKdata.readCDR('DB');
EBdata = BKdata.readCDR('EB'); HDdata = BKdata.readCDR('HD');
HKdata = BKdata.readCDR('HK');

BKdata.readCDR('BI');
for i=1:length(BKdata.cdr.BI)
    bidata = BKdata.cdr.BI{i};
    if bidata.lbl.MRO_FRAME_RATE{1} == frame_rate
        BIdata = bidata;
    end
end

GHdata = BKdata.readCDR('GH');
LCdata = TRRIFdata.readCDR('LC');

[RT14g_bkgd,BKdata_o] = calculate_Bkgd_wDF(DFdata,bkgd_robust,...
        PPdata,BSdata,DBdata,EBdata,HDdata,HKdata,BIdata,GHdata,LCdata);

end