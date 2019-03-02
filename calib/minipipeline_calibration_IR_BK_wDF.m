function [RT14g_bkgd,BKdata_o] = minipipeline_calibration_IR_BK_wDF(...
    DFdata,bkgd_robust,PPdata,BSdata,DBdata,EBdata,HDdata,HKdata,BIdata,...
    GHdata,LCdata)
% [RT14g_bkgd,BKdata_o] = minipipeline_calibration_IR_BK_wDF(...
%     DFdata,bkgd_robust,PPdata,BSdata,DBdata,EBdata,HDdata,HKdata,BIdata,...
%     GHdata,LCdata)
%  re-calculate Background from DF image.
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
% Read associated EDR DF image
%basenameDF = BKdata.basenames_SOURCE_OBS.DF;
%DFdata = CRISMdata(basenameDF,'');
DN12_df = DFdata.readimg();

frame_rate = DFdata.lbl.MRO_FRAME_RATE{1};
binx = DFdata.lbl.PIXEL_AVERAGING_WIDTH;

%-------------------------------------------------------------------------%
% flag saturated pixels
DFmask = (DN12_df==4095);

%-------------------------------------------------------------------------%
% first step (DN12 --> DN14)
% PPdata = BKdata.readCDR('PP');
rownum_table_df = DFdata.read_ROWNUM_TABLE();
[ DN14_df ] = DN12toDN14( DN12_df,PPdata,rownum_table_df );


%-------------------------------------------------------------------------%
% second step (subtract bias)
% BSdata = BKdata.readCDR('BS'); DBdata = BKdata.readCDR('DB');
% EBdata = BKdata.readCDR('EB'); HDdata = BKdata.readCDR('HD');
% HKdata = BKdata.readCDR('HK');
hkt_df = DFdata.readHKT();
hkt_dfc = correctHKTwithHD(hkt_df,HDdata);
hkt_dfcc = correctHKTwithHK(hkt_dfc,HKdata);
% select bias
% BKdata.readCDR('BI');
% for i=1:length(BKdata.cdr.BI)
%     bidata = BKdata.cdr.BI{i};
%     if bidata.lbl.MRO_FRAME_RATE{1} == frame_rate
%         BIdata = bidata;
%     end
% end
[ DN14a_df,BI_m_df ] = subtract_bias( DN14_df,BIdata,BSdata,DBdata,EBdata,hkt_dfcc,rownum_table_df,'BINX',binx );

%-------------------------------------------------------------------------%
% the third step (remove detector quadrant electronics ghost)
% GHdata = BKdata.readCDR('GH');
[ DN14b_df ] = remove_quadrantGhost( DN14a_df,GHdata,hkt_df,'BINX',binx );

%-------------------------------------------------------------------------%
% replace saturated pixel
DN14b_df(DFmask) = nan;

%-------------------------------------------------------------------------%
% fourth step (nonlinearity correction)
% LCdata = TRRIFdata.readCDR('LC');
[ DN14g_df ] = nonlinearity_correction( DN14b_df,LCdata,hkt_dfcc,'BINX',binx );

%-------------------------------------------------------------------------%
% fifth step (division by exposure time)
[ RT14g_df ] = divide_by_integrationTime( DN14g_df,hkt_dfcc );

% RT14g_df1(DF1mask) = nan;

%-------------------------------------------------------------------------%
% finally taking a mean.
BKdata_o = CRISMdata(DFdata.basename,DFdata.dirpath);
switch bkgd_robust
    case 0
        RT14g_bkgd = nanmean(RT14g_df(:,:,:),1);
    case 1
        RT14g_bkgd = robust_v2('mean',RT14g_df,1,'NOutliers',2);
    otherwise
        error('Undefined bkgd_robust=%d',bkgd_robust);
end
BKdata_o.img = RT14g_bkgd;

end