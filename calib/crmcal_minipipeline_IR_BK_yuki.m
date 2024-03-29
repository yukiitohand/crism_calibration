function [RT14g_bkgd,BKdata_o,RT14g_df_all] = crmcal_minipipeline_IR_BK_yuki(...
    DFdata,PPdata,BSdata,DBdata,EBdata,HDdata,HKdata,BIdata,DMdata,BPdata,GHdata,VLdata,LCdata,varargin)
% [RT14g_bkgd,BKdata_o,RT14g_df_all] = crmcal_minipipeline_IR_BK_yuki(...
%     DFdata,PPdata,BSdata,DBdata,EBdata,HDdata,HKdata,BIdata,DMdata,BPdata,GHdata,VLdata,LCdata,varargin)
%  re-calculate Background from DF image.
%   INPUTS
%    DFdata: CRISMdata object of the DF data
%    PPdara:
%    BSdata:
%    DBdata:
%    EBdata:
%    HDdata:
%    HKdata:
%    BIdata:
%    DMdata:
%    BPdata:
%    GHdata:
%    VLdata:
%    LCdata:
%   OUTPUTS
%    RT14g_bkgd: produced background image [1,S,B] (S: samples,B: bands) 
%    BKdata_o: BKdata that stores processed image at img. Band inverse is
%               not performed.
%    RT14g_df_all: produced background image [L,S,B] non averaged
%   OPTIONAL PARAMETERS
%   'SATURATiON_RMVL': integer, how to perform replacement of saturated
%           pixles {0,1,2}
%           0: no removal
%           1: digital saturation is removed
%           2: analogue saturation is also removed
%           (default) 1
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

saturation_rmvl = 2;
mean_robust = 1;
bprmvl = false;
mean_DN14 = true;

frame_rate = DFdata.lbl.MRO_FRAME_RATE.value;
rate_id = crism_get_frame_rate_id(frame_rate);
binx = DFdata.lbl.PIXEL_AVERAGING_WIDTH;
binx_bk = binx;
binning_bk = crism_get_binning_id(binx_bk);

if (rem(length(varargin),2)==1)
    error('Optional parameters should always go by pairs');
else
    for i=1:2:(length(varargin)-1)
        switch upper(varargin{i})
            case 'BINNING_BK'
                binning_bk = varargin{i+1};
                binx_bk = crism_get_binx(binning_bk);
            case 'BINX_BK'
                binx_bk = varargin{i+1};
                binning_bk = crism_get_binning_id(binx_bk);
            case 'SATURATION_RMVL'
                saturation_rmvl = varargin{i+1};
            case 'BPRMVL'
                bprmvl = varargin{i+1};
            case 'MEAN_ROBUST'
                mean_robust = varargin{i+1};
            case 'MEAN_DN14'
                mean_DN14 = varargin{i+1};
            otherwise
                error('Unrecognized option: %s', varargin{i});
        end
    end
end

%-------------------------------------------------------------------------%
% Read associated EDR DF image
DN12_df = DFdata.readimg();

% bin the image
if binx>1 && binx_bk>binx
    error('Do you want to bin an already binned image?');
end
DN12_df = crism_bin_image_frames(DN12_df,'binning',binning_bk);

%-------------------------------------------------------------------------%
% flag saturated pixels
if saturation_rmvl
    DFmask4095 = (DN12_df==4095);
end

%-------------------------------------------------------------------------%
% first step (DN12 --> DN14)
rownum_table_df = DFdata.read_ROWNUM_TABLE();
[ DN14_df ] = crmcal_DN12toDN14( DN12_df,PPdata,rownum_table_df );

%-------------------------------------------------------------------------%
% second step (subtract bias)
hkt_df = DFdata.readHKT();
hkt_dfc  = crism_correctHKTwithHD(hkt_df,HDdata);
hkt_dfcc = crism_correctHKTwithHK(hkt_dfc,HKdata);
[ DN14a_df,BI_m_df ] = crmcal_subtract_bias( DN14_df,BIdata,BSdata,DBdata,EBdata,...
    hkt_dfcc,rownum_table_df,'BINX',binx_bk );

%-------------------------------------------------------------------------%
% the third step (remove detector quadrant electronics ghost)
[ DN14b_df ] = crmcal_remove_quadrantGhost( DN14a_df,GHdata,hkt_df,'BINX',binx_bk );

%-------------------------------------------------------------------------%
% replace saturated pixel
% now uses VLdata
switch saturation_rmvl
    case 0
        % no saturation removal
        DN14c_df = DN14b_df;
    case 1
        % only digital saturation is dealt with
        DN14c_df = DN14b_df;
        DN14c_df(DFmask4095) = nan;
    case 2
        % analogue saturation is also dealt with
        [DN14c_df,mask_saturation] = crmcal_saturation_removal(DN14b_df,VLdata,DFmask4095,...
        'binx',binx_bk,'rate_id',rate_id);
    otherwise
        error('Saturation option %d is not defined',saturation_rmvl);
end

%-------------------------------------------------------------------------%
% bad pixel removal
if bprmvl
    [ DN14d_df,BP ] = crmcal_apriori_badpixel_removal( DN14c_df,BPdata,BPdata,DMdata,'InterpOpt',1 );
else
    DN14d_df = DN14c_df;
end

%-------------------------------------------------------------------------%
% taking a mean before nonlinearity correction
if mean_DN14
    switch mean_robust
        case 0
            DN14e_df = nanmean(DN14d_df(:,:,:),1);
        case 1
            DN14e_df = robust_v2('mean',DN14d_df,1,'NOutliers',4);
        otherwise
            error('Undefined bkgd_robust=%d',bkgd_robust);
    end
else
    DN14e_df = DN14d_df;
end

%-------------------------------------------------------------------------%
% fourth step (nonlinearity correction)
[ DN14g_df ]     = crmcal_nonlinearity_correction( DN14e_df,LCdata,hkt_dfcc,'BINX',binx_bk );
[ DN14g_df_all ] = crmcal_nonlinearity_correction( DN14d_df,LCdata,hkt_dfcc,'BINX',binx_bk );

%-------------------------------------------------------------------------%
% fifth step (division by exposure time)
[ RT14g_df ]     = crmcal_divide_by_integrationTime( DN14g_df,hkt_dfcc );
[ RT14g_df_all ] = crmcal_divide_by_integrationTime( DN14g_df_all,hkt_dfcc );

%-------------------------------------------------------------------------%
% taking mean last
if ~mean_DN14
    switch mean_robust
        case 0
            RT14g_bkgd = nanmean(RT14g_df(:,:,:),1);
        case 1
            RT14g_bkgd = robust_v2('mean',RT14g_df,1,'NOutliers',4);
        otherwise
            error('Undefined bkgd_robust=%d',mean_robust);
    end
else
    RT14g_bkgd = RT14g_df;
end
BKdata_o = CRISMdata(DFdata.basename,DFdata.dirpath);
BKdata_o.img = RT14g_bkgd;

end