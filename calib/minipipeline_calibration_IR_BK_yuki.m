function [RT14g_bkgd,BKdata_o,RT14g_df_all] = minipipeline_calibration_IR_BK_yuki(...
    DFdata,PPdata,BSdata,DBdata,EBdata,HDdata,HKdata,BIdata,DMdata,BPdata,GHdata,VLdata,LCdata,varargin)
% [RT14g_bkgd,BKdata_o,RT14g_df_all] = minipipeline_calibration_IR_BK_yuki(...
%    DFdata,PPdata,BSdata,DBdata,EBdata,HDdata,HKdata,BIdata,DMdata,BPdata,GHdata,VLdata,LCdata,varargin)
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
%   'DN4095_RMVL': binary, whether or not to perform replacement of saturated
%                  pixels or not.
%                  (default) false
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

dn4095_rmvl = false;
mean_robust = 1;
bprmvl = false;
mean_DN14 = true;
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
            otherwise
                % Hmmm, something wrong with the parameter string
                error(['Unrecognized option: ''' varargin{i} '''']);
        end
    end
end

%-------------------------------------------------------------------------%
% Read associated EDR DF image
DN12_df = DFdata.readimg();

frame_rate = DFdata.lbl.MRO_FRAME_RATE{1};
rate_id = get_frame_rate_id(frame_rate);
binx = DFdata.lbl.PIXEL_AVERAGING_WIDTH;

%-------------------------------------------------------------------------%
% flag saturated pixels
if dn4095_rmvl
    DFmask4095 = (DN12_df==4095);
end

%-------------------------------------------------------------------------%
% first step (DN12 --> DN14)
rownum_table_df = DFdata.read_ROWNUM_TABLE();
[ DN14_df ] = DN12toDN14( DN12_df,PPdata,rownum_table_df );

%-------------------------------------------------------------------------%
% second step (subtract bias)
hkt_df = DFdata.readHKT();
hkt_dfc = correctHKTwithHD(hkt_df,HDdata);
hkt_dfcc = correctHKTwithHK(hkt_dfc,HKdata);
[ DN14a_df,BI_m_df ] = subtract_bias( DN14_df,BIdata,BSdata,DBdata,EBdata,...
    hkt_dfcc,rownum_table_df,'BINX',binx );

%-------------------------------------------------------------------------%
% the third step (remove detector quadrant electronics ghost)
[ DN14b_df ] = remove_quadrantGhost( DN14a_df,GHdata,hkt_df,'BINX',binx );

%-------------------------------------------------------------------------%
% replace saturated pixel
% now uses VLdata
if dn4095_rmvl
    [DN14c_df,mask_saturation,mask_dead] = saturation_removal(DN14b_df,VLdata,DFmask4095,...
    'binx',binx,'rate_id',rate_id,'is_sphere',false);
    % DN14b_df(DFmask4095) = nan;
end

%-------------------------------------------------------------------------%
% bad pixel removal
if bprmvl
    [ DN14d_df,BP ] = apriori_badpixel_removal( DN14c_df,BPdata,BPdata,DMdata,'InterpOpt',1 );
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
            DN14e_df = robust_v2('mean',DN14d_df,1,'NOutliers',2);
        otherwise
            error('Undefined bkgd_robust=%d',bkgd_robust);
    end
else
    DN14e_df = DN14d_df;
end

%-------------------------------------------------------------------------%
% fourth step (nonlinearity correction)
[ DN14g_df ] = nonlinearity_correction( DN14e_df,LCdata,hkt_dfcc,'BINX',binx );
[ DN14g_df_all ] = nonlinearity_correction( DN14d_df,LCdata,hkt_dfcc,'BINX',binx );

%-------------------------------------------------------------------------%
% fifth step (division by exposure time)
[ RT14g_df ] = divide_by_integrationTime( DN14g_df,hkt_dfcc );
[ RT14g_df_all ] = divide_by_integrationTime( DN14g_df_all,hkt_dfcc );

%-------------------------------------------------------------------------%
% taking mean last
if ~mean_DN14
    switch mean_robust
        case 0
            RT14g_bkgd = nanmean(RT14g_df(:,:,:),1);
        case 1
            RT14g_bkgd = robust_v2('mean',RT14g_df,1,'NOutliers',4);
        otherwise
            error('Undefined bkgd_robust=%d',bkgd_robust);
    end
else
    RT14g_bkgd = RT14g_df;
end
BKdata_o = CRISMdata(DFdata.basename,DFdata.dirpath);
BKdata_o.img = RT14g_bkgd;

end