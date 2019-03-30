function [SPdata_o,RT14j_woc,RT14j,RT14h2_bk1_o,RT14h2_bk2_o] = minipipeline_calibration_IR_SP_yuki(...
    EDRSPdata,DFdata1,DFdata2,BKdata1,BKdata2,BPdata1,BPdata2,BIdata,...
    PPdata,BSdata,DBdata,EBdata,HDdata,HKdata,GHdata,VLdata,DMdata,LCdata,LLdata,...
    bkoption,varargin)
% [RT14j_woc,RT14j,RT14h2_bk1_o,RT14h2_bk2_o] = minipipeline_calibration_IR_SP_yuki(...
%     EDRSPdata,DFdata1,DFdata2,BKdata1,BKdata2,BPdata1,BPdata2,BIdata,...
%     PPdata,BSdata,DBdata,EBdata,HDdata,HKdata,GHdata,VLdata,DMdata,LCdata,LLdata,...
%     bkoption,varargin)
%   Mini pipeline for the calibration of CDR SPdata
%  Input Parameters
%   EDRSPdata: EDR SPdata
%   DFdata1,DFdata2, prior and post DF measurements
%   BKdata1,BKdata2, prior and post Background measurements
%   PPdata,BSdata,DBdata,EBdata,HDdata,HKdata,GHdata,DMdata,LCdata,LLdata
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
%       'SATURATiON_RMVL': integer, how to perform replacement of saturated
%           pixles {0,1,2}
%           0: no removal
%           1: digital saturation is removed
%           2: analogue saturation is also removed
%           (default) 2
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
SPdata_ref = [];
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
            case 'MEAN_ROBUST'
                mean_robust = varargin{i+1};
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
            case 'SPDATA_REF'
                SPdata_ref = varargin{i+1};
            otherwise
                % Hmmm, something wrong with the parameter string
                error(['Unrecognized option: ''' varargin{i} '''']);
        end
    end
end

%-------------------------------------------------------------------------%
% actual processing
%-------------------------------------------------------------------------%
frame_rate = EDRSPdata.lbl.MRO_FRAME_RATE{1};
rate_id = get_frame_rate_id(frame_rate);
binx = EDRSPdata.lbl.PIXEL_AVERAGING_WIDTH;

DN = EDRSPdata.readimg();
if saturation_rmvl
    flg_dsat = (DN==4095); % added by Yuki Feb.18, 2019 
end

% it was expected that saturated pixels are taken as nan after the
% subtraction of background, but it seems that sometimes, it doesn't occur 
% maybe because of some are not saturated.
rownum_table = EDRSPdata.read_ROWNUM_TABLE();
%-------------------------------------------------------------------------%
% first step (DN12 --> DN14)
[ DN14 ] = DN12toDN14( DN,PPdata,rownum_table );

if save_mem
    clear DN;
end
%-------------------------------------------------------------------------%
% second step (subtract bias)
hkt = EDRSPdata.readHKT();
hkt = correctHKTwithHD(hkt,HDdata);
hkt = correctHKTwithHK(hkt,HKdata);
% using Temperature recorded in the label in TRR I/F data 
%[ DN14a,BI_m ] = subtract_bias_1( DN14,BIdata,BSdata,DBdata,EBdata,hkt,rownum_table,[],'BINX',binx);
[ DN14a,BI_m ] = subtract_bias( DN14,BIdata,BSdata,DBdata,EBdata,hkt,rownum_table,'BINX',binx);
if save_mem
    clear DN14;
end
%-------------------------------------------------------------------------%
% the third step (remove detector quadrant electronics ghost)
[ DN14b,sumGhost ] = remove_quadrantGhost( DN14a,GHdata,hkt,'BINX',binx );
if save_mem
    clear DN14a;
end

%-------------------------------------------------------------------------%
% flag saturated pixels
%VLdata = TRRIFdata.readCDR('VL');
% added by Yuki Itoh.
% saturation removal is performed after detector quadrant ghost removal.
% this way is manually defined by Yuki. Doesn't follow the direction in the
% crism_dpsis.pdf
% Updated Mar.05 2019 by using VLdata
switch saturation_rmvl
    case 0
        % no saturation removal
        DN14c = DN14b;
    case 1
        % only digital saturation is dealt with
        DN14c = DN14b;
        DN14c(flg_dsat) = nan;
    case 2
        % analogue saturation is also dealt with
        [DN14c,mask_saturation] = saturation_removal(DN14b,VLdata,flg_dsat,...
        'binx',binx,'rate_id',rate_id);
    otherwise
        error('Saturation option %d is not defined',saturation_rmvl);
end
% if bk_saturation_rmvl
% %     DN14d = DN14c;
% %     DN14d(flg_saturation) = nan;
%     [DN14c,mask_saturation] = saturation_removal(DN14b,VLdata,flg_dsat,...
%         'binx',binx,'rate_id',rate_id);
%     %DN14c_woc = DN14b;
%     %DN14d_woc(flg_saturation) = nan;
% end

%-------------------------------------------------------------------------%
% apply bad a priori pixel interpolation
switch upper(apbprmvl)
    case 'HIGHORD'
        [ DN14d,BP ] = apriori_badpixel_removal( DN14c,BPdata1,BPdata2,DMdata,'InterpOpt',1 );
    case 'NONE'
        DN14d = DN14c;
end

%-------------------------------------------------------------------------%
% step 3.5 (take a mean over the obtained DN14 data)
% This step is specific for SP processing pipeline
if mean_DN14
    switch mean_robust
        case 0
            DN14e_woc = nanmean(DN14c(:,:,:),1);
            DN14e = nanmean(DN14d(:,:,:),1);
        case 1
            DN14e_woc = robust_v2('mean',DN14c,1,'NOutliers',10);
            DN14e = robust_v2('mean',DN14d,1,'NOutliers',10);
        otherwise
            error('Undefined mean_robust=%d',mean_robust);
    end
else
    DN14e_woc = DN14c;
    DN14e = DN14d;
end

%-------------------------------------------------------------------------%
% fourth step (nonlinearity correction)
[ DN14g ] = nonlinearity_correction( DN14e,LCdata,hkt,'BINX',binx );
[ DN14g_woc ] = nonlinearity_correction( DN14e_woc,LCdata,hkt,'BINX',binx );
if save_mem
    clear DN14c DN14b;
end

%-------------------------------------------------------------------------%
% fifth step (division by exposure time)
[ RT14g ] = divide_by_integrationTime( DN14g,hkt );
[ RT14g_woc ] = divide_by_integrationTime( DN14g_woc,hkt );
% [ RT14g ] = divide_by_integrationTime( DN14c,hkt );
% [ RT14g_woc ] = divide_by_integrationTime( DN14b,hkt );
if save_mem
    clear DN14g DN14g_woc;
end


%%
%-------------------------------------------------------------------------%
% process darks from scratch
[RT14g_df1_1,BKdata1_o,RT14g_df1] = minipipeline_calibration_IR_BK_yuki(...
    DFdata1,PPdata,BSdata,DBdata,EBdata,HDdata,HKdata,BIdata,DMdata,...
    BPdata1,GHdata,VLdata,LCdata,'SATURATION_RMVL',bk_saturation_rmvl,'BPRMVL',bk_bprmvl,...
    'MEAN_ROBUST',bk_mean_robust,'MEAN_DN14',bk_mean_DN14);
[RT14g_df1_2,BKdata2_o,RT14g_df2] = minipipeline_calibration_IR_BK_yuki(...
    DFdata2,PPdata,BSdata,DBdata,EBdata,HDdata,HKdata,BIdata,DMdata,...
    BPdata2,GHdata,VLdata,LCdata,'SATURATION_RMVL',bk_saturation_rmvl,'BPRMVL',bk_bprmvl,...
    'MEAN_ROBUST',bk_mean_robust,'MEAN_DN14',bk_mean_DN14);
hkt_df1 = DFdata1.readHKT();
hkt_df1c = correctHKTwithHD(hkt_df1,HDdata);
hkt_df1cc = correctHKTwithHK(hkt_df1c,HKdata);
hkt_df2 = DFdata2.readHKT();
hkt_df2c = correctHKTwithHD(hkt_df2,HDdata);
hkt_df2cc = correctHKTwithHK(hkt_df2c,HKdata);
%%
%-------------------------------------------------------------------------%
% background subtraction
switch bkoption
    case 1
        [ RT14h,Bkgd ] = background_subtraction( RT14g,BKdata1,BKdata2,hkt );
        [ RT14h_woc ] = background_subtraction( RT14g_woc,BKdata1_o,BKdata2_o,hkt );
        [ RT14h_bk1_o ] = background_subtraction( RT14g_df1,BKdata1_o,BKdata2_o,hkt_df1cc );
        [ RT14h_bk2_o ] = background_subtraction( RT14g_df2,BKdata1_o,BKdata2_o,hkt_df2cc );
    case 2
        [ RT14h,Bkgd ] = background_subtraction_v2( RT14g,BKdata1,BKdata2,hkt );
        [ RT14h_woc ] = background_subtraction_v2( RT14g_woc,BKdata1_o,BKdata2_o,hkt );
        [ RT14h_bk1_o ] = background_subtraction_v2( RT14g_df1,BKdata1_o,BKdata2_o,hkt_df1cc );
        [ RT14h_bk2_o ] = background_subtraction_v2( RT14g_df2,BKdata1_o,BKdata2_o,hkt_df2cc );
end
if save_mem
    clear RT14g RT14g_woc;
end
%-------------------------------------------------------------------------%
% dark column subtract
% I think additional dark column subtract is applied somewhere
[RT14h2,dc] = dark_column_subtract(RT14h,DMdata);
[RT14h2_woc,dc] = dark_column_subtract(RT14h_woc,DMdata);
[RT14h2_bk1_o,dc_bk1] = dark_column_subtract(RT14h_bk1_o,DMdata);
[RT14h2_bk2_o,dc_bk2] = dark_column_subtract(RT14h_bk2_o,DMdata);
% if save_mem
%     clear RT14h RT14h_woc;
% end

%-------------------------------------------------------------------------%
% second order light removal
[RT14j,K] = subtract_highorderlight(RT14h2,LLdata);
RT14j_woc = RT14h2_woc - K;
% [RT14j,K] = subtract_highorderlight(RT14h,LLdata);
% RT14j_woc = RT14h_woc - K;
if save_mem
    clear RT14h2 RT14h2_woc;
end

%-------------------------------------------------------------------------%
% dead pixel removal
[RT14jj,mask_dead] = deadpixel_removal(RT14j,VLdata,DMdata,'binx',binx,'rate_id',rate_id);
[RT14jj_woc,mask_dead] = deadpixel_removal(RT14j_woc,VLdata,DMdata,'binx',binx,'rate_id',rate_id);

%-------------------------------------------------------------------------%
% mean if 
if ~mean_DN14
    switch mean_robust
        case 0
            RT14j_SP_woc = nanmean(RT14jj_woc(:,:,:),1);
            RT14j_SP= nanmean(RT14jj(:,:,:),1);
        case 1
            RT14j_SP_woc = robust_v2('mean',RT14jj_woc,1,'NOutliers',10);
            RT14j_SP = robust_v2('mean',RT14jj,1,'NOutliers',10);
        otherwise
            error('Undefined mean_robust=%d',mean_robust);
    end
else
    RT14j_SP_woc = RT14j_woc;
    RT14j_SP = RT14j;
end


SPdata_o = CRISMdata(EDRSPdata.basename,'');
SPdata_o.img = RT14j_SP_woc;

end