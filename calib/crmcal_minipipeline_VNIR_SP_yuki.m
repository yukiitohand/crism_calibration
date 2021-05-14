function [SPdata_o,RT14j_woc,RT14j] = crmcal_minipipeline_VNIR_SP_yuki(...
    EDRSPdata,DFdata1,DFdata2,...
    PPdata,DBdata,EBdata,HDdata,HKdata,GHdata,VLdata,DMdata,LCdata,varargin)
%   Pipeline for the calibration of the CRISM VNIR images for producing 
%   SPdata
%  INPUTS
%   EDRSPdata: EDR SP data
%   DFdata1:  prior Dark measurement
%   DFdata2:  post Dark measurement
%   PPdata,DBdata,EBdata,HDdata,HKdata,GHdata,VLdata,DMdata,LCdata
%  OUTPUTS
%   SPdata_o
%   RT14j_woc
%   RT14j
%  OPTIONAL PARAMETERS
%   'SAVE_MEMORY'
%      saving memory or not. (default) false
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
%
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
apbprmvl = 'None';
saturation_rmvl = 2;
mean_robust = 1;
bk_mean_robust = 1;
% mean_DN14 = true;

binning_sp = 0;
binx_sp = 1;

if (rem(length(varargin),2)==1)
    error('Optional parameters should always go by pairs');
else
    for i=1:2:(length(varargin)-1)
        switch upper(varargin{i})
            case 'BINNING_SP'
                binning_sp = varargin{i+1};
                binx_sp = crism_get_binning(binning_sp);
            case 'BINX_SP'
                binx_sp = varargin{i+1};
                binning_sp = crism_get_binning_id(binx_sp);
            case 'SAVE_MEMORY'
                save_mem = varargin{i+1};
            %case 'APBPRMVL'
            %    apbprmvl = varargin{i+1};
            %    if ~any(strcmpi(apbprmvl,{'HighOrd','None'}))
            %        error('apbprmvl (%s) should be either {"HighOrd","None"}',apbprmvl);
            %    end
            % case 'MEAN_DN14'
            %    mean_DN14 = varargin{i+1};
            case 'SATURATION_RMVL'
                saturation_rmvl = varargin{i+1};
            case 'MEAN_ROBUST'
                mean_robust = varargin{i+1};
            case 'BK_MEAN_ROBUST'
                mean_robust = varargin{i+1};
            otherwise
                error('Unrecognized keyword: %s',varargin{i});
        end
    end
end
frame_rate = EDRSPdata.lbl.MRO_FRAME_RATE.value;
rate_id = crism_get_frame_rate_id(frame_rate);
% binx_sp = EDRSPdata.lbl.PIXEL_AVERAGING_WIDTH;

DN = EDRSPdata.readimg();

% apply binning this early
DN = crism_bin_image_frames(DN,'binning',binning_sp);

if saturation_rmvl
    flg_dsat = (DN==4095);
end

% it was expected that saturated pixels are taken as nan after the
% subtraction of background, but it seems that sometimes, it doesn't occur 
% maybe because of some are not saturated.
rownum_table = EDRSPdata.read_ROWNUM_TABLE();
%-------------------------------------------------------------------------%
% first step (DN12 --> DN14)
% PPdata = TRRIFdata.readCDR('PP');
[ DN14 ] = crmcal_DN12toDN14_VNIR( DN,PPdata,rownum_table );

if save_mem
    clear DN;
end

%-------------------------------------------------------------------------%
% process darks
DFdata1.readimg();
DN12_df1 = crism_bin_image_frames(DFdata1.img,'binning',binning_sp);
% BI is before the ghost correction, so I think saturation removal is not
% necessary.
% if saturation_rmvl
%     flg_dsat_df1 = (DFdata1.img==4095);
%     flg_dsat_df2 = (DFdata2.img==4095);
% end
DFdata2.readimg();
DN12_df2 = crism_bin_image_frames(DFdata2.img,'binning',binning_sp);
[ DN14_df1 ] = crmcal_DN12toDN14_VNIR( DN12_df1,PPdata,rownum_table );
[ DN14_df2 ] = crmcal_DN12toDN14_VNIR( DN12_df2,PPdata,rownum_table );

switch bk_mean_robust
    case 0
        DN14a_df1 = nanmean(DN14_df1(:,:,:),1);
        DN14a_df2 = nanmean(DN14_df2(:,:,:),1);
    case 1
        DN14a_df1 = robust_v2('mean',DN14_df1,1,'NOutliers',4);
        DN14a_df2 = robust_v2('mean',DN14_df2,1,'NOutliers',4);
    otherwise
        error('Undefined bkgd_robust=%d',bk_mean_robust);
end

BIdata1_o = DFdata1;
BIdata2_o = DFdata2;
BIdata1_o.img = DN14a_df1;
BIdata2_o.img = DN14a_df2;
hkt_df1 = DFdata1.readHKT();
hkt_df1 = crism_correctHKTwithHD(hkt_df1,HDdata);
hkt_df1 = crism_correctHKTwithHK(hkt_df1,HKdata);
BIdata1_o.hkt = hkt_df1;
BIdata1_o.lbl.MRO_DETECTOR_TEMPERATURE = nanmean(cat(1,BIdata1_o.hkt.data.VNIR_DETECTOR_TEMP1));
BIdata1_o.lbl.MRO_FPE_TEMPERATURE = nanmean(cat(1,BIdata1_o.hkt.data.VNIR_FPU_BOARD_TEMP));

hkt_df2 = DFdata2.readHKT();
hkt_df2 = crism_correctHKTwithHD(hkt_df2,HDdata);
hkt_df2 = crism_correctHKTwithHK(hkt_df2,HKdata);
BIdata2_o.hkt = hkt_df2;
BIdata2_o.lbl.MRO_DETECTOR_TEMPERATURE = nanmean(cat(1,BIdata2_o.hkt.data.VNIR_DETECTOR_TEMP1));
BIdata2_o.lbl.MRO_FPE_TEMPERATURE = nanmean(cat(1,BIdata2_o.hkt.data.VNIR_FPU_BOARD_TEMP));

%-------------------------------------------------------------------------%
% second step (subtract bias/dark)
% TRRIFdata.readCDR('BI');
% for i=1:length(TRRIFdata.cdr.BI)
%     bidata = TRRIFdata.cdr.BI(i);
%     if bidata.lbl.MRO_FRAME_RATE{1} == frame_rate
%         BIdata = bidata;
%     end
% end
% BSdata = TRRIFdata.readCDR('BS');
% DBdata = TRRIFdata.readCDR('DB');
% EBdata = TRRIFdata.readCDR('EB');
% HDdata = TRRIFdata.readCDR('HD');
% HKdata = TRRIFdata.readCDR('HK');
% TRRIFdata.readHKT(); hkt = TRRIFdata.hkt;
hkt = EDRSPdata.readHKT();
hkt = crism_correctHKTwithHD(hkt,HDdata);
hkt = crism_correctHKTwithHK(hkt,HKdata);
% using Temperature recorded in the label in TRR I/F data 
% [ DN14a,BI_m ] = crmcal_subtract_bias_VNIR_1( DN14,BIdata1,BIdata2,DBdata,EBdata,hkt,rownum_table,TRRIFdata.lbl );
[ DN14a,BI_m ] = crmcal_subtract_bias_VNIR_1( DN14,BIdata1_o,BIdata2_o,DBdata,EBdata,hkt,rownum_table,EDRSPdata.lbl,'binx',binx_sp);
% [ DN14a,BI_m ] = crmcal_subtract_bias( DN14,BIdata,BSdata,DBdata,EBdata,hkt,rownum_table,'BINX',binx);
% if save_mem
%     clear DN14;
% end

%-------------------------------------------------------------------------%
% the third step (remove detector quadrant electronics ghost)
% GHdata = TRRIFdata.readCDR('GH');
[ DN14b,sumGhost ] = crmcal_remove_quadrantGhost_VNIR( DN14a,GHdata,hkt,'BINX',binx_sp );
if save_mem
    clear DN14a;
end

%-------------------------------------------------------------------------%
% fourth step (dark column subtract)
% DMdata = TRRIFdata.readCDR('DM');
[DN14bb,dc] = crmcal_dark_column_subtract(DN14b,DMdata);
if save_mem
    clear DN14b;
end

%-------------------------------------------------------------------------%
% flag saturated pixels
% added by Yuki Itoh.
% saturation removal is performed after detector quadrant ghost removal.
% this way is manually defined by Yuki. Doesn't follow the direction in the
% crism_dpsis.pdf

% VLdata = TRRIFdata.readCDR('VL');
switch saturation_rmvl
    case 0
        % no saturation removal
        DN14c = DN14bb;
    case 1
        % only digital saturation is dealt with
        DN14c = DN14bb;
        DN14c(flg_dsat) = nan;
    case 2
        % analogue saturation is also dealt with
        [DN14c,mask_saturation] = crmcal_saturation_removal_VNIR(DN14bb,VLdata,flg_dsat,...
        'binx',binx_sp,'rate_id',rate_id);
    otherwise
        error('Saturation option %d is not defined',saturation_rmvl);
end

% dead pixel removal
[DN14cc,mask_dead] = crmcal_deadpixel_removal_VNIR(DN14c,VLdata,DMdata,'binx',binx_sp,'rate_id',rate_id);

%-------------------------------------------------------------------------%
% bad pixel removal
switch upper(apbprmvl)
    case 'HIGHORD'
        if length(TRRIFdata.cdr.BP)==2
            [ DN14d,BP ] = crmcal_apriori_badpixel_removal( DN14cc,BPdata1,BPdata2,DMdata,'InterpOpt',1 );
        else
            error('Please implement for a priori bad pixel removal for BP more than 2');
        end
    case 'NONE'
        DN14d = DN14cc;
end

%-------------------------------------------------------------------------%
% taking mean
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

%-------------------------------------------------------------------------%
% fifth step (nonlinearity correction)
% LCdata = TRRIFdata.readCDR('LC');
[ DN14g ]     = crmcal_nonlinearity_correction( DN14e,LCdata,hkt,'BINX',binx_sp );
[ DN14g_woc ] = crmcal_nonlinearity_correction( DN14e_woc,LCdata,hkt,'BINX',binx_sp );
if save_mem
    clear DN14c DN14b;
end


%-------------------------------------------------------------------------%
% sixth step (division by exposure time)
[ RT14g ]     = crmcal_divide_by_integrationTime( DN14g,hkt );
[ RT14g_woc ] = crmcal_divide_by_integrationTime( DN14g_woc,hkt );
if save_mem
    clear DN14g DN14g_woc;
end

%-------------------------------------------------------------------------%
% stray light correction
[ RT14j,SL ] = crmcal_straylight_correction( RT14g,DMdata,rownum_table,'BINX',binx_sp );
RT14j_woc = RT14g_woc - SL;

SPdata_o = CRISMdata(EDRSPdata.basename,'');
SPdata_o.img = RT14j_woc;

end