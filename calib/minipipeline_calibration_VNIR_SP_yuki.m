function [SPdata_o,RT14j_woc,RT14j] = minipipeline_calibration_VNIR_SP_yuki( EDRSPdata,DFdata1,DFdata2,...
    PPdata,DBdata,EBdata,HDdata,HKdata,GHdata,VLdata,DMdata,LCdata,varargin)
%   Pipeline for the calibration of the CRISM images for producing SPdata
%  Input Parameters

%   Optional Parameters

save_mem = false;
apbprmvl = 'None';
saturation_rmvl = 0;
mean_robust = 1;
bk_mean_robust = 1;
bk_meanDN14 = false;
bk_saturation_rmvl = 1;
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
            case 'SATURATION_RMVL'
                saturation_rmvl = varargin{i+1};
            case 'BK_SATURATION_RMVL'
                bk_saturation_rmvl = varargin{i+1};
            case {'ROBUST','MEAN_ROBUST'}
                mean_robust = varargin{i+1};
                if bkoption==1
                    error('no effect of "Bkgd_robust" when bkoption=%d',bkoption);
                end
            case 'BK_MEAN_DN14'
                bk_meanDN14 = varargin{i+1};
            otherwise
                % Hmmm, something wrong with the parameter string
                error(['Unrecognized option: ''' varargin{i} '''']);
        end
    end
end
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
% PPdata = TRRIFdata.readCDR('PP');
[ DN14 ] = DN12toDN14_VNIR( DN,PPdata,rownum_table );

if save_mem
    clear DN;
end

%-------------------------------------------------------------------------%
% process darks
DFdata1.readimg();
DFdata2.readimg();
[ DN14_df1 ] = DN12toDN14_VNIR( DFdata1.img,PPdata,rownum_table );
[ DN14_df2 ] = DN12toDN14_VNIR( DFdata2.img,PPdata,rownum_table );
switch bk_mean_robust
    case 0
        DN14a_df1 = nanmean(DN14_df1(:,:,:),1);
        DN14a_df2 = nanmean(DN14_df2(:,:,:),1);
    case 1
        DN14a_df1 = robust_v2('mean',DN14_df1,1,'NOutliers',4);
        DN14a_df2 = robust_v2('mean',DN14_df2,1,'NOutliers',4);
    otherwise
        error('Undefined bkgd_robust=%d',bkgd_robust);
end

BIdata1_o = DFdata1;
BIdata2_o = DFdata2;
BIdata1_o.img = DN14a_df1;
BIdata2_o.img = DN14a_df2;
hkt_df1 = DFdata1.readHKT();
hkt_df1 = correctHKTwithHD(hkt_df1,HDdata);
hkt_df1 = correctHKTwithHK(hkt_df1,HKdata);
BIdata1_o.hkt = hkt_df1;
BIdata1_o.lbl.MRO_DETECTOR_TEMPERATURE = nanmean(cat(1,BIdata1_o.hkt.data.VNIR_DETECTOR_TEMP1));
BIdata1_o.lbl.MRO_FPE_TEMPERATURE = nanmean(cat(1,BIdata1_o.hkt.data.VNIR_FPU_BOARD_TEMP));

hkt_df2 = DFdata2.readHKT();
hkt_df2 = correctHKTwithHD(hkt_df2,HDdata);
hkt_df2 = correctHKTwithHK(hkt_df2,HKdata);
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
hkt = correctHKTwithHD(hkt,HDdata);
hkt = correctHKTwithHK(hkt,HKdata);
% using Temperature recorded in the label in TRR I/F data 
% [ DN14a,BI_m ] = subtract_bias_VNIR_1( DN14,BIdata1,BIdata2,DBdata,EBdata,hkt,rownum_table,TRRIFdata.lbl );
[ DN14a,BI_m ] = subtract_bias_VNIR_1( DN14,BIdata1_o,BIdata2_o,DBdata,EBdata,hkt,rownum_table,EDRSPdata.lbl );
% [ DN14a,BI_m ] = subtract_bias( DN14,BIdata,BSdata,DBdata,EBdata,hkt,rownum_table,'BINX',binx);
% if save_mem
%     clear DN14;
% end

%-------------------------------------------------------------------------%
% the third step (remove detector quadrant electronics ghost)
% GHdata = TRRIFdata.readCDR('GH');
[ DN14b,sumGhost ] = remove_quadrantGhost_VNIR( DN14a,GHdata,hkt,'BINX',binx );
if save_mem
    clear DN14a;
end

%-------------------------------------------------------------------------%
% fourth step (dark column subtract)
% DMdata = TRRIFdata.readCDR('DM');
[DN14bb,dc] = dark_column_subtract(DN14b,DMdata);
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
        [DN14c,mask_saturation] = saturation_removal(DN14bb,VLdata,flg_dsat,...
        'binx',binx,'rate_id',rate_id);
    otherwise
        error('Saturation option %d is not defined',saturation_rmvl);
end

%-------------------------------------------------------------------------%
% apply bad a priori pixel interpolation
% TRRIFdata.readCDR('BP');
% switch EDRdata.lbl.OBSERVATION_TYPE
%     case {'FRT','HRL','HRS','MSP','HSP'}
%         for i=1:length(TRRIFdata.cdr.BP)
%             bpdata = TRRIFdata.cdr.BP(i);
%             if ~any(strcmpi(EDRdata.basename,bpdata.lbl.SOURCE_PRODUCT_ID))
%                 if any(strcmpi(DFdata1.basename,bpdata.lbl.SOURCE_PRODUCT_ID))
%                     BPdata1 = bpdata;
%                 elseif any(strcmpi(DFdata2.basename,bpdata.lbl.SOURCE_PRODUCT_ID))
%                     BPdata2 = bpdata;
%                 end
%             else
%                 BPdata_post = bpdata;
%             end
%         end
%     case {'FRS','ATO'}
%         % in case of FRS, DFdata1 and DFdata2 are same.
%         for i=1:length(TRRIFdata.cdr.BP)
%             bpdata = TRRIFdata.cdr.BP(i);
%             if ~any(strcmpi(EDRdata.basename,bpdata.lbl.SOURCE_PRODUCT_ID))
%                 if any(strcmpi(DFdata1.basename,bpdata.lbl.SOURCE_PRODUCT_ID))
%                     BPdata1 = bpdata; BPdata2 = bpdata;
%                 end
%             else
%                 BPdata_post = bpdata;
%             end
%         end
%     otherwise
%         error('Undefined observation type %s.',EDRdata.lbl.OBSERVATION_TYPE);
% end

% DMdata = TRRIFdata.readCDR('DM');
% [ DN14c,BP ] = apriori_badpixel_removal( DN14b,BPdata1,BPdata2,DMdata,'InterpOpt',1 );

switch upper(apbprmvl)
    case 'HIGHORD'
        if length(TRRIFdata.cdr.BP)==2
            [ DN14d,BP ] = apriori_badpixel_removal( DN14c,BPdata1,BPdata2,DMdata,'InterpOpt',1 );
        else
            error('Please implement for a priori bad pixel removal for BP more than 2');
        end
    case 'NONE'
        DN14d = DN14c;
end

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
[ DN14g ] = nonlinearity_correction( DN14e,LCdata,hkt,'BINX',binx );
[ DN14g_woc ] = nonlinearity_correction( DN14e_woc,LCdata,hkt,'BINX',binx );
if save_mem
    clear DN14c DN14b;
end


%-------------------------------------------------------------------------%
% sixth step (division by exposure time)
[ RT14g ] = divide_by_integrationTime( DN14g,hkt );
[ RT14g_woc ] = divide_by_integrationTime( DN14g_woc,hkt );
if save_mem
    clear DN14g DN14g_woc;
end

%-------------------------------------------------------------------------%
% stray light correction
[ RT14j,SL ] = straylight_correction( RT14g,DMdata,rownum_table );
RT14j_woc = RT14g_woc - SL;

SPdata_o = CRISMdata(EDRSPdata.basename,'');
SPdata_o.img = RT14j_woc;

end