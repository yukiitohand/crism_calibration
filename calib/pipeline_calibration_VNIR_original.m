function [RDn,RDn_woc] = pipeline_calibration_VNIR_original(TRRIFdata,EDRdata,DFdata1,DFdata2,BIdata1,BIdata2,bkoption,varargin)
%   Pipeline for the calibration of the CRISM images
%  Input Parameters
%   TRRIFdata, EDRdata: TRRIFdata is just for the information of CDR and
%                       other information
%   DFdata1,DFdata2, prior and post DF measurements
%   BKdata1,BKdata2, prior and post Background measurements
%   bkoption: option for the use of background image. {1,2}
%             1: linear background estimation using prior and post dark
%                measurements.
%             2: flat background estimation using only prior dark
%                measurements.
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
save_mem = false;
apbprmvl = 'HighOrd';
saturation_rmvl = 0;
bkgd_robust = false;
SPdata_o = [];
bk_meanDN14 = false;
bk_saturation_rmvl = 1;
if (rem(length(varargin),2)==1)
    error('Optional parameters should always go by pairs');
else
    for i=1:2:(length(varargin)-1)
        switch upper(varargin{i})
            case 'SPDATA_O'
                SPdata_o = varargin{i+1};
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
            case {'BKGD_ROBUST','BK_MEAN_ROBUST'}
                bkgd_robust = varargin{i+1};
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
frame_rate = TRRIFdata.lbl.MRO_FRAME_RATE{1};
rate_id = get_frame_rate_id(frame_rate);
binx = TRRIFdata.lbl.PIXEL_AVERAGING_WIDTH;

DN = EDRdata.readimg();
if saturation_rmvl
    flg_dsat = (DN==4095); % added by Yuki Feb.18, 2019 
end

% it was expected that saturated pixels are taken as nan after the
% subtraction of background, but it seems that sometimes, it doesn't occur 
% maybe because of some are not saturated.
rownum_table = EDRdata.read_ROWNUM_TABLE();
%-------------------------------------------------------------------------%
% first step (DN12 --> DN14)
PPdata = TRRIFdata.readCDR('PP');
[ DN14 ] = DN12toDN14_VNIR( DN,PPdata,rownum_table );

if save_mem
    clear DN;
end

%-------------------------------------------------------------------------%
% process darks
% DFdata1.readimg();
% DFdata2.readimg();
% [ DN14_df1 ] = DN12toDN14_VNIR( DFdata1.img,PPdata,rownum_table );
% [ DN14_df2 ] = DN12toDN14_VNIR( DFdata2.img,PPdata,rownum_table );
% switch bkgd_robust
%     case 0
%         DN14a_df1 = nanmean(DN14_df1(:,:,:),1);
%         DN14a_df2 = nanmean(DN14_df2(:,:,:),1);
%     case 1
%         DN14a_df1 = robust_v2('mean',DN14_df1,1,'NOutliers',4);
%         DN14a_df2 = robust_v2('mean',DN14_df2,1,'NOutliers',4);
%     otherwise
%         error('Undefined bkgd_robust=%d',bkgd_robust);
% end
% 
% BIdata1_o = BIdata1;
% BIdata2_o = BIdata2;
% BIdata1_o.img = DN14a_df1;
% BIdata2_o.img = DN14a_df2;


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
DBdata = TRRIFdata.readCDR('DB');
EBdata = TRRIFdata.readCDR('EB');
HDdata = TRRIFdata.readCDR('HD');
HKdata = TRRIFdata.readCDR('HK');
TRRIFdata.readHKT(); hkt = TRRIFdata.hkt;
%hkt = EDRdata.readHKT();
%hkt = correctHKTwithHD(hkt,HDdata);
%hkt = correctHKTwithHK(hkt,HKdata);
% using Temperature recorded in the label in TRR I/F data 
[ DN14a,BI_m ] = subtract_bias_VNIR_1( DN14,BIdata1,BIdata2,DBdata,EBdata,hkt,rownum_table,TRRIFdata.lbl );
% [ DN14a_o,BI_m_o ] = subtract_bias_VNIR_1( DN14,BIdata1_o,BIdata2_o,DBdata,EBdata,hkt,rownum_table,TRRIFdata.lbl );
% [ DN14a,BI_m ] = subtract_bias( DN14,BIdata,BSdata,DBdata,EBdata,hkt,rownum_table,'BINX',binx);
% if save_mem
%     clear DN14;
% end

%-------------------------------------------------------------------------%
% the third step (remove detector quadrant electronics ghost)
GHdata = TRRIFdata.readCDR('GH');
[ DN14b,sumGhost ] = remove_quadrantGhost_VNIR( DN14a,GHdata,hkt,'BINX',binx );
if save_mem
    clear DN14a;
end

%-------------------------------------------------------------------------%
% fourth step (dark column subtract)
DMdata = TRRIFdata.readCDR('DM');
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

VLdata = TRRIFdata.readCDR('VL');
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
        [DN14c,mask_saturation] = saturation_removal_VNIR(DN14bb,VLdata,flg_dsat,...
        'binx',binx,'rate_id',rate_id);
    otherwise
        error('Saturation option %d is not defined',saturation_rmvl);
end

%-------------------------------------------------------------------------%
% apply bad a priori pixel interpolation
TRRIFdata.readCDR('BP');
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
            [ DN14d,BP ] = apriori_badpixel_removal( DN14c,TRRIFdata.cdr.BP(1),TRRIFdata.cdr.BP(2),DMdata,'InterpOpt',1 );
        else
            error('Please implement for a priori bad pixel removal for BP more than 2');
        end
    case 'NONE'
        DN14d = DN14c;
end

%-------------------------------------------------------------------------%
% fifth step (nonlinearity correction)
LCdata = TRRIFdata.readCDR('LC');
[ DN14g ] = nonlinearity_correction( DN14d,LCdata,hkt,'BINX',binx );
[ DN14g_woc ] = nonlinearity_correction( DN14c,LCdata,hkt,'BINX',binx );
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
if strcmpi(TRRIFdata.prop.obs_class_type,'MSP') && (TRRIFdata.get_sclk_stop < 850327000)
    error('Stray light correction is not implemented for a special case: MSP & sclk<85032700');
else
    [ RT14j,SL ] = straylight_correction( RT14g,DMdata,rownum_table );
end
RT14j_woc = RT14g_woc - SL;
%%
%-------------------------------------------------------------------------%
% shutter mirror nonrepeatability correction
TRRIFdata.readCDR('SP');
for i=1:length(TRRIFdata.cdr.SP)
    spdata = TRRIFdata.cdr.SP(i);
    spdata_prop = getProp_basenameCDR4(spdata.basename);
    switch upper(spdata_prop.sensor_id)
%         case 'L'
%             SPdata = spdata;
        case 'S'
            SPdata = spdata;
        otherwise
            error('sensor_id %s is wrong',sensor_id);
    end
end



% [SR] = calculate_SR(SSdata,SPdata,SHdata,0);

%-------------------------------------------------------------------------%
% calculate spectroradiometric responsitivity
% [SPdataMP,SSdataMP,SHdataMP] = selectCDR4MP(SPdata);
if isempty(SPdata_o), SPdata_o = SPdata; end
SSdata = TRRIFdata.readCDR('SS');
SHdata = TRRIFdata.readCDR('SH');

% propSS8 = SSdata.prop;
% propSS8.version = 8;
% basenameSS8 = get_basenameCDR4_fromProp(propSS8);
% SSdata8 = CRISMdata(basenameSS8,'');
% SSdata8.readimg();

[MP] = calculate_MP(SPdata_o,SSdata,SHdata);
rowNumTableRSPj = SPdata_o.read_ROWNUM_TABLE();
TDdata = TRRIFdata.readCDR('TD');
[RSPj] = calculate_RSP_VNIR(SPdata_o,SSdata,SHdata,TDdata,0,rowNumTableRSPj);


%-------------------------------------------------------------------------%
% apply binning
DMdata = TRRIFdata.readCDR('DM');
[RSPl] = binning_RSP(RSPj,DMdata,rowNumTableRSPj,'BINX',binx);

%-------------------------------------------------------------------------%
% correct to radiance with the binned responsitivity and flat fielding
NUdata = TRRIFdata.readCDR('NU');
% propNU7 = NUdata.prop;
% propNU7.version = 3;
% propNU7.version = 7;
% basenameNU7 = get_basenameCDR4_fromProp(propNU7);
% NUdata7 = CRISMdata(basenameNU7,'');
[RDm,FF] = calculate_RD_VNIR(RT14j,RSPl,SPdata,SSdata,NUdata,rowNumTableRSPj);
[RDm_woc,FF_woc] = calculate_RD_VNIR(RT14j,RSPl,SPdata,SSdata,NUdata,rowNumTableRSPj);

% [RDm,FF] = calculate_RD(RT14j,RSPl,NUdata);
% [RDm_woc,FF] = calculate_RD(RT14j_woc,RSPl,NUdata);

%-------------------------------------------------------------------------%
% additional destriping step is added.
% [RDmm,coeff] = additional_stripe_removal_VNIR(RDm,DMdata,rownum_table);
% [RDmm_woc,coeff] = additional_stripe_removal_VNIR(RDm_woc,DMdata,rownum_table);

%[RDm_bk1_o,FF_bk1] = calculate_RD(RT14h2_bk1_o,RSPl,NUdata);
%[RDm_bk2_o,FF_bk2] = calculate_RD(RT14h2_bk2_o,RSPl,NUdata);
if save_mem
    clear RT14j RT14j_woc;
end

RDn = apply_DM(RDm,DMdata);
RDn_woc = apply_DM(RDm_woc,DMdata);
%RDn_bk1_o = apply_DM(RDm_bk1_o,DMdata);
%RDn_bk2_o = apply_DM(RDm_bk2_o,DMdata);

end