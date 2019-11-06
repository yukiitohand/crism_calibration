function [RDn,RDn_woc,RDn_bk1_o,RDn_bk2_o] = pipeline_calibration_IR_IF_yuki(...
    EDRdata,DFdata1,DFdata2,BKdata1,BKdata2,BPdata1,BPdata2,BIdata,...
    PPdata,BSdata,DBdata,EBdata,HDdata,HKdata,GHdata,VLdata,DMdata,LCdata,...
    LLdata,SPdata,SSdata,SHdata,NUdata,...
    SPdata_o,BKdata1_o,BKdata2_o,RT14g_df1,RT14g_df2,bkoption,MP,varargin)
%   Main Pipeline for the calibration of the CRISM IR IF data
% [RDn,RDn_woc,RDn_bk1_o,RDn_bk2_o] = pipeline_calibration_IR_IF_yuki(...
%     EDRdata,DFdata1,DFdata2,BKdata1,BKdata2,BPdata1,BPdata2,BIdata,...
%     PPdata,BSdata,DBdata,EBdata,HDdata,HKdata,GHdata,VLdata,DMdata,LCdata,...
%     LLdata,SPdata,SSdata,SHdata,NUdata,...
%     SPdata_o,BKdata1_o,BKdata2_o,RT14g_df1,RT14g_df2,bkoption,MP,varargin)
%  Input Parameters
%   EDRdata: TRRIFdata is just for the information of CDR and
%                       other information
%   DFdata1,DFdata2, prior and post DF measurements
%   BKdata1,BKdata2, prior and post CDR BK data
%   BPdata1,BPdata2 : CDR BP data
%   BIdata: CDR BI data
%   PPdata,BSdata,DBdata,EBdata,EBdata,HDdata,HKdata,GHdata,VLdata,DMdata,LCdata,...
%   LLdata,SPdata,SSdata,SHdata,NUdata : CDR data
%   
%   SPdata_o : manually processed SP data
%   BKdata1_o: manually processed BKdata1
%   BKdata2_o: manually processed BKdata2
%   RT14g_df1: non averaged dark (processed up to CDR BK stage)
%   RT14g_df2: non averaged dark (processed up to CDR BK stage)
%   
%   MP: scalar, shutter mirror correction parameter
%
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
%     'FLAT_FIELD' : 
%      boolean, whether or not to perform flat field correction using 
%      NUdata or not.
%      (default) true
save_mem = false;
apbprmvl = 'HighOrd';
saturation_rmvl = 2;
flat_field = true;
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
            case 'FLAT_FIELD'
                flat_field = varargin{i+1};
            otherwise
                % Hmmm, something wrong with the parameter string
                error(['Unrecognized option: ''' varargin{i} '''']);
        end
    end
end
frame_rate = EDRdata.lbl.MRO_FRAME_RATE{1};
rate_id = get_frame_rate_id(frame_rate);
binx = EDRdata.lbl.PIXEL_AVERAGING_WIDTH;

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
% PPdata = TRRIFdata.readCDR('PP');
[ DN14 ] = DN12toDN14( DN,PPdata,rownum_table );

if save_mem
    clear DN;
end
%-------------------------------------------------------------------------%
% second step (subtract bias)
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
hkt = EDRdata.readHKT();
hkt = correctHKTwithHD(hkt,HDdata);
hkt = correctHKTwithHK(hkt,HKdata);
% using Temperature recorded in the label in TRR I/F data 
% [ DN14a,BI_m ] = subtract_bias_1( DN14,BIdata,BSdata,DBdata,EBdata,hkt,rownum_table,TRRIFdata.lbl );
[ DN14a,BI_m ] = subtract_bias( DN14,BIdata,BSdata,DBdata,EBdata,hkt,rownum_table,'BINX',binx);
if save_mem
    clear DN14;
end
%-------------------------------------------------------------------------%
% the third step (remove detector quadrant electronics ghost)
% GHdata = TRRIFdata.readCDR('GH');
[ DN14b,sumGhost ] = remove_quadrantGhost( DN14a,GHdata,hkt,'BINX',binx );
if save_mem
    clear DN14a;
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


%-------------------------------------------------------------------------%
% apply bad a priori pixel interpolation
% TRRIFdata.readCDR('BP');
% switch EDRdata.lbl.OBSERVATION_TYPE
%     case {'FRT','HRL','HRS','MSP','HSP','FFC'}
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
        [ DN14d,BP ] = apriori_badpixel_removal( DN14c,BPdata1,BPdata2,DMdata,'InterpOpt',1 );
    case 'NONE'
        DN14d = DN14c;
end

%-------------------------------------------------------------------------%
% fourth step (nonlinearity correction)
% LCdata = TRRIFdata.readCDR('LC');
[ DN14g ] = nonlinearity_correction( DN14d,LCdata,hkt,'BINX',binx );
[ DN14g_woc ] = nonlinearity_correction( DN14c,LCdata,hkt,'BINX',binx );
if save_mem
    clear DN14c DN14b;
end

%-------------------------------------------------------------------------%
% fifth step (division by exposure time)
[ RT14g ] = divide_by_integrationTime( DN14g,hkt );
[ RT14g_woc ] = divide_by_integrationTime( DN14g_woc,hkt );
if save_mem
    clear DN14g DN14g_woc;
end

%-------------------------------------------------------------------------%
%%
%-------------------------------------------------------------------------%
% background subtraction
hkt_df1 = DFdata1.readHKT();
hkt_df1c = correctHKTwithHD(hkt_df1,HDdata);
hkt_df1cc = correctHKTwithHK(hkt_df1c,HKdata);
hkt_df2 = DFdata2.readHKT();
hkt_df2c = correctHKTwithHD(hkt_df2,HDdata);
hkt_df2cc = correctHKTwithHK(hkt_df2c,HKdata);
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
% DMdata = TRRIFdata.readCDR('DM');
[RT14h2,dc] = dark_column_subtract(RT14h,DMdata);
[RT14h2_woc,dc] = dark_column_subtract(RT14h_woc,DMdata);
[RT14h2_bk1_o,dc_bk1] = dark_column_subtract(RT14h_bk1_o,DMdata);
[RT14h2_bk2_o,dc_bk2] = dark_column_subtract(RT14h_bk2_o,DMdata);
if save_mem
    clear RT14h RT14h_woc;
end

%-------------------------------------------------------------------------%
% second order light removal
% LLdata = TRRIFdata.readCDR('LL');
[RT14j,K] = subtract_highorderlight(RT14h2,LLdata);
RT14j_woc = RT14h2_woc - K;
if save_mem
    clear RT14h2 RT14h2_woc;
end

%%
%-------------------------------------------------------------------------%
% shutter mirror nonrepeatability correction
% TRRIFdata.readCDR('SP');
% for i=1:length(TRRIFdata.cdr.SP)
%     spdata = TRRIFdata.cdr.SP(i);
%     spdata_prop = getProp_basenameCDR4(spdata.basename);
%     switch upper(spdata_prop.sensor_id)
%         case 'L'
%             SPdata = spdata;
%         case 'S'
%             SPdataVNIR = spdata;
%         otherwise
%             error('sensor_id %s is wrong',sensor_id);
%     end
% end

% SSdata = TRRIFdata.readCDR('SS');
% SHdata = TRRIFdata.readCDR('SH');
if isempty(SPdata_o)
    SPdata_o = SPdata;
    [SR] = calculate_SR(SSdata,SPdata,SHdata,0);
else
    [SR] = calculate_SR(SSdata,SPdata,SHdata,MP);
end

%-------------------------------------------------------------------------%
% calculate spectroradiometric responsitivity
% if isempty(SPdata_o), SPdata_o = SPdata; end
[RSPj] = calculate_RSP(SPdata_o,SR);
rowNumTableRSPj = SPdata.read_ROWNUM_TABLE();

%-------------------------------------------------------------------------%
% apply binning
% DMdata = TRRIFdata.readCDR('DM');
[RSPl] = binning_RSP(RSPj,DMdata,rowNumTableRSPj,'BINX',binx);

%-------------------------------------------------------------------------%
% correct to radiance with the binned responsitivity and flat fielding
% NUdata = TRRIFdata.readCDR('NU');
[RDm,FF] = calculate_RD(RT14j,RSPl,NUdata,'FLAT_FIELD',flat_field);
[RDm_woc,FF_woc] = calculate_RD(RT14j_woc,RSPl,NUdata,'FLAT_FIELD',flat_field);
[RDm_bk1_o,FF_bk1] = calculate_RD(RT14h2_bk1_o,RSPl,NUdata,'FLAT_FIELD',flat_field);
[RDm_bk2_o,FF_bk2] = calculate_RD(RT14h2_bk2_o,RSPl,NUdata,'FLAT_FIELD',flat_field);
if save_mem
    clear RT14j RT14j_woc;
end

RDn = apply_DM(RDm,DMdata);
RDn_woc = apply_DM(RDm_woc,DMdata);
RDn_bk1_o = apply_DM(RDm_bk1_o,DMdata);
RDn_bk2_o = apply_DM(RDm_bk2_o,DMdata);

end