function [RDn,RDn_woc,RDn_bk1_o,RDn_bk2_o] = minipipeline_calibration_IR_SP_yuki(...
    SPdata,DFdata1,DFdata2,BKdata1,BKdata2,BPdata1,BPdata2,bkoption,varargin)
%   Mini pipeline for the calibration of the CRISM images
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
dn4095_rmvl = false;
bkgd_robust = false;
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
            case 'DN4095_RMVL'
                dn4095_rmvl = varargin{i+1};
            case 'BKGD_ROBUST'
                bkgd_robust = varargin{i+1};
                if bkoption==1
                    error('no effect of "Bkgd_robust" when bkoption=%d',bkoption);
                end
            otherwise
                % Hmmm, something wrong with the parameter string
                error(['Unrecognized option: ''' varargin{i} '''']);
        end
    end
end

if isempty(SPdata.basenamesCDR), SPdata.load_basenamesCDR(); end
if isempty(SPdata.basenames_SOURCE_OBS)
    SPdata.load_basenames_SOURCE_OBS(); 
end


% get EDRSPdata from SPdata
basenameEDRSP = SPdata.basenames_SOURCE_OBS.SP;
if iscell(basenameEDRSP)
    error('Mulitple EDR SP is used. not supported.');
end
EDRSPdata = CRISMdata(basenameEDRSP,'');

% get DFdata from SPdata
propEDRDF1 = create_propOBSbasename();
propEDRDF1.obs_counter = sprintf('%02s',dec2hex(hex2dec(EDRSPdata.prop.obs_counter)-1));
propEDRDF1.obs_class_type = EDRSPdata.prop.obs_class_type;
propEDRDF1.obs_id = EDRSPdata.prop.obs_id;
propEDRDF1.product_type = EDRSPdata.prop.product_type;
propEDRDF1.sensor_id = EDRSPdata.prop.sensor_id;
propEDRDF1.activity_id = 'DF';

[~,~,~,basenameEDRDF1] = get_dirpath_observation_fromProp(propEDRDF1,'dwld',1);

DFdata1 = CRISMdata(basenameEDRDF1,'');
sclk_df1 = DFdata1.get_sclk_stop();


% get BKdata using sclk of the DFdata



% get BPdata using sclk of the DFdata



frame_rate = SPdata.lbl.MRO_FRAME_RATE{1};
binx = SPdata.lbl.PIXEL_AVERAGING_WIDTH;

DN = EDRdata.readimg();
if dn4095_rmvl
    flg_saturation = (DN==4095); % added by Yuki Feb.18, 2019 
end

% it was expected that saturated pixels are taken as nan after the
% subtraction of background, but it seems that sometimes, it doesn't occur 
% maybe because of some are not saturated.
rownum_table = EDRdata.read_ROWNUM_TABLE();
%-------------------------------------------------------------------------%
% first step (DN12 --> DN14)
PPdata = SPdata.readCDR('PP');
[ DN14 ] = DN12toDN14( DN,PPdata,rownum_table );
if save_mem
    clear DN;
end
%-------------------------------------------------------------------------%
% second step (subtract bias)
SPdata.readCDR('BI');
for i=1:length(SPdata.cdr.BI)
    bidata = SPdata.cdr.BI{i};
    if bidata.lbl.MRO_FRAME_RATE{1} == frame_rate
        BIdata = bidata;
    end
end
BSdata = SPdata.readCDR('BS');
DBdata = SPdata.readCDR('DB');
EBdata = SPdata.readCDR('EB');
HDdata = SPdata.readCDR('HD');
HKdata = SPdata.readCDR('HK');
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
GHdata = SPdata.readCDR('GH');
[ DN14b,sumGhost ] = remove_quadrantGhost( DN14a,GHdata,hkt,'BINX',binx );
if save_mem
    clear DN14a;
end

%-------------------------------------------------------------------------%
% apply bad a priori pixel interpolation
SPdata.readCDR('BP');
switch EDRdata.lbl.OBSERVATION_TYPE
    case {'FRT','HRL','HRS','MSP','HSP'}
        for i=1:length(SPdata.cdr.BP)
            bpdata = SPdata.cdr.BP{i};
            if ~any(strcmpi(EDRdata.basename,bpdata.lbl.SOURCE_PRODUCT_ID))
                if any(strcmpi(DFdata1.basename,bpdata.lbl.SOURCE_PRODUCT_ID))
                    BPdata1 = bpdata;
                elseif any(strcmpi(DFdata2.basename,bpdata.lbl.SOURCE_PRODUCT_ID))
                    BPdata2 = bpdata;
                end
            else
                BPdata_post = bpdata;
            end
        end
    case {'FRS','ATO'}
        % in case of FRS, DFdata1 and DFdata2 are same.
        for i=1:length(SPdata.cdr.BP)
            bpdata = SPdata.cdr.BP{i};
            if ~any(strcmpi(EDRdata.basename,bpdata.lbl.SOURCE_PRODUCT_ID))
                if any(strcmpi(DFdata1.basename,bpdata.lbl.SOURCE_PRODUCT_ID))
                    BPdata1 = bpdata; BPdata2 = bpdata;
                end
            else
                BPdata_post = bpdata;
            end
        end
    otherwise
        error('Undefined observation type %s.',EDRdata.lbl.OBSERVATION_TYPE);
end

DMdata = SPdata.readCDR('DM');
% [ DN14c,BP ] = apriori_badpixel_removal( DN14b,BPdata1,BPdata2,DMdata,'InterpOpt',1 );

switch upper(apbprmvl)
    case 'HIGHORD'
        [ DN14c,BP ] = apriori_badpixel_removal( DN14b,BPdata1,BPdata2,DMdata,'InterpOpt',1 );
    case 'NONE'
        DN14c = DN14b;
end
%-------------------------------------------------------------------------%
% flag saturated pixels
% VLdata = TRRIFdata.readCDR('VL');
% added by Yuki Itoh.
% saturation removal is performed after detector quadrant ghost removal.
% this way is manually defined by Yuki. Doesn't follow the direction in the
% crism_dpsis.pdf
if dn4095_rmvl
    DN14b(flg_saturation) = nan;
end

%-------------------------------------------------------------------------%
% fourth step (nonlinearity correction)
LCdata = SPdata.readCDR('LC');
[ DN14g ] = nonlinearity_correction( DN14c,LCdata,hkt,'BINX',binx );
[ DN14g_woc ] = nonlinearity_correction( DN14b,LCdata,hkt,'BINX',binx );
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


%%
%-------------------------------------------------------------------------%
% process darks from scratch
[~,BKdata1_o] = calculate_Bkgd_wDF(DFdata1,bkgd_robust,PPdata,BSdata,...
    DBdata,EBdata,HDdata,HKdata,BIdata,GHdata,LCdata);
[~,BKdata2_o] = calculate_Bkgd_wDF(DFdata2,bkgd_robust,PPdata,BSdata,...
    DBdata,EBdata,HDdata,HKdata,BIdata,GHdata,LCdata);
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
DMdata = SPdata.readCDR('DM');
[RT14h2,dc] = dark_column_subtract(RT14h,DMdata);
[RT14h2_woc,dc] = dark_column_subtract(RT14h_woc,DMdata);
[RT14h2_bk1_o,dc_bk1] = dark_column_subtract(RT14h_bk1_o,DMdata);
[RT14h2_bk2_o,dc_bk2] = dark_column_subtract(RT14h_bk2_o,DMdata);
if save_mem
    clear RT14h RT14h_woc;
end

%-------------------------------------------------------------------------%
% second order light removal
LLdata = SPdata.readCDR('LL');
[RT14j,K] = subtract_highorderlight(RT14h2,LLdata);
RT14j_woc = RT14h2_woc - K;
if save_mem
    clear RT14h2 RT14h2_woc;
end