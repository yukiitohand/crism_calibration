function [RDn,RDn_woc] = pipeline_calibration_IR_original(TRRIFdata,EDRdata,DFdata1,DFdata2,BKdata1,BKdata2,dcoption,varargin)
% perform calibration of CRISM data 
%  Input Parameters
%   TRRIFdata, EDRdata: TRRIFdata is just for the information of CDR and
%                       other information
%   DFdata1,DFdata2, prior and post DF measurements
%   BKdata1,BKdata2, prior and post Background measurements
%   dcoption: {0,1}
%             0: no dark column subtract
%             1: dark column subtract is performed
%             recommended to use 1.
%  Optional Parameters
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


isdebug = 0;
RAdata = [];
save_mem = false;
apbprmvl = 'HighOrd';
if (rem(length(varargin),2)==1)
    error('Optional parameters should always go by pairs');
else
    for i=1:2:(length(varargin)-1)
        switch upper(varargin{i})
            case 'DEBUG'
                isdebug = varargin{i+1};
            case 'RADATA'
                RAdata = varargin{i+1};
            case 'APBPRMVL'
                apbprmvl = varargin{i+1};
            case 'SAVE_MEMORY'
                save_mem = varargin{i+1};
        end
    end
end


frame_rate = TRRIFdata.lbl.MRO_FRAME_RATE{1};
binx = TRRIFdata.lbl.PIXEL_AVERAGING_WIDTH;

DN = EDRdata.readimg();
rownum_table = EDRdata.read_ROWNUM_TABLE();

% first step (DN12 --> DN14)
PPdata = TRRIFdata.readCDR('PP');
[ DN14 ] = DN12toDN14( DN,PPdata,rownum_table );
if save_mem
    clear DN;
end

% second step (subtract bias)
TRRIFdata.readCDR('BI');
for i=1:length(TRRIFdata.cdr.BI)
    bidata = TRRIFdata.cdr.BI{i};
    if bidata.lbl.MRO_FRAME_RATE{1} == frame_rate
        BIdata = bidata;
    end
end
BSdata = TRRIFdata.readCDR('BS');
DBdata = TRRIFdata.readCDR('DB');
EBdata = TRRIFdata.readCDR('EB');
HDdata = TRRIFdata.readCDR('HD');
HKdata = TRRIFdata.readCDR('HK');
TRRIFdata.readHKT(); hkt = TRRIFdata.hkt;
% using Temperature recorded in the label in TRR I/F data 
[ DN14a,BI_m ] = subtract_bias_1( DN14,BIdata,BSdata,DBdata,EBdata,hkt,rownum_table,TRRIFdata.lbl,'BINX',binx);
% [ DN14a,BI_m ] = subtract_bias( DN14,BIdata,BSdata,DBdata,EBdata,hkt,rownum_table);
if save_mem
    clear DN14;
end

% DMdata = TRRIFdata.readCDR('DM');
% switch dcoption
%     case 1
%         [DN14a,dc] = dark_column_subtract(DN14a,DMdata);
% %         [RT14h2_woc,dc] = dark_column_subtract(RT14h_woc,DMdata);
%     case 0
%         DN14a = DN14a;
% %         RT14h2_woc= RT14h_woc;
% end

% the third step (remove detector quadrant electronics ghost)
GHdata = TRRIFdata.readCDR('GH');
[ DN14b,sumGhost ] = remove_quadrantGhost( DN14a,GHdata,hkt,'BINX',binx );
if save_mem
    clear DN14a;
end

%%
% bad pixel removal

% apply bad a priori pixel interpolation
TRRIFdata.readCDR('BP');
for i=1:length(TRRIFdata.cdr.BP)
    bpdata = TRRIFdata.cdr.BP{i};
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
DMdata = TRRIFdata.readCDR('DM');
switch upper(apbprmvl)
    case 'HIGHORD'
        [ DN14c,BP ] = apriori_badpixel_removal( DN14b,BPdata1,BPdata2,DMdata,'InterpOpt',1 );
    case 'NONE'
        DN14c = DN14b;
end
%%
% fourth step (nonlinearity correction)
LCdata = TRRIFdata.readCDR('LC');
[ DN14g ] = nonlinearity_correction( DN14c,LCdata,hkt,'BINX',binx );
[ DN14g_woc ] = nonlinearity_correction( DN14b,LCdata,hkt,'BINX',binx );
if save_mem
    clear DN14c DN14b;
end

% fifth step (division by exposure time)
[ RT14g ] = divide_by_integrationTime( DN14g,hkt );
[ RT14g_woc ] = divide_by_integrationTime( DN14g_woc,hkt );
if save_mem
    clear DN14g DN14g_woc;
end

% sixth step (background subtraction)
% TRRIFdata.readCDR('BK');
% switch EDRdata.lbl.OBSERVATION_TYPE
%     case 'FRT'
%         for i=1:length(TRRIFdata.cdr.BK)
%             bkdata = TRRIFdata.cdr.BK{i};
%             if strcmpi(bkdata.get_obsid, obs_id_scene)
%                 if strcmpi(bkdata.get_obs_number,'06')
%                     BKdata1 = bkdata;
%                 elseif strcmpi(bkdata.get_obs_number,'08')
%                     BKdata2 = bkdata;
%                 end
%             end
%         end
%     otherwise
%         error(' calibration for OBSERVATION TYPE %s is not implemented',EDRdata.lbl.OBSERVATION_TYPE);
% end
% I think it makes sense to interpolate background here, too.
[ RT14h,Bkgd ] = background_subtraction( RT14g,BKdata1,BKdata2,hkt );
[ RT14h_woc,Bkgd ] = background_subtraction( RT14g_woc,BKdata1,BKdata2,hkt );
if save_mem
    clear RT14g RT14g_woc;
end

% I think additional dark column subtract is applied
%-------------------------------------------------------------------------%
% dark column subtract
% I think additional dark column subtract is applied somewhere
DMdata = TRRIFdata.readCDR('DM');
switch dcoption
    case 1
        [RT14h2,dc] = dark_column_subtract(RT14h,DMdata);
        [RT14h2_woc,dc] = dark_column_subtract(RT14h_woc,DMdata);
    case 0
        RT14h2 = RT14h;
%         RT14h2_woc= RT14h_woc;
end
if save_mem
    clear RT14h RT14h_woc;
end

% DMdata = TRRIFdata.readCDR('DM');
% DMdata.readimg();
% [L,S,B] = size(DN14);
% darkMask = double(DMdata.img == 2);
% darkMask(darkMask==0) = nan;
% imgDark = RT14h;
% for l=1:L
%     imgDark(l,:,:) = imgDark(l,:,:) .* darkMask;
% end
% imgDark = reshape(permute(imgDark,[2,3,1]),[S*B,L]);
% imgDark = nanmedian(imgDark,1);
% imgDark = squeeze(imgDark)';
% RT14h2 = RT14h - repmat(imgDark,[1,S,B]);

% second order light removal
LLdata = TRRIFdata.readCDR('LL');
[RT14j,K] = subtract_highorderlight(RT14h2,LLdata);
% [~,K] = subtract_highorderlight(RT14h2_woc,LLdata);
RT14j_woc = RT14h2_woc - K;
if save_mem
    clear RT14h2 RT14h2_woc;
end

%%
% shutter mirror nonrepeatability correction
TRRIFdata.readCDR('SP');
for i=1:length(TRRIFdata.cdr.SP)
    spdata = TRRIFdata.cdr.SP{i};
    spdata_prop = getProp_basenameCDR4(spdata.basename);
    switch upper(spdata_prop.sensor_id)
        case 'L'
            SPdata = spdata;
        case 'S'
            SPdataVNIR = spdata;
        otherwise
            error('sensor_id %s is wrong',sensor_id);
    end
end

[SPdataMP,SSdataMP,SHdataMP] = selectCDR4MP(SPdataVNIR);
[MP] = calculate_MP(SPdataMP,SSdataMP,SHdataMP);
SSdata = TRRIFdata.readCDR('SS');
SHdata = TRRIFdata.readCDR('SH');

% unfortunately MP=0 leads to reasonable result
[SR] = calculate_SR(SSdata,SPdata,SHdata,0);
%[SR] = calculate_SR(SSdata,SPdata,SHdata,MP);

% calculate spectroradiometric responsitivity
[RSPj] = calculate_RSP(SPdata,SR);
rowNumTableRSPj = SPdata.read_ROWNUM_TABLE();

% apply binning
DMdata = TRRIFdata.readCDR('DM');
[RSPl] = binning_RSP(RSPj,DMdata,rowNumTableRSPj);

% correct to radiance with the binned responsitivity and flat fielding
NUdata = TRRIFdata.readCDR('NU');
[RDm,FF] = calculate_RD(RT14j,RSPl,NUdata);
[RDm_woc,FF_woc] = calculate_RD(RT14j_woc,RSPl,NUdata);

if save_mem
    clear RT14j RT14j_woc;
end

RDn = apply_DM(RDm,DMdata);
RDn_woc = apply_DM(RDm_woc,DMdata);

if isdebug && ~isempty(RAdata)
    RAdata.readimg();
    [L,S,B] = size(RAdata.img);
    [RT14j_r] = RD2RT14(RAdata.img,RSPl,NUdata);
    RT14h_r = RT14j_r + K;
%     RT14g1_r = RT14h_r + repmat(dc,[1,S,B]);
    RT14g_r = RT14h_r + Bkgd;
    [DN14g_r] = multiply_by_integrationTime(RT14g_r,hkt);
    BPdata1.readimg(); BPdata2.readimg();
    BP = or(BPdata1.img,BPdata2.img);

    BPdata_post.readimg();
    BP_all = double(BPdata_post.readimg());
    BP_all(BP_all==0) = nan;
    GP_all = double(BPdata_post.img==0);
    GP_all(GP_all==0) = nan;
    BP_pri = double(BP);
    d = DN14g-DN14g_r;
    d_nan = d;
    for i=1:size(d_nan,1)
        d_nan_l = d_nan(i,:,:);
        d_nan_l(isnan(GP_all)) = nan;
        d_nan(i,:,:) = d_nan_l;
    end
    figure; 
    c = 500;
    plot(4:438,squeeze(d(:,c,4:438))'.*squeeze(GP_all(1,c,4:438)),'X');
end

end