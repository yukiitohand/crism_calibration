function [ DN14a,BI_m ] = subtract_bias_VNIR_1( DN14,BIdata1,BIdata2,DBdata,EBdata,hkt,rownum_table,lbl_TRR3,varargin )
% [ DN14a ] = subtract_bias_1( DN14,BIdata,BSdata,DBdata,EBdata )
%  This is for replicating calibration
%  The second step of the calibration: subtract Bias for VNIR image
%  Input parameters:
%    DN14   : 14bit DN image (L,S,B)
%    BIdata1 : CRISMdata obj, CDR BI data derived from the prior dark
%              measurement
%    BIdata2 : CRISMdata obj, CDR BI data derived from the post dark
%              measurement
%    DBdata : CRISMdata obj, CDR DB data
%    EBdata : CRISMdata obj, CDR EB data
%    hkt    : housekeeping table data (from TRR3)
%    rownum_table : ROWNUM_TABLE
%  Output parameters
%    DN14a    : processed DN data
%    BI_m   : calculated bias
%  Optional Parameters
%    'BINNING' : ID of binning mode {0,1,2,3}
%                (default) 0
%    'BINX'    : binning size (PIXEL_AVERAGING_WIDTH in LBL)
%                (default) 1 
%     *Note: The relationship between BINNING and BINX
%        BINNING   BINX
%              0      1
%              1      2
%              2      5
%              3     10
%      Please specify only one of 'BINNING', and 'BINX'
%
%  * Detail *
%
%                   ( alpha_ea * BI(1) + alpha_ad * BI(2) ) 
%  DN14a = DN14 - -------------------------------------------
%                             alpha_ea + alpha_ad
%
%                    ( alpha_ea * ( TaV - TdV ) + alpha_ad * ( TaV - TeV ) 
%        - beta_v * -------------------------------------------------------
%                                      alpha_ea + alpha_ad
%
%                    ( alpha_ea * ( TaW - TdW ) + alpha_ad * ( TaW - TeW ) 
%        - beta_w * -------------------------------------------------------
%                                      alpha_ea + alpha_ad
%
%   where
%       DN: 
%           DN value of the image
%       BI: 
%           CDR BI data, BI(1) and BI(2) are derived from the prior and
%           post dark measurements, respectively
%       alpha_ad =  t_a1 - t_d1:
%           t_a1: time_tag (approax. sclk) of the scene/sphere measurement,
%                 is this mean (so constant for each image frame)?? I
%                 assume this is defined different for each frame.
%           t_d1: time_tag (approax. sclk) of the prior dark measurement.
%       alpha_ea =  t_el - t_al:
%           t_e1: time_tag (approax. sclk) of the post dark measurement.
%       So, alpha_ea + alpha_ad = t_e1-t_d1
%       beta_v: 
%           coefficient (CDR DB)
%       T*V:
%           MRO_DETECTOR_TEMPERATURE in lbl ?
%           VNIR_DETECTOR_TEMP1, VNIR_DETECTOR_TEMP2 in HKT ?
%           VNIR detector temperature for the image
%           * can be {a,d,e} where "a", "d", and "e" indicate scene, prior
%           dark, and post dark, respectively.
%           TaV <- hkt
%           TdV and TeV <- CDR lbl
%       beta_w:
%           coefficient (CDR EB)
%       T*W:
%           MRO_FPE_TEMPERATURE in lbl ?
%           VNIR_FPU_BOARD_TEMP in HKT ?
%           VNIR focal plane board temperature for the image
%           * can be {a,d,e} where "a", "d", and "e" indicate scene, prior
%           dark, and post dark, respectively.
%           TaW <- hkt
%           TdW and TeW <- CDR lbl
%
%       
% if isempty(EDRdata.img),EDRdata.readimg(); end
% if isempty(EDRdata.hkt),EDRdata.readHKT(); end
if isempty(BIdata1.img),BIdata1.readimg(); end
if isempty(BIdata2.img),BIdata2.readimg(); end
% if isempty(BSdata.tab),BSdata.readTAB(); end
if isempty(DBdata.tab),DBdata.readTAB(); end
if isempty(EBdata.tab),EBdata.readTAB(); end

binx = 1;
if (rem(length(varargin),2)==1)
    error('Optional parameters should always go by pairs');
else
    for i=1:2:(length(varargin)-1)
        switch upper(varargin{i})
            case 'BINNING'
                binning = varargin{i+1};
                binx = get_binning(binning);
            case 'BINX'
                binx = varargin{i+1};
            otherwise
                % Hmmm, something wrong with the parameter string
                error(['Unrecognized option: ''' varargin{i} '''']);   
        end
    end
end

[L,S,B] = size(DN14);
rate_id = cat(1,hkt.data.RATE);

[t_a1] = get_frame_sclk_mean_fromHKT(hkt);
% t_a1 = 
% lazy for time stamps of CDR BI data
t_d1 = mean([BIdata1.get_sclk_start(),BIdata1.get_sclk_stop()]);
t_e1 = mean([BIdata2.get_sclk_start(),BIdata2.get_sclk_stop()]);

alpha_ea = abs(t_e1 - t_a1);
alpha_ad = abs(t_a1 - t_d1); % supports FRS for same input to BIdata1 and BIdata2
den = t_e1 - t_d1;

term2 = (alpha_ea .* BIdata1.img + alpha_ad .* BIdata2.img) ./ den;


beta_v = rateQuadrantTABformatter(rate_id,DBdata.tab,'A','BINX',binx);
%TaV = lbl_TRR3.MRO_DETECTOR_TEMPERATURE;
TdV = BIdata1.lbl.MRO_DETECTOR_TEMPERATURE;
TeV = BIdata2.lbl.MRO_DETECTOR_TEMPERATURE;
TaV = cat(1,hkt.data.VNIR_DETECTOR_TEMP1);
%TdV = nanmean(cat(1,BIdata1.hkt.data.VNIR_DETECTOR_TEMP1));
%TeV = nanmean(cat(1,BIdata2.hkt.data.VNIR_DETECTOR_TEMP1));


% term3 = beta_v.* (alpha_ea.*(TaV - TdV) + alpha_ad.*(TeV - TaV)) / den;
term3 = beta_v.* (alpha_ea.*(TaV - TdV) + alpha_ad.*(TaV - TeV)) / den;

beta_w = rateQuadrantTABformatter(rate_id,EBdata.tab,'A','BINX',binx);
% TaW = lbl_TRR3.MRO_FPE_TEMPERATURE;
TdW = BIdata1.lbl.MRO_FPE_TEMPERATURE;
TeW = BIdata2.lbl.MRO_FPE_TEMPERATURE;
TaW = cat(1,hkt.data.VNIR_FPU_BOARD_TEMP);
% TdW = nanmean(cat(1,BIdata1.hkt.data.VNIR_FPU_BOARD_TEMP));
% TeW = nanmean(cat(1,BIdata2.hkt.data.VNIR_FPU_BOARD_TEMP));


% term4 = beta_w.* (alpha_ea.*(TaW - TdW) + alpha_ad.*(TeW - TaW)) / den;
term4 = beta_w.* (alpha_ea.*(TaW - TdW) + alpha_ad.*(TaW - TeW)) / den;

BI_m = term2 + term3 + term4;

DN14a = DN14 - BI_m;

end