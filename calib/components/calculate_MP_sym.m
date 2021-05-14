function [MP] = calculate_MP_sym(SPdataVNIR,SSdataVNIR,SHdataVNIR)
% [MP] = calculate_MP_v2(SSdata,SPdataVNIR,SHdata)
%  calculate mirror parameter (MP) for Shutter mirror nonrepeatability
%  correction. Symmetric version, not exactly improving from calculate_MP,
%  not sure this is correct or not.
%   Input
%     SPdataVNIR: CRISMdata obj, CDR SP data VNIR
%     SSdataVNIR: CRISMdata obj, corresponding to the CDR SS data used
%                   for computing MP
%     SHdataVNIR: CRISMdata obj, corresponding to the CDR SH data used
%                   for computing MP
%   Output Parameters
%     MP: shutter mirror parameter, scalar
%   * Detail *
%    See Section 3.4 of "CRISM_DPSIS_Appendix_L_v5_2016-09-01.pdf"
%    VNIR data is used, no matter what


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
propSH = crism_getProp_basenameCDR4(SHdataVNIR.basename);
propSP = crism_getProp_basenameCDR4(SPdataVNIR.basename);
propSS = crism_getProp_basenameCDR4(SSdataVNIR.basename);

if isempty(SSdataVNIR.img), SSdataVNIR.readimg(); end;
if isempty(SPdataVNIR.img), SPdataVNIR.readimg(); end;
if isempty(SHdataVNIR.img), SHdataVNIR.readimg(); end;

DGS14 = SHdataVNIR.img(1,:,:);
exptm = get_integrationTime(SHdataVNIR.lbl.MRO_EXPOSURE_PARAMETER,...
                            SHdataVNIR.lbl.MRO_FRAME_RATE{1},'Hz');
%RGS14 = DGS14/exptm*1000;
RGS14 = DGS14;
Tg3 = SHdataVNIR.lbl.MRO_SPHERE_TEMPERATURE;
SHdataVNIR.read_ROWNUM_TABLE();

SPdataVNIR.read_ROWNUM_TABLE();
RTS14 = SPdataVNIR.img;
Tc3 = SPdataVNIR.lbl.MRO_SPHERE_TEMPERATURE;

SSdataVNIR.read_ROWNUM_TABLE();
theta = SSdataVNIR.img(1,:,:);
sigma = SSdataVNIR.img(2,:,:);
tau = SSdataVNIR.img(3,:,:);
SR_Tc3 = theta + sigma * Tc3 + tau + tau * Tc3.^2;
SR_Tg3 = theta + sigma * Tg3 + tau + tau * Tg3.^2;

SR_isnan = isnan(SR_Tc3);

RGS14(SR_isnan) = nan;
RTS14(SR_isnan) = nan;


switch propSP.side
    case 1
        rownum_A = 223:225;
        rownum_B = 232:234;    
    case 2
        rownum_A = 245:247;
        rownum_B = 235:237;
    otherwise
        error('invalid side %s',side);
end

idxA_RTS = arrayfun(@(x) find(SPdataVNIR.ROWNUM_TABLE==x),rownum_A);
idxB_RTS = arrayfun(@(x) find(SPdataVNIR.ROWNUM_TABLE==x),rownum_B);
idxA_RGS = arrayfun(@(x) find(SHdataVNIR.ROWNUM_TABLE==x),rownum_A);
idxB_RGS = arrayfun(@(x) find(SHdataVNIR.ROWNUM_TABLE==x),rownum_B);
idxA_SS = arrayfun(@(x) find(SSdataVNIR.ROWNUM_TABLE==x),rownum_A);
idxB_SS = arrayfun(@(x) find(SSdataVNIR.ROWNUM_TABLE==x),rownum_B);
% RGSA = RGS14(1,:,idxA_RGS).*SR_Tc3(1,:,idxA_SS)./SR_Tg3(1,:,idxA_SS);
% RGSA = squeeze(nansum(nansum(RGSA,2),3));
% RGSB = RGS14(1,:,idxB_RGS).*SR_Tc3(1,:,idxB_SS)./SR_Tg3(1,:,idxB_SS);
% RGSB = squeeze(nansum(nansum(RGSB,2),3));
% MP = squeeze(nansum(nansum(RTS14(1,:,idxB_RTS),3),2)) ./ RGSB ...
%               - squeeze(nansum(nansum(RTS14(1,:,idxA_RTS),3),2)) ./ RGSA;

RGSA = RGS14(1,:,idxA_RGS)./SR_Tg3(1,:,idxA_SS);
RGSA = squeeze(nansum(nansum(RGSA,2),3));
RTSA = RTS14(1,:,idxA_RTS)./SR_Tc3(1,:,idxA_SS);
RTSA = squeeze(nansum(nansum(RTSA,2),3));
RGSB = RGS14(1,:,idxB_RGS)./SR_Tg3(1,:,idxB_SS);
RGSB = squeeze(nansum(nansum(RGSB,2),3));
RTSB = RTS14(1,:,idxB_RTS)./SR_Tc3(1,:,idxB_SS);
RTSB = squeeze(nansum(nansum(RTSB,2),3));
MP = RTSB/RGSB - RTSA/RGSA;
% MP = squeeze(nansum(nansum(RTS14(1,:,idxB_RTS),3),2)) ./ RGSB ...
%               - squeeze(nansum(nansum(RTS14(1,:,idxA_RTS),3),2)) ./ RGSA;

end
