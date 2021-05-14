function [RSPj] = calculate_RSP_VNIR(SPdata,SSdata,SHdata,TDdata,MP,rowNumTableRSPj)
% [RSPj] = calculate_RSP_VNIR(SPdata,SSdata,SHdata,TDdata,MP,rownum_table)
%  calculate the spectralradiometric responsitivity using the sphere 
%  radiance on which shutter mirror non-repeatability correction is applied
%  for vnir.
%   Input
%     SPdata: CRISMdata obj, CDR SP data
%     SSdata: CRISMdata obj, CDR SS data
%     SHdata: CRISMdata obj, CDR SH data
%     TDdata: CRISMdata obj, CDR TD data, not applied so far.
%     MP: scalar, parameter for shutter mirror non-repeatability correction
%     rowNumTableRSPj: rownumtable for RSP
%   Output Parameters
%     RSP: spectral radiometric responsitivity

if isempty(SPdata.img), SPdata.readimg(); end
if isempty(TDdata.img), TDdata.readimg(); end
if isempty(SSdata.img), SSdata.readimg(); end
if isempty(SHdata.img), SHdata.readimg(); end

% propSS8 = SSdata.prop;
% propSS8.version = 8;
% basenameSS8 = crism_get_basenameCDR4_fromProp(propSS8);
% SSdata8 = CRISMdata(basenameSS8,'');
% SSdata8.readimg();

RTS14 = SPdata.img;

SC = SHdata.img(2,:,:);

Tc3 = SPdata.lbl.MRO_SPHERE_TEMPERATURE;

theta = SSdata.img(1,:,:);
sigma = SSdata.img(2,:,:);
tau = SSdata.img(3,:,:);

SR_Tc3 = theta + sigma * Tc3 + tau + tau * Tc3.^2;

SR = SR_Tc3 ./ (1 + MP.* SC);

RSPj = RTS14 ./ SR;

% For detector rows 185-215, SSdata is only applied
[~,blt563] = ismember((185:215)',rowNumTableRSPj);
blt563eq0 = (blt563==0);
% remove any wavelength doesn't match
if any(blt563eq0), blt563 = blt563(blt563eq0); end
RSPj(:,:,blt563) = SR_Tc3(:,:,blt563);

end