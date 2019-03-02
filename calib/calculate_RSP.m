function [RSP] = calculate_RSP(SPdata,SR)
% [RSP] = calculate_RSP(SPdata,SR)
%  calculate the spectralradiometric responsitivity using the sphere 
%  radiance on which shutter mirror non-repeatability correction is applied
%   Input
%     SPdata: CRISMdata obj, CDR SP data
%     SR    : spectral radiance
%   Output Parameters
%     RSP: spectral radiometric responsitivity

if isempty(SPdata.img), SPdata.readimg(); end;

RTS14 = SPdata.img;

RSP = RTS14 ./ SR;

end