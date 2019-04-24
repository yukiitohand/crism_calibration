function [SR] = calculate_SR(SSdata,SPdata,SHdata,MP)
% [SR] = calculate_SR(SSdata,SPdataVNIR,SPdata,SHdata)
%  calculate the sphere radiance on which shutter mirror non-repeatability
%  correction is applied
%   Input
%     SSdata: CRISMdata obj, CDR SS data
%     SPdataVNIR: CRISMdata obj, CDR SP data VNIR (used for obtaining MP)
%     SPdata: CRISMdata obj, CDR SP data
%     SHdata: CRISMdata obj, CDR SH data.
%     MP: scalar, shutter mirror parameter
%   Output Parameters
%     SR: spectral radiance

if isempty(SSdata.img), SSdata.readimg(); end
if isempty(SPdata.img), SPdata.readimg(); end
if isempty(SHdata.img), SHdata.readimg(); end


SC = SHdata.img(2,:,:);

Tc3 = SPdata.lbl.MRO_SPHERE_TEMPERATURE;

theta = SSdata.img(1,:,:);
sigma = SSdata.img(2,:,:);
tau = SSdata.img(3,:,:);

SR_Tc3 = theta + sigma * Tc3 + tau * Tc3.^2;

SR = SR_Tc3 .* (1 + MP.* SC);

end
    
    

