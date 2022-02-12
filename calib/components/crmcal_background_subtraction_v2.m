function [ RT14h,Bkgd ] = crmcal_background_subtraction_v2( RT14g,BKdata1,BKdata2,hkt )
% [ RT14h,Bkgd ] = crmcal_background_subtraction_v2( RT14g,BKdata1,BKdata2,hkt )
%  The sixth step of the calibration:  Subtract dark current
%  Input Parameters
%    RT14g  : 14bit DN image (L,S,B) divided by milliseconds
%    BKdata1 : CRISMdata obj, CDR BK data, before the measurement
%    BKdata2 : CRISMdata obj, CDR BK data, after the measurement
%              can be empty
%    hkt    : housekeeping table data (from TRR3)
%  Output parameters
%    RT14h    : processed data [L,S,B]
%    Bkgd     : processed Bkgd [L,S,B]
%
%  2019/11/06: YUKI ITOH: empty BKdata2 is supported.

if isempty(BKdata1.img), BKdata1.readimg(); end
if ~isempty(BKdata2)
    if isempty(BKdata2.img),BKdata2.readimg(); end
end

[L,S,Bands] = size(RT14g);

bk1_integ_t =BKdata1.lbl.MRO_EXPOSURE_PARAMETER;
bk1_rateHz = BKdata1.lbl.MRO_FRAME_RATE.value;
[bk1_exptm] = crism_get_integrationTime(bk1_integ_t,bk1_rateHz,'Hz');
% RTD1 = BKdata1.img ./ bk1_exptm;
RTD1 = BKdata1.img;
t_bk1 = mean([BKdata1.get_sclk_start(),BKdata1.get_sclk_stop()]);

if ~isempty(BKdata2)
    bk2_integ_t =BKdata2.lbl.MRO_EXPOSURE_PARAMETER;
    bk2_rateHz = BKdata2.lbl.MRO_FRAME_RATE.value;
    [bk2_exptm] = crism_get_integrationTime(bk2_integ_t,bk2_rateHz,'Hz');
    % RTD2 = BKdata2.img ./ bk2_exptm;
    RTD2 = BKdata2.img;
    t_bk2 = mean([BKdata2.get_sclk_start(),BKdata2.get_sclk_stop()]);
end

t_scene = [hkt.data.EXPOSURE_SCLK_S]+[hkt.data.EXPOSURE_SCLK_SS]/(2.^16);
t_scene = t_scene(:);

% dt = t_bk2-t_bk1;
% Bkgd = zeros(size(RT14g));
% for l=1:L
%     Bkgdl =  RTD1*(t_bk2-t_scene(l)) + RTD2*(t_scene(l)-t_bk1);
%     Bkgd(l,:,:) = Bkgdl;
% end
% Bkgd = Bkgd ./ dt;
% Bkgd = repmat(RTD1,[L,1,1]); % just use the prior dark measurement

Bkgd = RTD1;

RT14h = RT14g - Bkgd;

end
