function [hkt_corr] = crism_correctHKTwithHD(hkt,HDdata)
% 

if isempty(HDdata.tab), HDdata.readTAB(); end;
rate = cat(1,hkt.data.RATE);

hk_names = unique({HDdata.tab.data.HK_NAME});

hkt_corr = hkt;

for i=1:length(hk_names)
    hk_name = hk_names{i};
    if strcmpi(hk_name,'IR_SPECTRAL_CAVITY_TEMP')
        hk_name_hkt = 'SPECT_CAVITY_TEMP_IR';
    elseif strcmpi(hk_name,'IR_SPHERE_TEMP')
        hk_name_hkt = 'SPHERE_TEMP_IR';
    else
        hk_name_hkt = upper(hk_name);
    end
    y = cat(1,hkt.data.(hk_name_hkt));
    tab_data_sub = searchby('HK_NAME',hk_name,HDdata.tab.data,'Comp_func','strcmpi');
    tab_sub = [];
    tab_sub.data = tab_data_sub;
    tab_sub.colinfo = HDdata.tab.colinfo;
    
    % VNIR_SPHERE_COEFFICIENT (col 3)
    % Add this coefficient if hkt.data.VNIR_SPHERE_PWR
    c3 = crism_frame_rateTABformatter(rate,tab_sub,'VNIR_SPHERE_COEFFICIENT',...
                                'COLUMNS',1);
    vnir_sphere_pwr = cat(1,hkt.data.VNIR_SPHERE_PWR);
    y3 = y + vnir_sphere_pwr .* c3;

    % VNIR_FLOOD1_COEFFICIENT (col 4)
    % Add hkt.data.VNIR_Flood1_Level / (this coefficient)
    % if hkt.data.VNIR_FLOOD1_PWR
    c4 = crism_frame_rateTABformatter(rate,tab_sub,'VNIR_FLOOD1_COEFFICIENT',...
                                'COLUMNS',1);
    c4(c4==0) = inf;
    vnir_flood1_pwr = cat(1,hkt.data.VNIR_FLOOD1_PWR);
    vnir_flood1_level = cat(1,hkt.data.VNIR_FLOOD1_LEVEL);
    y4 = y3 + vnir_flood1_pwr .* (vnir_flood1_level ./ c4);

    % VNIR_FLOOD2_COEFFICIENT (col 5)
    % Add hkt.data.VNIR_Flood2_Level / (this coefficient)
    % if hkt.data.VNIR_FLOOD2_PWR
    c5 = crism_frame_rateTABformatter(rate,tab_sub,'VNIR_FLOOD2_COEFFICIENT',...
                                'COLUMNS',1);
    c5(c5==0) = inf;
    vnir_flood2_pwr = cat(1,hkt.data.VNIR_FLOOD2_PWR);
    vnir_flood2_level = cat(1,hkt.data.VNIR_FLOOD2_LEVEL);
    y5 = y4 + vnir_flood2_pwr .* (vnir_flood2_level ./ c5);

    % IR_SPHERE_COEFFICIENT (col 6)
    % Add this coefficient if hkt.data.IR_SPHERE_PWR
    c6 = crism_frame_rateTABformatter(rate,tab_sub,'IR_SPHERE_COEFFICIENT',...
                                'COLUMNS',1);
    
    ir_sphere_pwr = cat(1,hkt.data.IR_SPHERE_PWR);
    y6 = y5 + ir_sphere_pwr .* c6;

    % IR_FLOOD1_COEFFICIENT (col 7)
    % Add hkt.data.IR_Flood1_Level / (this coefficient)
    % if hkt.data.IR_FLOOD1_PWR
    c7 = crism_frame_rateTABformatter(rate,tab_sub,'IR_FLOOD1_COEFFICIENT',...
                                'COLUMNS',1);
    c7(c7==0) = inf;
    ir_flood1_pwr = cat(1,hkt.data.IR_FLOOD1_PWR);
    ir_flood1_level = cat(1,hkt.data.IR_FLOOD1_LEVEL);
    y7 = y6 + ir_flood1_pwr .* (ir_flood1_level ./ c7);

    % IR_FLOOD2_COEFFICIENT (col 8)
    % Add hkt.data.IR_Flood2_Level / (this coefficient)
    % if hkt.data.IR_FLOOD2_PWR
    c8 = crism_frame_rateTABformatter(rate,tab_sub,'IR_FLOOD2_COEFFICIENT',...
                                'COLUMNS',1);
    c8(c8==0) = inf;
    ir_flood2_pwr = cat(1,hkt.data.IR_FLOOD2_PWR);
    ir_flood2_level = cat(1,hkt.data.IR_FLOOD2_LEVEL);
    y8 = y7 + ir_flood2_pwr .* (ir_flood2_level ./ c8);

    % COOLER1_LEVEL_COEFFICIENT (col 9)
    % Add hkt.data.IR_COOL1_LEVEL / (this coefficient)
    % if hkt.data.IR_COOL1_PWR
    c9 = crism_frame_rateTABformatter(rate,tab_sub,'COOLER1_LEVEL_COEFFICIENT',...
                                'COLUMNS',1);
    c9(c9==0) = inf;
    ir_cool1_pwr = cat(1,hkt.data.IR_COOL1_PWR);
    ir_cool1_level = cat(1,hkt.data.IR_COOL1_LEVEL);
    y9 = y8 + ir_cool1_pwr .* (ir_cool1_level ./ c9);

    % COOLER2_LEVEL_COEFFICIENT (col 10)
    % Add hkt.data.IR_COOL2_LEVEL / (this coefficient)
    % if hkt.data.IR_COOL2_PWR
    c10 = crism_frame_rateTABformatter(rate,tab_sub,'COOLER2_LEVEL_COEFFICIENT',...
                                'COLUMNS',1);
    c10(c10==0) = inf;
    ir_cool2_pwr = cat(1,hkt.data.IR_COOL2_PWR);
    ir_cool2_level = cat(1,hkt.data.IR_COOL2_LEVEL);
    y10 = y9 + ir_cool2_pwr .* (ir_cool2_level ./ c10);

    % COOLER3_LEVEL_COEFFICIENT (col 11)
    % Add hkt.data.IR_COOL3_LEVEL / (this coefficient)
    % if hkt.data.IR_COOL3_PWR
    c11 = crism_frame_rateTABformatter(rate,tab_sub,'COOLER3_LEVEL_COEFFICIENT',...
                                'COLUMNS',1);
    c11(c11==0) = inf;
    ir_cool3_pwr = cat(1,hkt.data.IR_COOL3_PWR);
    ir_cool3_level = cat(1,hkt.data.IR_COOL3_LEVEL);
    y11 = y10 + ir_cool3_pwr .* (ir_cool3_level ./ c11);

    % VNIR_EXPOSURE_COEFFICIENT (col 12)
    % Add (1-480) hkt.data.VNIR_EXPOSE / (this coefficient)
    c12 = crism_frame_rateTABformatter(rate,tab_sub,'VNIR_EXPOSURE_COEFFICIENT',...
                                'COLUMNS',1);
    c12(c12==0) = inf;
    vnir_expose = cat(1,hkt.data.VNIR_EXPOSE);
    y12 = y11 + (vnir_expose ./ c12);

    % IR_EXPOSURE_COEFFICIENT (col 13)
    % Add (1-480) hkt.data.IR_EXPOSE / (this coefficient)
    c13 = crism_frame_rateTABformatter(rate,tab_sub,'IR_EXPOSURE_COEFFICIENT',...
                                'COLUMNS',1);
    c13(c13==0) = inf;
    ir_expose = cat(1,hkt.data.IR_EXPOSE);
    y13 = y12 + (ir_expose ./ c13);

    % FRAMERATE_COMPENSATION_OFFSET (col 14)
    % Add this offset to counts
    offset = crism_frame_rateTABformatter(rate,tab_sub,'FRAMERATE_COMPENSATION_OFFSET',...
                                'COLUMNS',1);
    y14 = y13 + offset;

    % FRAMERATE_COMPENSATION_SLOPE (col 15)
    % Multiply counts by this coefficient
    slope = crism_frame_rateTABformatter(rate,tab_sub,'FRAMERATE_COMPENSATION_SLOPE',...
                                'COLUMNS',1);
    y15 = y14 .* slope;

    % FRAMERATE_COMPENSATION_CURV (col 16)
    % Multiply counts^2 by this coefficient and add to counts
    curv_cff = crism_frame_rateTABformatter(rate,tab_sub,'FRAMERATE_COMPENSATION_CURV',...
                                'COLUMNS',1);
%     y16 = y15 + curv_cff .* y15.^2;
    y16 = offset + y13 .* slope + curv_cff .* y13.^2;
    y_corr = y16;
%     y_corr = offset + y.*slope + y.^2 .* curv_cff ...
%              + vnir_sphere_pwr .* c3 ...
%              + vnir_flood1_pwr .* (vnir_flood1_level ./ c4) ...
%              + vnir_flood2_pwr .* (vnir_flood2_level ./ c5) ...
%              + ir_sphere_pwr .* c6 ...
%              + ir_flood1_pwr .* (ir_flood1_level ./ c7) ...
%              + ir_flood2_pwr .* (ir_flood2_level ./ c8) ...
%              + ir_cool1_pwr .* (ir_cool1_level ./ c9) ...
%              + ir_cool2_pwr .* (ir_cool2_level ./ c10) ...
%              + ir_cool3_pwr .* (ir_cool3_level ./ c11) ...
%              + (vnir_expose ./ c12) ...
%              + (ir_expose ./ c13);

    y_corr = num2cell(y_corr);
    [hkt_corr.data.(upper(hk_name_hkt))] = y_corr{:};
    
end

end






