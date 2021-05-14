function [hkt_corrHK] = crism_correctHKTwithHK(hkt_corr,HKdata)

if isempty(HKdata.tab), HKdata.readTAB(); end;
hkt_corrHK = hkt_corr;

for i=1:length(HKdata.tab.data)
    hk_name = HKdata.tab.data(i).HK_NAME;
    switch upper(hk_name)
        case 'IR_SPECTRAL_CAVITY_TEMP'
            hk_name_hkt = 'SPECT_CAVITY_TEMP_IR';
        case 'IR_SPHERE_TEMP'
            hk_name_hkt = 'SPHERE_TEMP_IR';
        case 'DPU_POWERBOARD_TEMP'
            hk_name_hkt = 'DPU_BOARD_TEMP';
        case 'DPU_PLUS5_CURRENT'
            hk_name_hkt = 'DPU_P5_CURRENT';
        case 'DPU_SCAN_MOTOR_CURRENT'
            hk_name_hkt = 'SCAN_MOT_CURR';
        case 'IR_POWERBOARD_TEMP'
            hk_name_hkt = 'IR_BOARD_TEMP';
        case 'IR_HEATER_CURRENT'
            hk_name_hkt = 'HEATER_34_CURRENT';
        case 'VNIR_POWERBOARD_TEMP'
            hk_name_hkt = 'VNIR_BOARD_TEMP';
        case 'VNIR_HEATER_CURRENT'
            hk_name_hkt = 'HEATER_12_CURRENT';
        case 'COOLER_POWERBOARD_TEMP'
            hk_name_hkt = 'COOL_BOARD_TEMP';
        case 'COOLER_CURRENT'
            hk_name_hkt = 'COOL_CURRENT';
        case 'COOLER_HOP_HEATER_CURRENT'
            hk_name_hkt = 'HOP_HEATER_12_CURR';
        case 'IR_SHUTTER_MOTOR_CURRENT'
            hk_name_hkt = 'SHUT_MOTOR_CURR_IR';
        case 'IR_COOLER_TEMP1'
            hk_name_hkt = 'COOLER_TEMP1';
        case 'IR_COOLER_TEMP2'
            hk_name_hkt = 'COOLER_TEMP2';
        case 'IR_COOLER_TEMP3'
            hk_name_hkt = 'COOLER_TEMP3';
        case 'IR_SHUTTER_MOTOR_TEMP'
            hk_name_hkt = 'SHUT_MOTOR_TEMP_IR';
        case 'IR_SPIDER_TEMP'
            hk_name_hkt = 'SPIDER_TEMP_IR';
        case 'IR_OSU_CAVITY_TEMP'
            hk_name_hkt = 'OSU_CAVITY_TEMP_IR';
        case 'IR_SPHERE_LAMP_CURRENT'
            hk_name_hkt = 'IR_SPHERE_CURR';
        case 'IR_HOP_TEMPERATURE'
            hk_name_hkt = 'HOP_TEMP_IR';
        case 'VNIR_SHUTTER_MOTOR_CURRENT'
            hk_name_hkt = 'SHUT_MOTOR_CURR_VNIR';    
        case 'VNIR_TELESCOPE_TEMP'
            hk_name_hkt = 'TELESCOPE_TEMP';
        case 'VNIR_OPTICAL_BENCH_TEMP'
            hk_name_hkt = 'OPTICAL_BENCH_TEMP';
        case 'VNIR_SPIDER_TEMP'
            hk_name_hkt = 'SPIDER_TEMP_VNIR';
        case 'VNIR_SPECTRAL_CAVITY_TEMP'
            hk_name_hkt = 'SPECT_CAVITY_TEMP_VNIR';
        case 'VNIR_OSU_CAVITY_TEMP'
            hk_name_hkt = 'OSU_CAVITY_TEMP_VNIR';
        case 'VNIR_BAFFLE_TEMP'
            hk_name_hkt = 'BAFFLETEMP';
        otherwise
            hk_name_hkt = upper(hk_name);
    end
    y = cat(1,hkt_corr.data.(hk_name_hkt));
    
    a = HKdata.tab.data(i).HK_A;
    b = HKdata.tab.data(i).HK_B;
    c = HKdata.tab.data(i).HK_C;
    d = HKdata.tab.data(i).HK_D;
    
    if d==0
         y_corr = a .* y.^2 + b.*y + c;
    else
        y_corr = a .* y.^2 + b .* y ./ (d - y) + c;
    end
    
    y_corr = num2cell(y_corr);
    [hkt_corrHK.data.(upper(hk_name_hkt))] = y_corr{:};
    
end

end