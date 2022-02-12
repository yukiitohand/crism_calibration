function [ DN14 ] = crmcal_DN12toDN14( DN,PPdata,rownum_table )
% [ DN14 ] = crmcal_DN12toDN14( DN,PPdata,rownum_table )
%   convert 12bit DN data to 14bit DNdata
%  Input Parameters
%   DN     : 12bit DN image (L,S,B)
%   PPdata : CDR PP data, storing gain and offset for conversion
%   rownum_table : ROWNUM_TABLE
%  Output Parameters
%   DN14   : 14bit DN image (L,S,B)
% 
%  DN14 = OFFSET_lambda + DN12/Gain_lambda

if isempty(PPdata.tab), PPdata.readTAB(); end

[L,S,B] = size(DN);

IR_ROW = [PPdata.tab.data.IR_ROW];
IR_GAIN = [PPdata.tab.data.IR_GAIN];
IR_OFFSET = [PPdata.tab.data.IR_OFFSET];

% match row
idx = arrayfun(@(x) find(x==IR_ROW),rownum_table);

ir_gain = IR_GAIN(idx);
ir_offset = IR_OFFSET(idx);

ir_gain = reshape(ir_gain,[1,1,B]);
ir_offset = reshape(ir_offset,[1,1,B]);

% DN14 = repmat(ir_offset,[L,S,1]) + DN ./ repmat(ir_gain,[L,S,1]);

DN14 = ir_offset + DN ./ ir_gain;

end









