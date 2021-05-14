function [ DN14 ] = crmcal_DN12toDN14_VNIR( DN,PPdata,rownum_table )
% [ DN14 ] = crmcal_DN12toDN14_VNIR( DN,PPdata,rownum_table )
%   convert 12bit DN data to 14bit DNdata for VNIR
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

VNIR_ROW = [PPdata.tab.data.VNIR_ROW];
VNIR_GAIN = [PPdata.tab.data.VNIR_GAIN];
VNIR_OFFSET = [PPdata.tab.data.VNIR_OFFSET];

% match row
idx = arrayfun(@(x) find(x==VNIR_ROW),rownum_table);

vnir_gain = VNIR_GAIN(idx);
vnir_offset = VNIR_OFFSET(idx);

vnir_gain = reshape(vnir_gain,[1,1,B]);
vnir_offset = reshape(vnir_offset,[1,1,B]);

DN14 = repmat(vnir_offset,[L,S,1]) + DN ./ repmat(vnir_gain,[L,S,1]);

end