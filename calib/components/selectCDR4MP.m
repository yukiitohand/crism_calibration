function [SPdata,SSdata,SHdata] = selectCDR4MP(SPdata_input,varargin)

propSP_input = crism_getProp_basenameCDR4(SPdata_input.basename);

% if the data is not VNIR, then raise error
if ~strcmpi(propSP_input.sensor_id,'S')
    error('The input must be VNIR data');
end
    

if propSP_input.binning==0 && propSP_input.wavelength_filter==0
    % if no binning and wavelength filter is applied, then just use the
    % input SPdata
    SPdata = SPdata_input; propSP = propSP_input;
elseif propSP_input.binning>0 && propSP_input.wavelength_filter==0
    % if binning is perfomed and wavelength filter is 0, then select SP
    % data with binning==0. If there is no such data, then use SPdata_input
    % directly.
    propSP_new = setfield(propSP_input,'binning',0);
    [dir_info,basenameSP_new] = crism_search_cdr_fromProp(propSP_new,varargin{:});
    dirpathSP_new = dir_info.dirfullpath_local;
%     [basenamePtrnSP] = crism_get_basenameCDR4_fromProp(propSP_new);
%     [dir_info,basenameSP] = crism_get_dirpath_cdr(SPdata_input.basename,varargin{:});
%     dirpath_cdr = dir_info.dirfullpath_local;
%     remote_subdir = dir_info.subdir_remote;
%     [basenameSP] = readDownloadBasename_v3(basenamePtrnSP,...
%                                     dirpath_cdr,remote_subdir,varargin{:});
    if ~isempty(basenameSP_new)
        SPdata = CRISMdata(basenameSP_new,dirpathSP_new);
        propSP = crism_getProp_basenameCDR4(SPdata.basename);
    else
        SPdata = SPdata_input;
        propSP = propSP_input;
    end
else
    error('not implemented for wavelength_filter(%d)>0',propSP_input.wavelength_filter);
end

%VNIR SS data (as of 2018.01.15, 8 is the latest version.)
% 2 is the right one.
propSS = crism_create_propCDR4basename('Acro','SS','Binning',propSP.binning,...
    'Wavelength_filter',propSP.wavelength_filter,'Side',propSP.side,...
    'SENSOR_ID',propSP.sensor_id,'Version',2);
% basenameSS = crism_get_basenameCDR4_fromProp(propSS);
[dir_info,basenameSS] = crism_search_cdr_fromProp(propSS,varargin{:});
dirpathSS = dir_info.dirfullpath_local;
% basenameSSPtrn = crism_get_basenameCDR4_fromProp(propSS);
% basenameSS = readDownloadBasename_v3(basenameSSPtrn,dirpathSS,remote_subdirSS);
if ~isempty(dirpathSS)
    % This will throw an error.
    basenameSS
    error('Please fix this');
    SSdata = CRISMdata(basenameSS,dirpathSS);
else
    error('%s cannot be found',basenameSS);
end

%VNIR SH data (as of 2018.01.15, 8 is the latest version.)
propSH = crism_create_propCDR4basename('Acro','SH','Binning',propSP.binning,...
    'Wavelength_filter',propSP.wavelength_filter,'Side',propSP.side,...
    'SENSOR_ID',propSP.sensor_id,'Version',4);
[dir_info,basenameSH] = crism_search_cdr_fromProp(propSH);
dirpathSH = dir_info.dirfullpath_local;
% basenameSHPtrn = crism_get_basenameCDR4_fromProp(propSH);
% basenameSH = readDownloadBasename_v3(basenameSHPtrn,dirpathSH,remote_subdirSH,varargin{:});
if ~isempty(dirpathSH)
    SHdata = CRISMdata(basenameSH,dirpathSH);
else
    error('%s cannot be found',basenameSH);
end

end
