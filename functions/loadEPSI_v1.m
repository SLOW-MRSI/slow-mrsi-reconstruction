function  [data, Info,twix_obj,mdh]  = loadEPSI_v1(fileName, dataPath, rawDataPath, opt)
% Load data from EPSI_fid sequence using the mapVbVd function from Siemens
% Created by Dr. Chao Ma @mgh, 2018-07-13
% SLOW-MRSI package maintenance: Dr. Guodong Weng, University of Bern, 2026-07-04
%

    mdh = [];

    if exist([dataPath,fileName],'file') && opt.reload == 0	% if the file exist 
        load([dataPath,fileName],'twix_obj','mdh');
    else 									% The .mat file not specified, read raw data
        if isfield(opt,'measID')
            measID                      = opt.measID;
        else
            measID                      = [];
            if isfield(opt,'path_to_data')
                path_to_data            = opt.path_to_data;
            else
                [fileNameRaw, pathName] = uigetfile({[rawDataPath,'*.dat']},...
                    'Select the .dat data file to be read','MultiSelect','on');
                if isequal(fileNameRaw,0)
                    disp('User pressed cancel')
                    return;
                end
                path_to_data            = [pathName,fileNameRaw];
            end
        end

        % obtain the twix_obj
        if ~isempty(measID)
            tpath                       = pwd;
            cd(rawDataPath);
            [twix_obj]                  = mapVBVD(measID); % only keep the first output
            cd(tpath);
        else
            [twix_obj,mdh]              = mapVBVD(path_to_data);
        end

        % save data
%       save([dataPath,fileName],'twix_obj','mdh');
    end

	% obtain Info
    if(length(twix_obj)>1)
        twix_obj                        = twix_obj{end};
    end
	Info.params                         = twix_obj.hdr.Meas;
	Info.Nc                             = twix_obj.image.NCha;
    Info.Nx                             = Info.params.lBaseResolution;
	Info.Ny                             = Info.params.lPhaseEncodingLines;
	Info.Nz                             = Info.params.lPartitions;
    if(isfield(opt,'Nave'))
        Info.Nave                       = opt.Nave;
    else
        Info.Nave                       = twix_obj.image.NAve;
    end
    if(isfield(opt,'Nrep'))
        Info.Nrep                       = opt.Nrep;
    else
        Info.Nrep                       = twix_obj.image.NRep;
    end
    if(~isfield(opt,'B0Nav'))
        opt.B0Nav                       = 0;%twix_obj.hdr.MeasYaps.sWipMemBlock.alFree{9};
        if (isempty(opt.B0Nav))
            opt.B0Nav                   = 0;
        end
    end

	% obtain unsorted data
	data                                = twix_obj.image.unsorted();
    
    % data_1                              = data(:,:,1:2:end); % ful
    % data_2                              = data(:,:,2:2:end); % par
    % 
    % Info.Nt                             = size(data_1,1);
    % tNPE                                = size(data_1,3)./Info.Nave/Info.Nrep;
    % data_1                                = reshape(data_1,[Info.Nt,Info.Nc,Info.Nave,tNPE,Info.Nrep]); 
    % data_2                                = reshape(data_2,[Info.Nt,Info.Nc,Info.Nave,tNPE,Info.Nrep]); 
    % 
    % if(opt.B0Nav>0)
    %     ind_B0Nav                       = 1:(opt.B0Nav+1):tNPE;
    %     if ind_B0Nav(end)<tNPE
    %         ind_B0Nav                   = [ind_B0Nav,tNPE];
    %     end
    %     ind_image                       = [];
    %     for tind = 2:length(ind_B0Nav)
    %         ind_image                   = [ind_image,(ind_B0Nav(tind-1)+1:ind_B0Nav(tind)-1)];
    %     end
    %     Info.B0NavData                  = data_1(:,:,:,ind_B0Nav,:);
    %     NPE                             = length(ind_image);
    % else
    %     ind_image                       = 1:tNPE;
    %     NPE                             = tNPE;
    % end
    % 
	% % k-t space dataPath    
	% Info.PreCutoff                      = mdh.sCutOff(1,1);
	% Info.PostCutoff                     = mdh.sCutOff(1,2);
    % if NPE == Info.Nx*Info.Ny*Info.Nz
    %     ktData                          = reshape(data_1(:,:,:,ind_image,:),[Info.Nt,Info.Nc,Info.Nave,...
    %                                                                        Info.Nx,Info.Ny,Info.Nz,Info.Nrep]);
    % else
    %     if strcmp(twix_obj.image.softwareVersion,'vb')
    %         temp                        = reshape(twix_obj.image.Phs(:),[Info.Nave,tNPE,Info.Nrep]);
    %         PE1                         = reshape(temp(1,ind_image,:),[NPE,Info.Nrep]);
    %         temp                        = reshape(twix_obj.image.Lin(:),[Info.Nave,tNPE,Info.Nrep]);
    %         PE2                         = reshape(temp(1,ind_image,:),[NPE,Info.Nrep]);
    %         temp                        = reshape(twix_obj.image.Seg(:),[Info.Nave,tNPE,Info.Nrep]);
    %         PE3                         = reshape(temp(1,ind_image,:),[NPE,Info.Nrep]);
    %     else
    %         temp                        = reshape(twix_obj.image.Lin(1:2:end),[Info.Nave,tNPE,Info.Nrep]); % GW
    %         PE1                         = reshape(temp(1,ind_image,:),[NPE,Info.Nrep]);
    %         temp                        = reshape(twix_obj.image.Seg(1:2:end),[Info.Nave,tNPE,Info.Nrep]);
    %         PE2                         = reshape(temp(1,ind_image,:),[NPE,Info.Nrep]);
    %         temp                        = reshape(twix_obj.image.Par(1:2:end),[Info.Nave,tNPE,Info.Nrep]);
    %         PE3                         = reshape(temp(1,ind_image,:),[NPE,Info.Nrep]);
    %     end
    %     ktData_1                          = zeros([Info.Nt,Info.Nc,Info.Nave,Info.Nx,Info.Ny,Info.Nz,Info.Nrep]);
    %     ktData_2                          = zeros([Info.Nt,Info.Nc,Info.Nave,Info.Nx,Info.Ny,Info.Nz,Info.Nrep]);
    %     Info.PE1                        = PE1;
    %     Info.PE2                        = PE2;
    %     Info.PE3                        = PE3;
    %     data_1                          = data_1(:,:,:,ind_image,:);
    %     data_2                          = data_2(:,:,:,ind_image,:);
    %     for indRep = 1:Info.Nrep
    %         for indPE = 1:NPE
    %             ktData_1(:,:,:,PE1(indPE,indRep),PE2(indPE,indRep),PE3(indPE,indRep),:) = data_1(:,:,:,indPE,:); 
    %             ktData_2(:,:,:,PE1(indPE,indRep),PE2(indPE,indRep),PE3(indPE,indRep),:) = data_2(:,:,:,indPE,:); 
    %             Info.kyz_mask(PE1(indPE,indRep),PE2(indPE,indRep),PE3(indPE,indRep),indRep) = true;
    %         end
    %     end
    % end
    % clear data;
    % 
    % % re-arrange data
    % if(length(size(ktData_1))==5)
    %     if (Info.Nz>1)
    % 
    %     else
    %         ktData_1                      = permute(ktData_1,[4,5,1,2,3]);
    %         ktData_2                      = permute(ktData_2,[4,5,1,2,3]);
    %     end
    % elseif (length(size(ktData_1))==6)
    %     if (Info.Nz>1)
    %         ktData_1                      = permute(ktData_1,[4,5,6,1,2,3]);    
    %         ktData_2                      = permute(ktData_2,[4,5,6,1,2,3]);  
    %     end
    % end
    % 
    % ktData_1                              = flipdim(ktData_1,2); % flip data in the readout direction
    % ktData_2                              = flipdim(ktData_2,2); % flip data in the readout direction
end
