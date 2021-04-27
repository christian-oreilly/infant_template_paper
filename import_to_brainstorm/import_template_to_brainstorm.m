
% Create a new protocol
ProtocolName='infant_templates';
gui_brainstorm('DeleteProtocol', ProtocolName);
gui_brainstorm('CreateProtocol', ProtocolName, 0, 0);

% Get EEG defaults folder
RefChannelPath = bst_fullfile(bst_get('BrainstormDefaultsDir'), 'eeg', 'ICBM152');

study   = bst_get('Study');
ProtocolSubjects = bst_get('ProtocolSubjects');   
ProtocolInfo = bst_get('ProtocolInfo');

files = dir("../Templates/ANTS*3T");

montages = {'10-5', '10-10', '10-20', 'HGSN128'};

userDir = char(java.lang.System.getProperty('user.home'));
delete [userDir '/.brainstorm/defaults/anatomy/*']

for noSubject=1:1 %size(files, 1)
    
    subjectName = files(noSubject).name;    
    FsDir = [files(noSubject).folder '/' subjectName];        
               
    nbSubjects = length(ProtocolSubjects.Subject);
    iSubject = noSubject;
    sSubject = db_template('Subject');
    sSubject.Name     = subjectName;
    sSubject.FileName = bst_fullfile(subjectName, 'brainstormsubject.mat');
    
    % Get subjects defaults from protocol
    sSubject.UseDefaultAnat    = ProtocolInfo.UseDefaultAnat;
    sSubject.UseDefaultChannel = ProtocolInfo.UseDefaultChannel;    
    
    sSubject = db_add_subject(sSubject, iSubject);
   

    [BstMriFile, sMri] = import_mri(noSubject, [FsDir '/mri/T1.mgz']);

    % Read the fiducials
    fiducialFile = [FsDir '/fiducials.tsv'];
    fiducials = tdfread(fiducialFile, '\t');    
    fiducials_world = [fiducials.x  fiducials.y  fiducials.z];      
    fiducials_mri = cs_convert(sMri, 'world', 'mri', fiducials_world ./ 1000) .* 1000;
    fiducials_mri = fiducials_mri + repmat(transpose(sMri.Header.Pxyz_c), ...
                                           size(fiducials_mri, 1), 1);
           
    % in m
	sFid = struct('NAS', fiducials_mri(3, :), ...
                  'LPA', fiducials_mri(5, :), ...
                  'RPA', fiducials_mri(6, :), ...
                  'AC',  fiducials_mri(1, :), ...
                  'PC',  fiducials_mri(2, :), ...
                  'IH',  fiducials_mri(15, :));  

    import_anatomy_fs(iSubject, FsDir, 15000, 0, sFid);
    
    %sMri = bst_get('MriFile', ProtocolSubjects.Subject(iSubject).Anatomy.FileName);

    % Import head surface
    BemFiles = {...
        bst_fullfile(FsDir, 'bem', 'inner_skull.surf'), ...
        bst_fullfile(FsDir, 'bem', 'outer_skull.surf'), ...
        bst_fullfile(FsDir, 'bem', 'outer_skin.surf')};
    [iNewSurfaces, BstBemFiles] = import_surfaces(iSubject, BemFiles, 'FS', 0);    

    % Create folder for channel file
    for iMontage = 1:5
        switch iMontage
            case 1
                MontageName = '10-5';
                RefChannelFile = 'channel_ASA_10-05_343.mat';
            case 2
                MontageName = '10-10';
                RefChannelFile = 'channel_10-10_65.mat';
            case 3
                MontageName = '10-20';
                RefChannelFile = 'channel_10-20_19.mat';
            case 4
                MontageName = 'HGSN128';
                RefChannelFile = 'channel_GSN_HydroCel_128_E1.mat';
            case 5
                MontageName = 'HGSN128';
                RefChannelFile = 'channel_GSN_HydroCel_128_E001.mat';
        end
        % Load reference channel file
		ChannelMat = in_bst_channel(bst_fullfile(RefChannelPath, RefChannelFile));
		FolderName = file_standardize(str_remove_parenth(ChannelMat.Comment));

		% Create new folder
		iStudy = db_add_condition(subjectName, FolderName);
		sStudy = bst_get('Study', iStudy);

        % in mm
        ChannelMat.SCS.LPA = sFid.LPA; % / 1000.0; 
        ChannelMat.SCS.NAS = sFid.NAS; % / 1000.0;
        ChannelMat.SCS.RPA = sFid.RPA; % / 1000.0;            
        
        %ChannelMat.SCS.LPA = sFid.LPA / 1000.0; 
        %ChannelMat.SCS.NAS = sFid.NAS / 1000.0;
        %ChannelMat.SCS.RPA = sFid.RPA / 1000.0;            
                
        
        % Read channel file
        montage = [MontageName '_electrodes.tsv'];
        channels = tdfread([FsDir '/montages/' montage], '\t');  
        channels_world = [channels.x channels.y channels.z];  

        channels_mri = cs_convert(sMri, 'world', 'mri', channels_world ./ 1000) .* 1000;
        channels_mri = channels_mri + repmat(transpose(sMri.Header.Pxyz_c), ...
                                               size(channels_mri, 1), 1);        
        
        for iChannel=1:size(channels_mri, 1)
            if iMontage >= 4
                assert(str2num(ChannelMat.Channel(iChannel).Name(2:end)) == ...
                       str2num(channels.name(iChannel, 2:end)));   
                ChannelMat.Channel(iChannel).Loc = transpose(channels_mri(iChannel, :));
            else                
                ChannelMat.Channel(iChannel) = ChannelMat.Channel(1);
                ChannelMat.Channel(iChannel).Loc = transpose(channels_mri(iChannel, :));
                ChannelMat.Channel(iChannel).Name = strtrim(channels.name(iChannel, :));
            end
        end

        % Setting the head points
        N = size(channels_mri, 1) + size(fiducials.name, 1);
        ChannelMat.HeadPoints = struct('Loc', [zeros(3, N)], 'Label', ...
                                       {cell(1, N)}, 'Type', {cell(1, N)});       
        iHeadPoint = 1;
        for iFid=1:size(fiducials.name, 1)
            label = strtrim(fiducials.name(iFid, :));
            ChannelMat.HeadPoints.Label(iHeadPoint) = {label};
            if strcmp(label, 'NAS') ||  strcmp(label, 'LPA') ||  strcmp(label, 'RPA')
                ChannelMat.HeadPoints.Type(iHeadPoint) = {'CARDINAL'};
            else
                ChannelMat.HeadPoints.Type(iHeadPoint) = {'EXTRA'};                
            end
            Loc =  transpose(fiducials_mri(iFid, :));
            ChannelMat.HeadPoints.Loc(:, iHeadPoint) = Loc;
            iHeadPoint = iHeadPoint+1;
        end
        
        for iCh=1:size(channels.name, 1)
            label = strtrim(channels.name(iCh, :));
            ChannelMat.HeadPoints.Label(iHeadPoint) = {label};
            ChannelMat.HeadPoints.Type(iHeadPoint) = {'EEG'};
            Loc =  transpose(channels_mri(iCh, :));
            ChannelMat.HeadPoints.Loc(:, iHeadPoint) = Loc;
            iHeadPoint = iHeadPoint+1;
        end        
       
        % Re-Register based on NAS/LPA/RPA
        %ChannelMat = channel_detect_type(ChannelMat, 1, 0);    
        
        
        
        
        
        
        % Compute transformation
        transfSCS = cs_compute(ChannelMat, 'scs');
        ChannelMat.SCS.R      = transfSCS.R;
        ChannelMat.SCS.T      = transfSCS.T;
        ChannelMat.SCS.Origin = transfSCS.Origin;
        % Convert the fiducials positions
        %   NOTE: The division/multiplication by 1000 is to compensate the T/1000 applied in the cs_convert().
        %         This hack was added becaue cs_convert() is intended to work on sMri structures, 
        %         in which NAS/LPA/RPA/T fields are in millimeters, while in ChannelMat they are in meters.
        ChannelMat.SCS.NAS = cs_convert(ChannelMat, 'mri', 'scs', ChannelMat.SCS.NAS ./ 1000) .* 1000;
        ChannelMat.SCS.LPA = cs_convert(ChannelMat, 'mri', 'scs', ChannelMat.SCS.LPA ./ 1000) .* 1000;
        ChannelMat.SCS.RPA = cs_convert(ChannelMat, 'mri', 'scs', ChannelMat.SCS.RPA ./ 1000) .* 1000;
        % Process each sensor
        for i = 1:length(ChannelMat.Channel)
            % Converts the electrodes locations to SCS (subject coordinates system)
            ChannelMat.Channel(i).Loc = cs_convert(ChannelMat, 'mri', 'scs', ChannelMat.Channel(i).Loc')' ./ 1000;
        end

        ChannelMat.HeadPoints.Loc = cs_convert(ChannelMat, 'mri', 'scs', ChannelMat.HeadPoints.Loc')'  ./ 1000;

        % Add to the list of transformation
        ChannelMat.TransfMeg{end+1} = [ChannelMat.SCS.R, ChannelMat.SCS.T; 0 0 0 1];
        ChannelMat.TransfMegLabels{end+1} = 'Native=>Brainstorm/CTF';
        ChannelMat.TransfEeg{end+1} = [ChannelMat.SCS.R, ChannelMat.SCS.T; 0 0 0 1];
        ChannelMat.TransfEegLabels{end+1} = 'Native=>Brainstorm/CTF';        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        

        % Save channel file
        ChannelFile = bst_fullfile(bst_fileparts(file_fullpath(sStudy.FileName)), RefChannelFile);
        bst_save(ChannelFile, ChannelMat, 'v7');
		db_reload_studies(iStudy);
    end    

    bst_memory('UnloadAll', 'Forced');
    db_reload_subjects(iSubject);
    panel_protocols('UpdateTree');
    db_save();    
    
    export_default_anat(iSubject, subjectName, 1);    
end

