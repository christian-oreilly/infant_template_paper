% Script for generating Brainstorm anatomy templates for O'Reilly infant atlases.
% Reference: https://github.com/christian-oreilly/infant_template_paper
% 
% Author: Francois Tadel, 2021

atlas_dir = 'C:\Work\RawData\Test\Oreilly\2021';
output_dir = 'C:\Work\RawData\Test\Oreilly\bst_templates';

% Atlases to import:   OrigName, TemplateName, NAS, LPA, RPA, IC, PC, IH
AtlasNames = {...
    'ANTS2-0Weeks3T',   'Oreilly_0.5m_2021',  [132.00, 176.00, 134.00], [93.00, 130.00, 112.00], [171.00, 134.00, 110.00], [133.00, 128.00, 137.00], [133.00, 109.00, 137.00], [133.00, 118.00, 162.00]; ...
    'ANTS1-0Months3T',  'Oreilly_1m_2021',    [128.00, 176.00, 137.00], [86.00, 133.00, 114.00], [169.00, 136.00, 111.00], [129.00, 128.00, 141.00], [129.00, 109.00, 141.00], [129.00, 118.00, 171.00]; ...
    'ANTS2-0Months3T',  'Oreilly_2m_2021',    [131.00, 179.00, 140.00], [90.00, 130.00, 116.00], [172.00, 132.00, 114.00], [131.00, 129.00, 142.00], [131.00, 109.00, 145.00], [131.00, 114.00, 179.00]; ...
    'ANTS3-0Months3T',  'Oreilly_3m_2021',    [130.00, 187.00, 110.00], [79.00, 126.00, 100.00], [184.00, 128.00, 101.00], [132.00, 135.00, 135.00], [132.00, 115.00, 137.00], [132.68, 130.41, 175.72]; ...
    'ANTS4-5Months3T',  'Oreilly_4.5m_2021',  [129.00, 190.00, 110.00], [76.00, 128.00, 101.00], [181.00, 128.00, 101.00], [128.00, 137.00, 139.00], [128.00, 114.00, 139.00], [128.26, 125.69, 182.14]; ...
    'ANTS6-0Months3T',  'Oreilly_6m_2021',    [125.00, 191.00, 105.00], [76.00, 125.00,  97.00], [177.00, 127.00,  92.00], [129.00, 140.00, 136.00], [129.00, 118.00, 140.00], [129.14, 129.02, 179.91]; ...
    'ANTS7-5Months3T',  'Oreilly_7.5m_2021',  [130.00, 194.00, 111.00], [75.00, 126.00,  97.00], [181.00, 126.00,  97.00], [], [], []; ...
    'ANTS9-0Months3T',  'Oreilly_9m_2021',    [129.00, 194.00, 109.00], [76.00, 126.00, 101.00], [183.00, 126.00, 101.00], [], [], []; ...
    'ANTS10-5Months3T', 'Oreilly_10.5m_2021', [129.00, 197.00, 101.00], [74.00, 124.00, 91.00], [187.00, 124.00, 99.00], [], [], []; ...
    'ANTS12-0Months3T', 'Oreilly_12m_2021',   [126.00, 198.00, 113.00], [72.00, 126.00, 97.00], [184.00, 130.00, 98.00], [], [], []; ...
    'ANTS15-0Months3T', 'Oreilly_15m_2021',   [125.00, 198.00, 114.00], [71.00, 124.00, 95.00], [184.00, 128.00, 96.00], [], [], []; ...
    'ANTS18-0Months3T', 'Oreilly_18m_2021',   [129.00, 197.00, 102.00], [71.00, 122.00, 97.00], [185.00, 123.00, 96.00], [], [], []; ...
    'ANTS2-0Years3T',   'Oreilly_24m_2021',   [126.00, 197.00, 104.00], [70.00, 119.00, 95.00], [185.00, 122.00, 95.00], [], [], []};


% Start brainstorm without the GUI
if ~brainstorm('status')
    brainstorm nogui
end
% Create output folder if doesn't exist
if ~file_exist(output_dir)
    mkdir(output_dir);
end

% Import atlas in database
for iAtlas = 1:size(AtlasNames,1)
    % Delete existing subject
    SubjectName = AtlasNames{iAtlas,1};
    [sSubject, iSubject] = bst_get('Subject', SubjectName);
    if ~isempty(sSubject)
        db_delete_subjects(iSubject);
    end
    
    % FreeSurfer folder
    FsDir = fullfile(atlas_dir, SubjectName);
    
    % Read T1 MRI, to get the world transformation
    T1File = fullfile(FsDir, 'mri', 'T1.mgz');
    sMri = in_mri(T1File, 'MGH', 0);

    % Process: Import FreeSurfer folder
    bst_process('CallProcess', 'process_import_anatomy', [], [], ...
        'subjectname', SubjectName, ...
        'mrifile',     {FsDir, 'FreeSurfer+Thick'}, ...
        'nvertices',   15000, ...
        'nas', AtlasNames{iAtlas,3}, ...
        'lpa', AtlasNames{iAtlas,4}, ...
        'rpa', AtlasNames{iAtlas,5}, ...
        'ac', AtlasNames{iAtlas,6}, ...
        'pc', AtlasNames{iAtlas,7}, ...
        'ih', AtlasNames{iAtlas,8});

    % Get subject
    [sSubject, iSubject] = bst_get('Subject', SubjectName);
    % Delete unwanted files
    AnatDir = bst_fileparts(file_fullpath(sSubject.FileName));
    MidFiles = file_find(AnatDir, 'tess_cortex_mid_*.mat', 1, 0);
    file_delete(MidFiles, 1, 3);
    db_reload_subjects(iSubject);
    
    % Import BEM surfaces
    BemFiles = {fullfile(FsDir, 'bem', 'inner_skull.surf'), ...
                fullfile(FsDir, 'bem', 'outer_skull.surf'), ...
                fullfile(FsDir, 'bem', 'outer_skin.surf')};
    import_surfaces(iSubject, BemFiles, 'FS', 1);
    
    % List all electrodes files
    TsvFiles = file_find(FsDir, '*_electrodes.tsv', 2, 0);
    % Import each channel file
    for iTsv = 1:length(TsvFiles)
        % Split filename
        [fPath, MontageName, fExt] = bst_fileparts(TsvFiles{iTsv});
        MontageName = strrep(MontageName, '_electrodes', '');
        % Create folder
        iStudy = db_add_condition(SubjectName, MontageName);
        % Import channel file
        ChannelFile = import_channel(iStudy, TsvFiles{iTsv}, 'FREESURFER-TSV', 0, 0, 1, 0, 0);
        % Update file
        if ~isempty(ChannelFile)
            % Get full path
            ChannelFile = file_fullpath(ChannelFile);
            % Update file comment
            ChannelMat = load(ChannelFile);
            ChannelMat.Comment = [MontageName, sprintf(' (%d)', length(ChannelMat.Channel))];
            bst_save(ChannelFile, ChannelMat, 'v7');
            % Rename file based on input file name
            [fPath, fBase, fExt] = bst_fileparts(ChannelFile);
            movefile(ChannelFile, fullfile(fPath, ['channel_', file_standardize(MontageName), '.mat']));
            % Reload folder
            db_reload_studies(iStudy);
        end
    end
    
    % Rename subject
    db_rename_subject(SubjectName, AtlasNames{iAtlas,2});
    % Export as anatomy template
    ZipFile = export_default_anat(iSubject, AtlasNames{iAtlas,2}, 1);
    % Move zip file to output folder
    movefile(ZipFile, output_dir);
end
