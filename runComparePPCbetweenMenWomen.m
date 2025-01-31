% This GUI helps compare PPC vs freq between men and women for a
% given protocol
function runComparePPCbetweenMenWomen()
 
    % Create the GUI
    fig = figure('Name', 'PPC comparision between men and women', 'NumberTitle', 'off', ...
                 'Position', [50, 100, 800, 600]);

    % Panel for Controlpanel
    controlPanel = uipanel('Parent', fig, 'Position', [0, 0.92, 1, 0.15], ...
                           'Title', 'Controls');

    % Dropdown for Subject Group
    uicontrol('Parent', controlPanel, 'Style', 'text', 'Position', [10, 20, 120, 20], 'String', 'Subject Group:');
    groupDropdown = uicontrol('Parent', controlPanel, 'Style', 'popupmenu', 'Position', [130, 20, 120, 20], ...
                               'String', {'Meditators', 'Controls'});

    % Dropdown for Protocol
    uicontrol('Parent', controlPanel, 'Style', 'text', 'Position', [260, 20, 120, 20], 'String', 'Protocol :');
    protocolDropdown = uicontrol('Parent', controlPanel, 'Style', 'popupmenu', 'Position', [350, 20, 120, 20], ...
                                  'String', {'EO1', 'EC1', 'G1', 'M1', 'G2', 'EO2', 'EC2', 'M2' });
        % Dropdown for Analysis Choice
    uicontrol('Parent', controlPanel, 'Style', 'text', 'Position', [480, 20, 120, 20], 'String', 'Analysis Choice:');
    analysisDropdown = uicontrol('Parent', controlPanel, 'Style', 'popupmenu', 'Position', [600, 20, 120, 20], ...
                                  'String', {'Baseline', 'Stimulus'});

    % Dropdown for Reference Electrode Group
    uicontrol('Parent', controlPanel, 'Style', 'text', 'Position', [750, 20, 190, 20], 'String', 'Reference Electrode Group:');
    referenceDropdown = uicontrol('Parent', controlPanel, 'Style', 'popupmenu', 'Position', [950, 20, 120, 20], ...
                                   'String', {'Left Occipital','Right Occipital','Back Occipital','Left Frontal', 'Right Frontal'});
    % Dropdown for FreqRange
    uicontrol('Parent', controlPanel, 'Style', 'text', 'Position', [1100, 20, 190, 20], 'String', 'Freq Range:');
    uicontrol('Parent', controlPanel, 'Style', 'edit', 'Position', [1250, 20, 50, 20], 'String', '50', ...
              'Tag', 'freqMin');
    uicontrol('Parent', controlPanel, 'Style', 'edit', 'Position', [1300, 20, 50, 20], 'String', '60', ...
              'Tag', 'freqMax');
    % Plot Button
    uicontrol('Parent', controlPanel, 'Style', 'pushbutton', 'Position', [1400, 32, 80, 20], 'String', 'Plot', ...
              'Callback', @(~, ~) plotConnectivity());

% plot handles _______________________________________
        hConnPPC = getPlotHandles(1,8,[0.05 0.4 0.6 0.45], 0.02, 0.01, 1);
        hTopoRef = getPlotHandles(1,2,[0.66 0.55 0.32 0.3],0.002,0.002,1);
        hTopo  = getPlotHandles(1,3,[0.1 0.05 0.5 0.3],0.02,0,1);
        hConnDistance = getPlotHandles(1,1,[0.7 0.1 0.2 0.25],0.01,0.01,1);
      % Clear Button____________________________________
    uicontrol('Parent', controlPanel, 'Style', 'pushbutton', 'Position', [1400, 12, 80, 20], 'String', 'Clear', ...
              'Callback', @(~, ~) cla_Callback());

    % Plot Connectivity Function
    function plotConnectivity()
        % Get selected group, protocol, and analysis choice
        subjectsChoice = groupDropdown.Value;
        protocolName = protocolDropdown.Value;
          analysisChoice = analysisDropdown.Value;
        refElectrodes = referenceDropdown.Value;
        freqMin = str2double(findobj('Tag', 'freqMin').String);
        freqMax = str2double(findobj('Tag', 'freqMax').String);
        freqRange = [freqMin, freqMax];

        % Map selections to actual strings
        subjectGroups = {'med', 'con'};
        protocols = {'EO1', 'EC1', 'G1', 'M1', 'G2', 'EO2', 'EC2', 'M2' };
        refElectrodeGroups = {[14 44 47], [19 49 52], [16 17 18 48], [34 36 3],[62 63 30]};
        freqRanges = {[60 90], [15 25], [7 10]};
         analyses = {'bl', 'st'};

        subjectsChoice = subjectGroups{subjectsChoice};
        protocolName = protocols{protocolName};
        analysisChoice = analyses{analysisChoice};
        refElectrodes = refElectrodeGroups{refElectrodes};

      
        %refElectrodes = [16 17 18 48]; % back occipital
        groupType = 'rel';
        cutoffList = [2 25];
        displaySignificanceFlag = 1;
        pairedSubjectNameList = getPairedSubjectsBK1;            
        subjectNameLists{1} = pairedSubjectNameList(:,1); % med
        subjectNameLists{2} = pairedSubjectNameList(:,2); % con
        pairedDataFlag      = 0;
        % Sub-select Subjects based on Demographics
        [subjectNameList,~,~,ageListAllSub,genderListAllSub] = getDemographicDetails('BK1');
        % sub-select based on Gender
       
        maleSubjectNameList = subjectNameList(strcmpi(genderListAllSub, 'M'));
        femaleSubjectNameList = subjectNameList(strcmpi(genderListAllSub, 'F'));
       
            subjectNameListsM{1} = intersect(subjectNameLists{1},maleSubjectNameList,'stable'); % men meditators
            subjectNameListsM{2} = intersect(subjectNameLists{2},maleSubjectNameList,'stable'); % men controls
  
            subjectNameListsW{1} = intersect(subjectNameLists{1},femaleSubjectNameList,'stable'); % women med
            subjectNameListsW{2} = intersect(subjectNameLists{2},femaleSubjectNameList,'stable');
   

if ~exist('refElectrodes','var');         refElectrodes = [];           end
if ~exist('groupType','var');             groupType='rel';              end
if ~exist('connMethod','var');            connMethod = 'ppc';           end
if ~exist('badEyeCondition','var');       badEyeCondition='ep';         end
if ~exist('badTrialVersion','var');       badTrialVersion='v8';         end

if ~exist('freqRangeList','var')
    freqRangeList{1} = freqRange; % takes from GUI
%     freqRangeList{2} = [20 34]; % SG
%     freqRangeList{3} = [35 65]; % FG
end
numFreqRanges = length(freqRangeList);
freqRangeColors = copper(numFreqRanges);

if ~exist('axisRangeList','var')
    axisRangeList{1} = [0 100]; % freqLims
    axisRangeList{2} = [0 1]; % yLims for 
    axisRangeList{3} = [0 1]; % cLimsTopo
end
if ~exist('cutoffList','var')
    cutoffList = [2 30]; % elcs & trials
end

if ~exist('useMedianFlag','var');         useMedianFlag = 0;            end
if ~exist('pairedDataFlag','var');        pairedDataFlag = 1;           end
if ~exist('displayDataFlag','var');       displayDataFlag = 1;          end

saveFolderName = 'savedData';
capType = 'actiCap64_UOL';


[connDataM,connDataW, freqVals,connDataElectrodeGroupM,connDataElectrodeGroupW, electrodeGroupList,groupNameList,binnedCenters] = getConnDataAllSubjects(subjectNameListsM,subjectNameListsW, refElectrodes,groupType,connMethod,badEyeCondition,badTrialVersion,protocolName, analysisChoice,cutoffList,pairedDataFlag,saveFolderName,capType);

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Display options %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
displaySettings.fontSizeLarge = 10;
displaySettings.tickLengthMedium = [0.025 0];

titleStr{1} = 'Men ';
titleStr{2} = 'Women ';
titleStr{3} = 'W-M';
freqLims = axisRangeList{1};
yLims = axisRangeList{2};
cLimsTopo = axisRangeList{3};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Show Electrodes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if displayDataFlag
    % Display electrode positions
    if strcmp(groupType,'rel')  
        montageChanlocs = showElectrodeGroups(hTopoRef,capType,electrodeGroupList(1,:),groupNameList); % Just show for the first electrode
    else
        montageChanlocs = showElectrodeGroups(hTopoRef,capType,electrodeGroupList,groupNameList);
    end
end
%%%%%%%%%%%%%%%%%%%%%%% Get frequency positions %%%%%%%%%%%%%%%%%%%%%%%%%%%
freqPosList = cell(1,numFreqRanges);
lineNoiseRange = [46 54];
lineNoiseFreqVals = 46:54;
badFreqPos = intersect(find(freqVals>=lineNoiseRange(1)),find(freqVals<=lineNoiseRange(2)));
for i = 1:numFreqRanges
    freqPosList{i} = setdiff(intersect(find(freqVals>=freqRangeList{i}(1)),find(freqVals<=freqRangeList{i}(2))),badFreqPos);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%% Show Results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numElectrodeGroups = length(electrodeGroupList);

for m=1:numElectrodeGroups
    if strcmp(subjectsChoice, 'med')
    dataToPlot{1} = squeeze(connDataElectrodeGroupM{1}(:,m,:)); % med, St
    dataToPlot{2} = squeeze(connDataElectrodeGroupW{1}(:,m,:)); % med , Bl
    else
    dataToPlot{1} = squeeze(connDataElectrodeGroupM{2}(:,m,:)); % con, st
    dataToPlot{2} = squeeze(connDataElectrodeGroupW{2}(:,m,:)); % con , bl
    end
    
    if displayDataFlag
        % Connectivity as a function of frequency__________________
   displaySettings.colorNames(1,:) = [0 1 0]; % green for St
   displaySettings.colorNames(2,:) =  [1 0 0]; % red for Bl
   displaySettings.colorNames(3,:) =  [0 0 0]; % black for diff

        displayAndcompareData(hConnPPC(m), dataToPlot,freqVals,displaySettings,yLims,displaySignificanceFlag,useMedianFlag,~pairedDataFlag);
        sgtitle(["PPC of " + subjectsChoice + " during " + protocolName + " -- baseline and stimulus wrt elecs " + num2str(refElectrodes)]);

        title(groupNameList{m});
        xlim(hConnPPC(m),freqLims);
        if m==1
        ylabel('PPC'); xlabel('frequency');
        end
        for j=1:2
            text(10,yLims(2)-0.1*j,[titleStr{j} '(' num2str(size(dataToPlot{j},1)) ')'],'color',displaySettings.colorNames(j,:),'parent',hConnPPC(m));
            
        end
    end
end
 % Display topoplots______________________
topoplotDataToReturn = cell(2,numFreqRanges);
    if strcmp(subjectsChoice, 'med')
            xtemp{1} = connDataM{1}; % 
            xtemp{2} = connDataW{1}; %
            xtemp{3} = xtemp{1} -xtemp{2}; % diff
    else 
            xtemp{1} = connDataM{2}; % M
            xtemp{2} = connDataW{2}; % Women
            xtemp{3} = xtemp{2} -xtemp{1}; % W-M
    end
for i=1:numFreqRanges

    dataToCompare = cell(1,2);
    % Topoplots_________

    for j=1:3 % 
            x = squeeze(mean(xtemp{j}(:,:,freqPosList{i}),3, 'omitnan'));
            if useMedianFlag
                data = squeeze(median(x,1,'omitnan'));
            else
                data = squeeze(mean(x,1,'omitnan'));
            end
        
        topoplotDataToReturn{i,j} = data;

        if displayDataFlag
                axes(hTopo(i,j)); %#ok<*LAXES>
                colormap('jet');
                topoplot(data,montageChanlocs,'electrodes','on','plotrad',0.6,'headrad',0.6); colorbar;
                %(.., 'maplimits',cLimsTopo,..)
                if i==1
                title(titleStr{j},'color',displaySettings.colorNames(j,:));
                end
        end

        % Conn versus distance_______________________
        %x = squeeze(mean(connDataElectrodeGroupSt{j}(:,:,freqPosList{i}),3));
          if strcmp(subjectsChoice, 'med')
            dataToCompare{1} = squeeze(mean(connDataElectrodeGroupM{1}(:,:,freqPosList{i}),3));
            dataToCompare{2} = squeeze(mean(connDataElectrodeGroupW{1}(:,:,freqPosList{i}),3));
          else
            dataToCompare{1} = squeeze(mean(connDataElectrodeGroupM{2}(:,:,freqPosList{i}),3));
            dataToCompare{2} = squeeze(mean(connDataElectrodeGroupW{2}(:,:,freqPosList{i}),3));
          end
        
   end
%     flipdataToCompare = flip(dataToCompare);
    displayAndcompareData(hConnDistance(i), dataToCompare, binnedCenters, displaySettings,yLims,1,useMedianFlag,~pairedDataFlag);
    ylabel('PPC'); xlabel('cos\theta');
        if i==1
        title('PPC vs cos\theta');
       end
end
end 

% __________________________________________________________________________
 function [connDataM,connDataW, freqVals,connDataElectrodeGroupM,connDataElectrodeGroupW, electrodeGroupList,groupNameList,binnedCenters] = getConnDataAllSubjects(subjectNameListsM,subjectNameListsW, refElectrodes,groupType,connMethod,badEyeCondition,badTrialVersion,protocolName, analysisChoice,cutoffList,pairedDataFlag,saveFolderName,capType)

[electrodeGroupList,groupNameList,binnedCenters] = getElectrodeGroupsConn(groupType,refElectrodes,capType);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize
badSubjectList = cell(1,2); % common list for both protocol 1 & 2
connDataTMPm = cell(1,2);
connDataTMPw = cell(1,2);
connDataElectrodeGroupTMPm = cell(1,2);
connDataElectrodeGroupTMPw = cell(1,2);
numRefElectrodes = length(refElectrodes);
numElectrodeGroups = length(electrodeGroupList);

for i=1:2 % med and con
    for j=1:length(subjectNameListsW{i}) % no. of subjects in 'i'
        subjectNameM = subjectNameListsM{i}{j};
        subjectNameW = subjectNameListsW{i}{j};
        tmpDataM = load(fullfile(saveFolderName,subjectNameM,[protocolName '_' badEyeCondition '_' badTrialVersion '_' connMethod '.mat']));
        tmpDataW = load(fullfile(saveFolderName,subjectNameW,[protocolName '_' badEyeCondition '_' badTrialVersion '_' connMethod '.mat']));
        numGoodTrialsM = tmpDataM.numGoodTrials;
        numGoodTrialsW = tmpDataW.numGoodTrials;
        if numGoodTrialsM < cutoffList(2) || numGoodTrialsW < cutoffList(2)
            badSubjectList{i}(j) = 1;
        else
            if strcmp(analysisChoice,'bl')
                connDataTMPtmpM = tmpDataM.connPre(refElectrodes,:,:);
                connDataTMPtmpW = tmpDataW.connPre(refElectrodes,:,:);
                freqVals = tmpDataM.freqPre;
            elseif strcmp(analysisChoice,'st')
                connDataTMPtmpM = tmpDataM.connPost(refElectrodes,:,:);
                connDataTMPtmpW = tmpDataW.connPost(refElectrodes,:,:);
                freqVals = tmpDataM.freqPost;
            else
                connDataTMPtmpM = (tmpDataM.connPre(refElectrodes,:,:) + tmpDataM.connPost(refElectrodes,:,:))/2;
                connDataTMPtmpW = (tmpDataW.connPre(refElectrodes,:,:) + tmpDataW.connPost(refElectrodes,:,:))/2;
                freqVals = tmpDataM.freqPost;
            end

            numGoodElectodesM = trace(~isnan(squeeze(connDataTMPtmpM(:,refElectrodes,1))));
            numGoodElectodesW = trace(~isnan(squeeze(connDataTMPtmpW(:,refElectrodes,1))));

            if numGoodElectodesM >= cutoffList(1) && numGoodElectodesW >= cutoffList(1)
                badSubjectList{i}(j) = 0;


                connDataTMPm{i}{j} = squeeze(mean(connDataTMPtmpM,1,'omitnan'));
                connDataTMPw{i}{j} = squeeze(mean(connDataTMPtmpW,1,'omitnan'));

                connDataElectrodeGroupTMPtmpM = zeros(numRefElectrodes,numElectrodeGroups,length(freqVals));
                connDataElectrodeGroupTMPtmpW = zeros(numRefElectrodes,numElectrodeGroups,length(freqVals));
                for e=1:numRefElectrodes
                    for b=1:numElectrodeGroups
                        if strcmp(groupType,'rel')
                            connDataElectrodeGroupTMPtmpM(e,b,:) = squeeze(mean(connDataTMPtmpM(e,electrodeGroupList{e,b},:),2,'omitnan'));
                            connDataElectrodeGroupTMPtmpW(e,b,:) = squeeze(mean(connDataTMPtmpW(e,electrodeGroupList{e,b},:),2,'omitnan'));
                        else
                            connDataElectrodeGroupTMPtmpM(e,b,:) = squeeze(mean(connDataTMPtmpM(e,electrodeGroupList{b},:),2,'omitnan'));
                             connDataElectrodeGroupTMPtmpW(e,b,:) = squeeze(mean(connDataTMPtmpW(e,electrodeGroupList{b},:),2,'omitnan'));
                        end
                    end
                end
                connDataElectrodeGroupTMPm{i}{j} = squeeze(mean(connDataElectrodeGroupTMPtmpM,1,'omitnan'));
                connDataElectrodeGroupTMPw{i}{j} = squeeze(mean(connDataElectrodeGroupTMPtmpW,1,'omitnan'));
                else
                    badSubjectList{i}(j) = 1;
            end
        end
    end
end
    % if paired dataset, take common bad subjects
    if pairedDataFlag
        badSubjectList{1} = union(find(badSubjectList{1}),find(badSubjectList{2}));
        badSubjectList{2} = badSubjectList{1};
    else 
        for o=1:2
            badSubjectList{o}= find(badSubjectList{o});
        end
    end

% Remove bad subjects
connDataM = cell(1,3);
connDataW = cell(1,3);
connDataElectrodeGroupM = cell(1,3);
connDataElectrodeGroupW = cell(1,3);

for q=1:2
x1 = connDataTMPm{q};
x2 = connDataTMPw{q};
y1 = connDataElectrodeGroupTMPm{q};
y2 = connDataElectrodeGroupTMPw{q};
badSubjectPos = badSubjectList{q};
badSubjectPos = badSubjectPos(badSubjectPos<=length(connDataTMPm{q})); 
x1(badSubjectPos)=[];
x2(badSubjectPos)=[];
y1(badSubjectPos)=[]; 
y2(badSubjectPos)=[];  

numSubjects = length(x1);

    for p=1:numSubjects
        connDataM{q}(p,:,:) = x1{p};
        connDataW{q}(p,:,:) = x2{p};
        connDataElectrodeGroupM{q}(p,:,:) = y1{p};
        connDataElectrodeGroupW{q}(p,:,:) = y2{p};
    end
end
% connDataM{3} = connDataM{1} - connDataM{2}; % diff: med-con
% connDataW{3} = connDataW{1} - connDataW{2};
% connDataElectrodeGroupM{3} = connDataElectrodeGroupM{1} - connDataElectrodeGroupM{2};
% connDataElectrodeGroupW{3} = connDataElectrodeGroupW{1} - connDataElectrodeGroupW{2};
 end
function cla_Callback(~,~)
        claGivenPlotHandle(hConnPPC);
        claGivenPlotHandle(hTopoRef);
        claGivenPlotHandle(hTopo);
        claGivenPlotHandle(hConnDistance);

        function claGivenPlotHandle(plotHandles)
            [numRows,numCols] = size(plotHandles);
            for ii=1:numRows
                for j=1:numCols
                    cla(plotHandles(ii,j));
                end
            end
        end
    end
end




