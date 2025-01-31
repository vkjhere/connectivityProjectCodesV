% This GUI helps compare PPC vs freq between 2 protocols.
function runComparePPCbetweenProtocols()
    % Create the GUI
    fig = figure('Name', 'PPC comparision between two protocols', 'NumberTitle', 'off', ...
                 'Position', [50, 100, 800, 600]);

    % Panel for Controlpanel
    controlPanel = uipanel('Parent', fig, 'Position', [0, 0.92, 1, 0.15], ...
                           'Title', 'Controls');

    % Dropdown for Subject Group
    uicontrol('Parent', controlPanel, 'Style', 'text', 'Position', [10, 20, 120, 20], 'String', 'Subject Group:');
    groupDropdown = uicontrol('Parent', controlPanel, 'Style', 'popupmenu', 'Position', [130, 20, 120, 20], ...
                               'String', {'Meditators', 'Controls'});

    % Dropdown for Protocol1
    uicontrol('Parent', controlPanel, 'Style', 'text', 'Position', [260, 30, 120, 20], 'String', 'Protocol 1:');
    protocol1Dropdown = uicontrol('Parent', controlPanel, 'Style', 'popupmenu', 'Position', [380, 30, 120, 20], ...
                                  'String', {'EO1', 'EC1', 'G1', 'M1', 'G2', 'EO2', 'EC2', 'M2' });

    % Dropdown for Protocol2
    uicontrol('Parent', controlPanel, 'Style', 'text', 'Position', [260, 10, 120, 20], 'String', 'Protocol 2:');
    protocol2Dropdown = uicontrol('Parent', controlPanel, 'Style', 'popupmenu', 'Position',[380, 10, 120, 20], ...
                                  'String', {'EO1', 'EC1', 'G1', 'M1', 'G2', 'EO2', 'EC2', 'M2' });
    % Dropdown for Analysis Choice
    uicontrol('Parent', controlPanel, 'Style', 'text', 'Position', [520, 20, 120, 20], 'String', 'Analysis Choice:');
    analysisDropdown = uicontrol('Parent', controlPanel, 'Style', 'popupmenu', 'Position', [640, 20, 120, 20], ...
                                  'String', {'Baseline', 'Stimulus'});
    % Dropdown for Reference Electrode Group
    uicontrol('Parent', controlPanel, 'Style', 'text', 'Position', [780, 20, 190, 20], 'String', 'Reference Electrode Group:');
    referenceDropdown = uicontrol('Parent', controlPanel, 'Style', 'popupmenu', 'Position', [980, 20, 120, 20], ...
                                   'String', {'Left Occipital','Right Occipital','Back Occipital','Left Frontal', 'Right Frontal'});
  % Dropdown for FreqRange
    uicontrol('Parent', controlPanel, 'Style', 'text', 'Position', [1100, 20, 190, 20], 'String', 'Freq Range:');
    uicontrol('Parent', controlPanel, 'Style', 'edit', 'Position', [1250, 20, 50, 20], 'String', '50', ...
              'Tag', 'freqMin');
    uicontrol('Parent', controlPanel, 'Style', 'edit', 'Position', [1300, 20, 50, 20], 'String', '60', ...
              'Tag', 'freqMax');
    % Plot Button
    uicontrol('Parent', controlPanel, 'Style', 'pushbutton', 'Position', [1400, 10, 80, 20], 'String', 'Plot', ...
              'Callback', @(~, ~) plotConnectivity());
    % plot handles _______________________________________
        hConnPPC = getPlotHandles(1,8,[0.05 0.4 0.57 0.45], 0.02, 0.01, 1);
        hTopoRef = getPlotHandles(1,2,[0.66 0.55 0.32 0.3],0.002,0.002,1);
        hTopo  = getPlotHandles(1,3,[0.1 0.05 0.5 0.3],0.02,0,1);
        hConnDistance = getPlotHandles(1,1,[0.7 0.1 0.2 0.25],0.01,0.01,1);
      % Clear Button___________________________________
    uicontrol('Parent', controlPanel, 'Style', 'pushbutton', 'Position', [1400, 30, 80, 20], 'String', 'Clear', ...
              'Callback', @(~, ~) cla_Callback());

    % Plot Connectivity Function
    function plotConnectivity()
        % Get selected group, protocol, and analysis choice
        subjectsChoice = groupDropdown.Value;
        protocolName1 = protocol1Dropdown.Value;
        protocolName2 = protocol2Dropdown.Value;
        analysisChoice = analysisDropdown.Value;
        refElectrodes = referenceDropdown.Value;
        freqMin = str2double(findobj('Tag', 'freqMin').String);
        freqMax = str2double(findobj('Tag', 'freqMax').String);
        freqRange = [freqMin, freqMax];
        % Map selections to actual strings
        subjectGroups = {'med', 'con'};
        protocols = {'EO1', 'EC1', 'G1', 'M1', 'G2', 'EO2', 'EC2', 'M2' };
        analyses = {'bl', 'st'};
        refElectrodeGroups = {[14 44 47], [19 49 52], [16 17 18 48], [34 36 3],[62 63 30]};
        freqRanges = {[60 90], [15 25],[7 10]};

        subjectsChoice = subjectGroups{subjectsChoice};
        protocolName1 = protocols{protocolName1};
        protocolName2 = protocols{protocolName2};
        analysisChoice = analyses{analysisChoice};
        refElectrodes = refElectrodeGroups{refElectrodes};
        

        displaySignificanceFlag = 1;

      
        %refElectrodes = [16 17 18 48]; % back occipital
        groupType = 'rel';
        cutoffList = [1 20];
        pairedSubjectNameList = getPairedSubjectsBK1;            
        subjectNameLists{1} = pairedSubjectNameList(:,1);
        subjectNameLists{2} = pairedSubjectNameList(:,2);
        pairedDataFlag      = 1;
       

if ~exist('refElectrodes','var');         refElectrodes = [];           end
if ~exist('groupType','var');             groupType='rel';              end
if ~exist('connMethod','var');            connMethod = 'ppc';           end
if ~exist('badEyeCondition','var');       badEyeCondition='ep';         end
if ~exist('badTrialVersion','var');       badTrialVersion='v8';         end

if ~exist('freqRangeList','var')
    freqRangeList{1} = freqRange; % takes from GUI
    %freqRangeList{2} = [20 34]; % SG
    %freqRangeList{3} = [35 65]; % FG
end
numFreqRanges = length(freqRangeList);
if ~exist('axisRangeList','var')
    axisRangeList{1} = [0 100]; % freqLims
    axisRangeList{2} = [0 1]; % yLims for 
    axisRangeList{3} = [0 1]; % cLimsTopo
end
if ~exist('cutoffList','var')
    cutoffList = [2 30]; % elcs & trials
end

if ~exist('useMedianFlag','var');         useMedianFlag = 0;            end
if ~exist('pairedDataFlag','var');        pairedDataFlag = 0;           end
if ~exist('displayDataFlag','var');       displayDataFlag = 1;          end

saveFolderName = 'savedData';
capType = 'actiCap64_UOL';


[connData1,connData2, freqVals,connDataElectrodeGroup1,connDataElectrodeGroup2, electrodeGroupList,groupNameList,binnedCenters] = getConnDataAllSubjects(subjectNameLists,refElectrodes,groupType,connMethod,badEyeCondition,badTrialVersion,protocolName1,protocolName2, analysisChoice,cutoffList,pairedDataFlag,saveFolderName,capType);

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Display options %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
displaySettings.fontSizeLarge = 10;
displaySettings.tickLengthMedium = [0.025 0];

titleStr{1} = protocolName1;
titleStr{2} = protocolName2;
titleStr{3} = 'difference';
freqLims = axisRangeList{1};
yLims = axisRangeList{2};
cLimsTopo = axisRangeList{3};
numElectrodeGroups = length(electrodeGroupList);
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
%[~,~,~] = psdOnTopOfPPC(subjectNameLists,badEyeCondition,badTrialVersion,stRange, protocolPos, protocolPosRef,analysisChoice,badElectrodeRejectionFlag,cutoffList,pairedDataFlag, hAllPlots, displayDataFlag);

for m=1:numElectrodeGroups
    if strcmp(subjectsChoice, 'med')
    dataToPlot{1} = squeeze(connDataElectrodeGroup1{1}(:,m,:)); % med, protocol 1
    dataToPlot{2} = squeeze(connDataElectrodeGroup2{1}(:,m,:)); % med , protocol 2
    else
    dataToPlot{1} = squeeze(connDataElectrodeGroup1{2}(:,m,:)); % con, protocol 1
    dataToPlot{2} = squeeze(connDataElectrodeGroup2{2}(:,m,:)); % con , prot 2
    end
    
    if displayDataFlag
        % Connectivity as a function of frequency__________________
   displaySettings.colorNames(1,:) = [0 1 0]; % green for protocol 1
   displaySettings.colorNames(2,:) =  [1 0 0]; % red for protocol 2
   displaySettings.colorNames(3,:) =  [0 0 0]; % black for diff
        displayAndcompareData(hConnPPC(m), dataToPlot,freqVals,displaySettings,yLims,displaySignificanceFlag,useMedianFlag,~pairedDataFlag);
        sgtitle("PPC of " + subjectsChoice + " during " + protocolName1 + " and " + protocolName2 + " wrt elecs " + num2str(refElectrodes));
       
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
            xtemp{1} = connData1{1}; % prot 1
            xtemp{2} = connData2{1}; % prot 2
            xtemp{3} = xtemp{1} -xtemp{2}; % diff
    else 
            xtemp{1} = connData1{2}; % p1
            xtemp{2} = connData2{2}; % p2
            xtemp{3} = xtemp{1} -xtemp{2}; % diff
    end
for i=1:numFreqRanges

    dataToCompare = cell(1,2);
    % Topoplots_________

    for j=1:3
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
            dataToCompare{1} = squeeze(mean(connDataElectrodeGroup1{1}(:,:,freqPosList{i}),3));
            dataToCompare{2} = squeeze(mean(connDataElectrodeGroup2{1}(:,:,freqPosList{i}),3));
          else
            dataToCompare{1} = squeeze(mean(connDataElectrodeGroup1{2}(:,:,freqPosList{i}),3));
            dataToCompare{2} = squeeze(mean(connDataElectrodeGroup2{2}(:,:,freqPosList{i}),3));
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


    function [connData1,connData2, freqVals,connDataElectrodeGroup1,connDataElectrodeGroup2, electrodeGroupList,groupNameList,binnedCenters] = getConnDataAllSubjects(subjectNameLists,refElectrodes,groupType,connMethod,badEyeCondition,badTrialVersion,protocolName1, protocolName2,analysisChoice,cutoffList,pairedDataFlag,saveFolderName,capType)

[electrodeGroupList,groupNameList,binnedCenters] = getElectrodeGroupsConn(groupType,refElectrodes,capType);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize
badSubjectList = cell(1,2); % common list for both protocol 1 & 2

connDataTMP1 = cell(1,2);
connDataTMP2 = cell(1,2);
connDataElectrodeGroupTMP1 = cell(1,2);
connDataElectrodeGroupTMP2 = cell(1,2);
numRefElectrodes = length(refElectrodes);
numElectrodeGroups = length(electrodeGroupList);

for i=1:2 % med and con
    for j=1:length(subjectNameLists{i}) % no. of subjects in 'i'
        subjectName = subjectNameLists{i}{j};

        tmpData1 = load(fullfile(saveFolderName,subjectName,[protocolName1 '_' badEyeCondition '_' badTrialVersion '_' connMethod '.mat']));
        tmpData2 = load(fullfile(saveFolderName,subjectName,[protocolName2 '_' badEyeCondition '_' badTrialVersion '_' connMethod '.mat']));
        numGoodTrials1 = tmpData1.numGoodTrials;
        numGoodTrials2 = tmpData2.numGoodTrials;
        if numGoodTrials1 < cutoffList(2) || numGoodTrials2 < cutoffList(2)
            badSubjectList{i}(j) = 1;
        else
            if strcmp(analysisChoice,'bl')
                connDataTMPtmp1 = tmpData1.connPre(refElectrodes,:,:);
                connDataTMPtmp2 = tmpData2.connPre(refElectrodes,:,:);
                freqVals = tmpData1.freqPre;
            elseif strcmp(analysisChoice,'st')
                connDataTMPtmp1 = tmpData1.connPost(refElectrodes,:,:);
                connDataTMPtmp2 = tmpData2.connPost(refElectrodes,:,:);
                freqVals = tmpData1.freqPost;
            else
                connDataTMPtmp1 = (tmpData1.connPre(refElectrodes,:,:) + tmpData1.connPost(refElectrodes,:,:))/2;
                connDataTMPtmp2 = (tmpData2.connPre(refElectrodes,:,:) + tmpData2.connPost(refElectrodes,:,:))/2;
                freqVals = tmpData1.freqPost;
            end

            numGoodElectodes1 = trace(~isnan(squeeze(connDataTMPtmp1(:,refElectrodes,1))));
            numGoodElectodes2 = trace(~isnan(squeeze(connDataTMPtmp2(:,refElectrodes,1))));

            if numGoodElectodes1 >= cutoffList(1) && numGoodElectodes2 >= cutoffList(1)
                badSubjectList{i}(j) = 0;


                connDataTMP1{i}{j} = squeeze(mean(connDataTMPtmp1,1,'omitnan'));
                connDataTMP2{i}{j} = squeeze(mean(connDataTMPtmp2,1,'omitnan'));

                connDataElectrodeGroupTMPtmp1 = zeros(numRefElectrodes,numElectrodeGroups,length(freqVals));
                connDataElectrodeGroupTMPtmp2 = zeros(numRefElectrodes,numElectrodeGroups,length(freqVals));
                for e=1:numRefElectrodes
                    for b=1:numElectrodeGroups
                        if strcmp(groupType,'rel')
                            connDataElectrodeGroupTMPtmp1(e,b,:) = squeeze(mean(connDataTMPtmp1(e,electrodeGroupList{e,b},:),2,'omitnan'));
                            connDataElectrodeGroupTMPtmp2(e,b,:) = squeeze(mean(connDataTMPtmp2(e,electrodeGroupList{e,b},:),2,'omitnan'));
                        else
                            connDataElectrodeGroupTMPtmp1(e,b,:) = squeeze(mean(connDataTMPtmp1(e,electrodeGroupList{b},:),2,'omitnan'));
                             connDataElectrodeGroupTMPtmp2(e,b,:) = squeeze(mean(connDataTMPtmp2(e,electrodeGroupList{b},:),2,'omitnan'));
                        end
                    end
                end
                connDataElectrodeGroupTMP1{i}{j} = squeeze(mean(connDataElectrodeGroupTMPtmp1,1,'omitnan'));
                connDataElectrodeGroupTMP2{i}{j} = squeeze(mean(connDataElectrodeGroupTMPtmp2,1,'omitnan'));
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
connData1 = cell(1,3);
connData2 = cell(1,3);
connDataElectrodeGroup1 = cell(1,3);
connDataElectrodeGroup2 = cell(1,3);

for q=1:2
x1 = connDataTMP1{q};
x2 = connDataTMP2{q};
y1 = connDataElectrodeGroupTMP1{q};
y2 = connDataElectrodeGroupTMP2{q};
badSubjectPos = badSubjectList{q};
badSubjectPos = badSubjectPos(badSubjectPos<=length(connDataTMP1{q})); 
x1(badSubjectPos)=[];
x2(badSubjectPos)=[];
y1(badSubjectPos)=[]; 
y2(badSubjectPos)=[];  

numSubjects = length(x1);

    for p=1:numSubjects
        connData1{q}(p,:,:) = x1{p};
        connData2{q}(p,:,:) = x2{p};
        connDataElectrodeGroup1{q}(p,:,:) = y1{p};
        connDataElectrodeGroup2{q}(p,:,:) = y2{p};
    end
end
connData1{3} = connData1{1} - connData1{2}; % diff: med-con
connData2{3} = connData2{1} - connData2{2};
connDataElectrodeGroup1{3} = connDataElectrodeGroup1{1} - connDataElectrodeGroup1{2};
connDataElectrodeGroup2{3} = connDataElectrodeGroup2{1} - connDataElectrodeGroup2{2};
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

%%%%%%%%%%%%% functioms %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% gives and plots customized PSD data to be plotted on top of PPC curves -- called by displayConnDataAllSubjects
% the logPSD values are scaled down and shifted up to fit in the PPC plots
function [psdDataToReturn, powerDataToReturn, freqVals] = psdOnTopOfPPC(subjectNameLists,badEyeCondition,badTrialVersion,stRange, protocolPos, protocolPosRef,analysisChoice,badElectrodeRejectionFlag,cutoffList,pairedDataFlag, hAllPlots, displayDataFlag)
if ~exist('savePSDFolderName', 'var');    savePSDFolderName ='savedPSDDataTill200' ;    end
if ~exist('protocolName','var');          protocolName='G1';            end
if ~exist('analysisChoice','var');        analysisChoice='st';          end
if ~exist('refChoice','var');             refChoice='none';             end
if ~exist('groupType','var');             groupType='rel';              end % electrode groups based on relative distance from Seed (ref) electrodes
if ~exist('badEyeCondition','var');       badEyeCondition='ep';         end
if ~exist('badTrialVersion','var');       badTrialVersion='v8';         end
if ~exist('badElectrodeRejectionFlag','var'); badElectrodeRejectionFlag=1;  end

if ~exist('stRange','var');               stRange = [0.25 1.25];        end

if ~exist('freqRangeList','var')
    freqRangeList{1} = [8 13]; % alpha
    freqRangeList{2} = [20 34]; % SG
    freqRangeList{3} = [35 65]; % FG
end
if ~exist('axisRangeList','var')
    axisRangeList{1} = [0 100];
    axisRangeList{2} = [-2.5 2.5];
    axisRangeList{3} = [-1.5 1.5];
end
if ~exist('cutoffList','var')
    cutoffList = [2 30];
end
cutoffNumElectrodes = cutoffList(1);
cutoffNumTrials = cutoffList(2);

if ~exist('useMedianFlag','var');         useMedianFlag = 0;            end
if ~exist('hAllPlots','var');             hAllPlots = [];               end
if ~exist('pairedDataFlag','var');        pairedDataFlag = 1;           end
if ~exist('displayDataFlag','var');       displayDataFlag = 1;          end

numFreqRanges = length(freqRangeList);
freqRangeColors = copper(numFreqRanges);
refElectrodes = [16 17 18 48]; % back occipital group

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Display options %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
displaySettings.fontSizeLarge = 10;
displaySettings.tickLengthMedium = [0.025 0];
displaySettings.colorNames(1,:) = [0.8 0 0.8];      % Purple 
displaySettings.colorNames(2,:) = [0.25 0.41 0.88]; % Cyan
titleStr{1} = 'Meditators';
titleStr{2} = 'Controls';

freqLims = axisRangeList{1};
%yLimsPSD = axisRangeList{2};
yLimsPSD = [-1 1.1];
cLimsTopo = axisRangeList{3};

%%%%%%%%%%%%%%%%%%%%%%%% Get electrode groups %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gridType = 'EEG';
capType = 'actiCap64_UOL';
% savePSDFolderName = 'savedDataTill200';

[electrodeGroupList,groupNameList] = getElectrodeGroupsConn(gridType,refElectrodes, capType);
numGroups = length(electrodeGroupList);

   displaySettings.colorNames(1,:) = [0 1 0]; % green 
   displaySettings.colorNames(2,:) =  [1 0 0]; % red 
%%%%%%%%%%%%%%%%%%%%%%%%%%%% plot handles %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if displayDataFlag
    if isempty(hAllPlots)
        hPSD  = getPlotHandles(1,numGroups,[0.05 0.55 0.6 0.3],0.02,0.02,1);
        %hPower = getPlotHandles(numFreqRanges,numGroups,[0.05 0.05 0.6 0.45],0.02,0.02,0);
        hTopo0 = getPlotHandles(1,2,[0.675 0.7 0.3 0.15],0.02,0.02,1);
        %hTopo1 = getPlotHandles(1,3,[0.675 0.55 0.3 0.13],0.02,0.02,1);
        %hTopo2 = getPlotHandles(numFreqRanges,3,[0.675 0.05 0.3 0.45],0.02,0.02,1);
    else
        hPSD = hAllPlots.hConn1;
        % hPower = hAllPlots.hPower;
        hTopo0 = hAllPlots.hTopoRef;
        % hTopo1 = hAllPlots.hTopo1;
        % hTopo2 = hAllPlots.hTopo2;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% Show Electrodes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if displayDataFlag
    % Display electrode positions
    if strcmp(groupType,'rel')  
        montageChanlocs = showElectrodeGroups(hTopo0,capType,electrodeGroupList(1,:),groupNameList); % Just show for the first electrode
    else
        montageChanlocs = showElectrodeGroups(hTopo0,capType,electrodeGroupList,groupNameList);
    end
end
end

        pairedSubjectNameList = getPairedSubjectsBK1;            
        subjectNameLists{1} = pairedSubjectNameList(:,1);
        subjectNameLists{2} = pairedSubjectNameList(:,2);
      

goodSubjectNameLists = getGoodSubjectNameList(subjectNameLists,badEyeCondition,badTrialVersion,stRange,protocolPos,protocolPosRef,analysisChoice,badElectrodeRejectionFlag,cutoffNumTrials,pairedDataFlag,savePSDFolderName);
[powerData,powerDataRef,freqVals] = getPowerDataAllSubjects(goodSubjectNameLists,badEyeCondition,badTrialVersion,stRange,protocolPos,protocolPosRef,analysisChoice,badElectrodeRejectionFlag,cutoffNumTrials,savePSDFolderName);

psdDataToReturn = cell(1,numGroups);
powerDataToReturn = cell(numGroups,numFreqRanges);
goodSubjectNameListsToReturn = cell(numGroups,2);

for i=1:numGroups

    %%%%%%%%%%%%% Find bad subjects based on electrode cutoff %%%%%%%%%%%%%
    badSubjectPosList = cell(1,2);
    for j=1:2
        pData = powerData{j}(electrodeGroupList{i},:,:);
        numGoodElecs = length(electrodeGroupList{i}) - sum(isnan(squeeze(pData(:,1,:))),1);
        badSubjectPosList{j} = find(numGoodElecs<cutoffNumElectrodes);

        if ~isempty(protocolPosRef)
            pDataRef = powerDataRef{j}(electrodeGroupList{i},:,:);
            numGoodElecsRef = length(electrodeGroupList{i}) - sum(isnan(squeeze(pDataRef(:,1,:))),1);
            badSubjectPosRef = find(numGoodElecsRef<=cutoffNumElectrodes);
            badSubjectPosList{j} = unique(cat(2,badSubjectPosList{j},badSubjectPosRef));
        end
    end

    % if paired dataset, take common bad subjects
    if pairedDataFlag
        badSubjectPosList{1} = union(badSubjectPosList{1},badSubjectPosList{2});
        badSubjectPosList{2} = badSubjectPosList{1};
    end

    % Get data
    meanPSDData = cell(1,2);
    meanPSDDataRef = cell(1,2);
    logPSDData = cell(1,2);

    for j=1:2
        pData = powerData{j}(electrodeGroupList{i},:,:);
        if ~isempty(protocolPosRef)
            pDataRef = powerDataRef{j}(electrodeGroupList{i},:,:);
        end
        badSubjectPos = badSubjectPosList{j};

        tmp = goodSubjectNameLists{j};
        tmp(badSubjectPos) = [];
        goodSubjectNameListsToReturn{i,j} = tmp;

        if ~isempty(badSubjectPos)       
            disp([groupNameList{i} ', ' titleStr{j} ', '  'Not enough good electrodes for ' num2str(length(badSubjectPos)) ' subjects.']);
            pData(:,:,badSubjectPos)=[];
            if ~isempty(protocolPosRef)
                pDataRef(:,:,badSubjectPos)=[];
            end
        end
        meanPSDData{j} = squeeze(mean(pData,1,'omitnan'))';

        if isempty(protocolPosRef)
            logPSDData{j} = log10(meanPSDData{j});
            logPSDData{j} = logPSDData{j}./5; % scaled down by 1/5 to fit in PPC plots
            logPSDData{j} = logPSDData{j}+1;
        else
            meanPSDDataRef{j} = squeeze(mean(pDataRef,1,'omitnan'))';
            logPSDData{j} = 10*(log10(meanPSDData{j}) - log10(meanPSDDataRef{j}));
            logPSDData{j} = logPSDData{j}./2; % scaled down ny 1/5 to fit in PPC plots
            logPSDData{j} = logPSDData{j}+0.5;
        end

        if displayDataFlag
            text(10,yLimsPSD(2)-0.1*j,[titleStr{j} '(' num2str(size(meanPSDData{j},1)) ')'],'color',displaySettings.colorNames(j,:),'parent',hPSD(i));
        end
    end

    psdDataToReturn{i} = logPSDData;
    if displayDataFlag
        displayAndcompareData(hPSD(i),logPSDData,freqVals,displaySettings,yLimsPSD,0,useMedianFlag,~pairedDataFlag);
        title(groupNameList{i});
        xlim(hPSD(i),freqLims);
    end
end

end
%%%%%%%%%%%%%%% functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [powerData,powerDataRef,freqVals] = getPowerDataAllSubjects(subjectNameLists,badEyeCondition,badTrialVersion,stRange,protocolPos,protocolPosRef,analysisChoice,badElectrodeRejectionFlag,cutoffNumTrials, savePSDFolderName)

powerData = cell(1,2);
powerDataRef = cell(1,2);

for i=1:2
    powerDataTMP=[];
    powerDataRefTMP=[];

    for j=1:length(subjectNameLists{i})
        subjectName = subjectNameLists{i}{j};
        % using psd data up to 200Hz, shared by Srishty
        
        tmpData = load(fullfile(savePSDFolderName,[subjectName '_' badEyeCondition '_' badTrialVersion '_' num2str(1000*stRange(1)) '_' num2str(1000*stRange(2))]));
        freqVals = tmpData.freqVals;

        tmpPower = getPowerData(tmpData,protocolPos,analysisChoice,badElectrodeRejectionFlag,cutoffNumTrials);
        powerDataTMP = cat(3,powerDataTMP,tmpPower);

        if ~isempty(protocolPosRef)
            tmpPowerRef = getPowerData(tmpData,protocolPosRef,'bl',badElectrodeRejectionFlag,cutoffNumTrials);
            powerDataRefTMP = cat(3,powerDataRefTMP,tmpPowerRef);
        end
    end
    powerData{i} = powerDataTMP;
    powerDataRef{i} = powerDataRefTMP;
end
end
function tmpPower = getPowerData(tmpData,protocolPos,analysisChoice,badElectrodeRejectionFlag,cutoffNumTrials)

numTrials = tmpData.numTrials(protocolPos);
badElectrodes = getBadElectrodes(tmpData.badElectrodes,badElectrodeRejectionFlag,protocolPos);

if numTrials < cutoffNumTrials
    tmpPower = [];
else
    if strcmpi(analysisChoice,'st')
        tmpPower = tmpData.psdValsST{protocolPos};
    elseif strcmpi(analysisChoice,'bl')
        tmpPower = tmpData.psdValsBL{protocolPos};
    else
        tmpPower = (tmpData.psdValsST{protocolPos}+tmpData.psdValsBL{protocolPos})/2; % average
    end
    tmpPower(badElectrodes,:) = NaN;
end
end

function badElectrodes = getBadElectrodes(badElectrodeList,badElectrodeRejectionFlag,protocolPos)

if badElectrodeRejectionFlag==1 % Bad electrodes for the protocol
    badElectrodes = badElectrodeList{protocolPos};
elseif badElectrodeRejectionFlag==2 % common bad electrodes for all protocols
    badElectrodes=[];
    for i=1:length(badElectrodeList)
        badElectrodes=cat(1,badElectrodes,badElectrodeList{i});
    end
    badElectrodes = unique(badElectrodes);
elseif badElectrodeRejectionFlag==3 % Bad electrodes of G1
    badElectrodes = badElectrodeList{3};
end
end
function goodSubjectNameLists = getGoodSubjectNameList(subjectNameLists,badEyeCondition,badTrialVersion,stRange,protocolPos,protocolPosRef,analysisChoice,badElectrodeRejectionFlag,cutoffNumTrials,pairedDataFlag,savePSDFolderName)

% For unpaired case, subjects can be rejected if either data in analysis or
% ref period is bad. For paired case, a pair is rejected even if one of the two
% subjects in the pair is bad. Based on the condition, we get a new good
% subject list.

badSubjectIndex = cell(1,2);
badSubjectIndexRef = cell(1,2);

for i=1:2
    numSubjects = length(subjectNameLists{i});
    badSubjectIndexTMP = zeros(1,numSubjects);
    badSubjectIndexRefTMP = zeros(1,numSubjects);

    for j=1:numSubjects
        subjectName = subjectNameLists{i}{j};
        tmpData = load(fullfile(savePSDFolderName,[subjectName '_' badEyeCondition '_' badTrialVersion '_' num2str(1000*stRange(1)) '_' num2str(1000*stRange(2))]));

        if isempty(getPowerData(tmpData,protocolPos,analysisChoice,badElectrodeRejectionFlag,cutoffNumTrials))
            disp(['Not enough trials for subject: ' subjectName]);
            badSubjectIndexTMP(j)=1;
        end

        if ~isempty(protocolPosRef)
            if isempty(getPowerData(tmpData,protocolPosRef,'bl',badElectrodeRejectionFlag,cutoffNumTrials))
                disp(['Not enough trials in ref period for subject: ' subjectName]);
                badSubjectIndexRefTMP(j)=1;
            end
        end
    end
    badSubjectIndex{i} = badSubjectIndexTMP;
    badSubjectIndexRef{i} = badSubjectIndexRefTMP;
end

%%%%%%%%%%%%%%%%%%%%%%%%% Now find good subjects %%%%%%%%%%%%%%%%%%%%%%%%%%
goodSubjectNameLists = cell(1,2);

if ~pairedDataFlag
    for i=1:2
        subjectNameListTMP = subjectNameLists{i};
        if isempty(protocolPosRef)
            badPos = find(badSubjectIndex{i});
        else
            badPos = union(find(badSubjectIndex{i}),find(badSubjectIndexRef{i}));
        end
        subjectNameListTMP(badPos)=[];
        goodSubjectNameLists{i} = subjectNameListTMP;
    end
else
    
    if isempty(protocolPosRef)
        badPos = find(sum(cell2mat(badSubjectIndex')));
    else
        badPos = union(find(sum(cell2mat(badSubjectIndex'))),find(sum(cell2mat(badSubjectIndexRef'))));
    end

    for i=1:2
        subjectNameListTMP = subjectNameLists{i};
        subjectNameListTMP(badPos)=[];
        goodSubjectNameLists{i} = subjectNameListTMP;
    end
end
end