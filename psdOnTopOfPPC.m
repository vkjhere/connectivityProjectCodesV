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