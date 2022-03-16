%% aim
%{

%}
%% define qualifying areas
clear all;
cellUseAreas = {...
	'Primary visual area',...
	'posteromedial visual area',...
	};


%% select all neurons in LP and drifting grating stimuli
if ~exist('sAggStim','var') || isempty(sAggStim)
	[sAggStim,sAggNeuron]=loadDataNpx('','driftinggrating');
end
sAggNeuron(strcmpi({sAggNeuron.SubjectType},'DBA')) = [];

%% pre-allocate matrices
intAreas = numel(cellUseAreas);
matR_PP_All = nan(intAreas,intAreas,numel(sAggStim),2);
matR_OP_All = nan(intAreas,intAreas,numel(sAggStim),2);
matR_OO_All = nan(intAreas,intAreas,numel(sAggStim),2);
matR_PO_All = nan(intAreas,intAreas,numel(sAggStim),2);
mat_xR_All = nan(intAreas,intAreas,numel(sAggStim),2);
matDecPerf = [];

%% go through recordings
tic
for intRec=1:numel(sAggStim)
	% get matching recording data
	strRec = sAggStim(intRec).Exp;
	sThisRec = sAggStim(strcmpi(strRec,{sAggStim(:).Exp}));
	
	%remove stimulus sets that are not 24 stim types
	sThisRec.cellBlock(cellfun(@(x) x.intTrialNum/x.intNumRepeats,sThisRec.cellBlock) ~= 24) = [];
	% concatenate stimulus structures
	structStim = catstim(sThisRec.cellBlock);
	vecStimOnTime = structStim.vecStimOnTime;
	vecStimOffTime = structStim.vecStimOffTime;
	vecOrientation = cell2vec({structStim.sStimObject(structStim.vecTrialStimTypes).Orientation})';
	[vecOriIdx,vecUniqueOris,vecRepNum,cellSelect,vecTrialRepetition]  = label2idx(vecOrientation);
	indRem=vecTrialRepetition>min(vecRepNum);
	vecOrientation(indRem) = [];
	vecStimOnTime(indRem) = [];
	vecStimOffTime(indRem) = [];
	[vecOriIdx,vecUniqueOris,vecRepNum,cellSelect,vecTrialRepetition]  = label2idx(vecOrientation);
	intTrialNum = numel(vecStimOnTime);
	intOriNum = numel(unique(vecOrientation));
	intRepNum = intTrialNum/intOriNum;
	
	%remove neurons from other recs
	indQualifyingNeurons = contains({sAggNeuron.Exp},strRec);
	
	%remove neurons in incorrect areas
	indConsiderNeurons = contains({sAggNeuron.Area},cellUseAreas,'IgnoreCase',true);
	
	%remove bad clusters
	indGoodNeurons = (cell2vec({sAggNeuron.KilosortGood}) == 1) | (cell2vec({sAggNeuron.Contamination}) < 0.1);
	
	%subselect from total
	indUseNeurons = indQualifyingNeurons(:) & indConsiderNeurons(:) & indGoodNeurons(:);
	sUseNeuron = sAggNeuron(indUseNeurons);
	
	%% select area 1
	for intArea=1:intAreas
		strArea = cellUseAreas{intArea};
		indArea1Neurons = contains({sUseNeuron.Area},strArea,'IgnoreCase',true);
		if sum(indArea1Neurons) == 0, continue;end
		%% get orientation responses & single-trial population noise
		sArea1Neurons = sUseNeuron(indArea1Neurons);
		
		%% get spike times
		cellSpikeTimes = {sArea1Neurons.SpikeTimes};
		intNumN = numel(cellSpikeTimes);
		dblStimDur = roundi(median(vecStimOffTime - vecStimOnTime),1,'ceil');
		dblPreTime = 0.5;
		dblPostTime = 0.5;
		dblMaxDur = dblStimDur+dblPreTime+dblPostTime;
		dblBinWidth = 0.1;
		vecBinEdges = 0:dblBinWidth:dblMaxDur;
		vecStimTime = vecBinEdges(2:end)-dblBinWidth/2 - dblPreTime;
		intBinNum = numel(vecBinEdges)-1;
		matBNSR = nan(intBinNum,intNumN,intOriNum,intRepNum);
		matBNT = nan(intBinNum,intNumN,intTrialNum);
		
		vecRepCounter = zeros(1,intOriNum);
		for intN=1:intNumN
			[vecTrialPerSpike,vecTimePerSpike] = getSpikesInTrial(cellSpikeTimes{intN},vecStimOnTime-dblPreTime,dblMaxDur);
			for intTrial=1:intTrialNum
				if intN==1
					intTrialOriIdx = vecOriIdx(intTrial);
					vecRepCounter(intTrialOriIdx) = vecRepCounter(intTrialOriIdx) + 1;
					intRep = vecRepCounter(intTrialOriIdx);
				end
				vecSpikeT = vecTimePerSpike(vecTrialPerSpike==intTrial);
				vecSpikeC = histcounts(vecSpikeT,vecBinEdges);
				matBNSR(:,intN,intTrialOriIdx,intRep) = vecSpikeC;
				matBNT(:,intN,intTrial) = vecSpikeC;
			end
		end
		
		%% time progression
		matDecPerf(:,end+1)=nan;
		dblLambda = 1;
		intTypeCV = 2;
		vecOri180 = mod(vecOrientation,180)*2;
		vecRepNum180 = vecRepNum(1:12)*2;
		intOriNum180 = intOriNum/2;
		matDecConfusion = nan(intOriNum180,intOriNum180,intBinNum);
		vecSpikesPerBin = sum(sum(matBNT,2),3)';
		[varDataOut,vecUnique,vecPriorDistribution,cellSelect,vecRepetition] = label2idx(vecOri180);
		for intBinIdx=1:intBinNum
			intBinIdx
			[dblPerformanceCV,vecDecodedIndexCV,matPosteriorProbability,dblMeanErrorDegs,matConfusion,matWeights] = ...
				doCrossValidatedDecodingLR(squeeze(matBNT(intBinIdx,:,:)),vecOri180,intTypeCV,vecPriorDistribution,dblLambda);
			vecDecErr(intBinIdx) = dblMeanErrorDegs;
			matDecPerf(intBinIdx,end) = dblPerformanceCV;
			matDecConfusion(:,:,intBinIdx) = matConfusion;
		end
		%%
		dblChance = 1/numel(vecPriorDistribution);
		figure;maxfig;
		subplot(2,3,1)
		plot(vecStimTime,matDecPerf(:,end)');
		hold on
		plot([vecStimTime(1) vecStimTime(end)],[dblChance dblChance],'--','color',[0.5 0.5 0.5]);
		hold off
		title(sprintf('Dec perf; %s',strArea))
		xlabel('Time after onset (s)');
		ylabel('Fraction correct decoded');
		fixfig;
		
		subplot(2,3,2)
		plot(vecStimTime,vecSpikesPerBin)
		title('Spikes per bin')
		xlabel('Time after onset (s)');
		ylabel('Spikes per bin');
		fixfig;
		
		subplot(2,3,3)
		plot(vecStimTime,(matDecPerf(:,end)-dblChance)'./vecSpikesPerBin)
		title('Dec perf / spike')
		xlabel('Time after onset (s)');
		ylabel('Performance/spike');
		fixfig;
	end
end
toc
