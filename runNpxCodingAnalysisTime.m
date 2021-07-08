%% aim
%{
"distinct subnetworks in mouse visual system for visual and non-visual signals"

At single trial basis, split which info information limiting noise
(projected unto f') and non-limiting noise (orthogonal directions). Then
look at correlation of those noise types between areas. E.g., lgn diff
corrs correlate with v1, but SC non-diff corrs correlate with v1


Introduce with: spike count correlations are generally low, but important.
Show how to decompose noise correlations on single trials and that
orth>para. Show that inter-areal correlations of shuffles are 0, then show
real data: all red. Trial-to -trial variability in MD coding is highly
correlated despite pairwise noise correlations being low

%}
%% define qualifying areas
cellUseAreas = {...
	'Subiculum',...
	'Posterior complex of the thalamus',...
	'Anterior pretectal nucleus',...
	'Superior colliculus',...
	'Lateral posterior nucleus of the thalamus',...
	'Nucleus of the optic tract',...
	'Dorsal part of the lateral geniculate complex',...
	'Primary visual area',...
	'posteromedial visual area',...
	'Anteromedial visual area',...
	'Retrosplenial area',...
	'Anterior area',...
	};


%% select all neurons in LP and drifting grating stimuli
if ~exist('sAggStim','var') || isempty(sAggStim)
	[sAggStim,sAggNeuron]=loadDataNpx('','driftinggrating');
end

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
for intRec=10:numel(sAggStim)
	% get matching recording data
	strRec = sAggStim(intRec).Rec;
	sThisRec = sAggStim(strcmpi(strRec,{sAggStim(:).Rec}));
	
	%remove stimulus sets that are not 24 stim types
	sThisRec.cellStim(cellfun(@(x) x.structEP.intStimTypes,sThisRec.cellStim) ~= 24) = [];
	% concatenate stimulus structures
	structStim = catstim(sThisRec.cellStim);
	vecStimOnTime = structStim.vecStimOnTime;
	vecStimOffTime = structStim.vecStimOffTime;
	vecOrientation = structStim.Orientation;
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
	indQualifyingNeurons = contains({sAggNeuron.Rec},strRec);
	
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
		dblStimDur = roundi(min(vecStimOffTime - vecStimOnTime),1,'ceil');
		dblPreTime = 0.5;
		dblPostTime = 0.5;
		dblMaxDur = dblStimDur+dblPreTime+dblPostTime;
		vecBinEdges = 0:0.1:dblMaxDur;
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
		for intBinIdx=1:intBinNum
			intBinIdx
			[dblPerformance,vecDecodedIndexCV,matPosteriorProbability,matWeights,dblMeanErrorDegs,matConfusion] = ...
				doCrossValidatedDecodingLR(squeeze(matBNT(intBinIdx,:,:)),vecOri180,intTypeCV,dblLambda);
			vecDecErr(intBinIdx) = dblMeanErrorDegs;
			matDecPerf(intBinIdx,end) = dblPerformance;
			matDecConfusion(:,:,intBinIdx) = matConfusion;
		end
		%%
		figure
		plot(zscore(matDecPerf(:,end)./vecSpikesPerBin))
		hold on
		plot(zscore(vecSpikesPerBin))
		hold off
		title(sprintf('%s',strArea))
		%pause
	end
end
toc
