%% aim
%{
"distinct subnetworks in mouse visual system for visual and non-visual signals"

Predict noise in area 2 from noise in area 1, calculate angle between
prediction and real, and compare with shuffled. If real diff is lower, then
we found inter-areal structure of high-d noise. Repeat after removing
independent noise; does difference hold?

%}
%% aim 2
%{
"non-visual activity in visual system is globally invariant to stimulus orientation"

1) Predict noise in area y from noise in area x for stimulus z using stimulus
z; compare with using stimulus z+1,+2...
2) compare this baseline to shuffled where generalizability is destroyed (how?)
3) compare also to predictability where z+1 is rotated to z, using
difference in f-primes between z and z+1

Possibly need to reduce to predictive subspace + linear regression as in Semedo et al., fig S6

%}
%% to do
%{
it now predicts noise & performs decoding within and across areas, but:
i) we need more repetitions for more areas
ii) analysis needs to be extended to predict across stimuli
iii) analysis needs to be extended to predict across stimuli after rotation

=> we can compare predictability of noise between areas; e.g., is subiculum
driving non-visual noise?

%}
%% define qualifying areas
vecN = [];
vecT = [];
strFigDir = 'F:\Data\Results\NpxDims\';

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

%% go through recordings
tic
for intRec=1:numel(sAggStim)
	% get matching recording data
	strRec = sAggStim(intRec).Rec;
	sThisRec = sAggStim(strcmpi(strRec,{sAggStim(:).Rec}));
	
	%get short ID
	cellSplit = strsplit(strRec,'_');
	strRecID = strjoin(cellSplit(1:2),'_');
	
	%remove stimulus sets that are not 24 stim types
	sThisRec.cellStim(cellfun(@(x) x.structEP.intStimTypes,sThisRec.cellStim) ~= 24) = [];
	% concatenate stimulus structures
	structStim = catstim(sThisRec.cellStim);
	vecStimOnTime = structStim.vecStimOnTime;
	vecStimOffTime = structStim.vecStimOffTime;
	vecOrientation = structStim.Orientation;
	
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
	for intArea1=1:intAreas
		strArea1 = cellUseAreas{intArea1};
		indArea1Neurons = contains({sUseNeuron.Area},strArea1,'IgnoreCase',true);
		if sum(indArea1Neurons) < 40 || numel(vecStimOnTime) < 500, continue;end
		vecN(end+1) = sum(indArea1Neurons);
		vecT(end+1) = numel(vecStimOnTime);
		
		%% get orientation responses & single-trial population noise
		sArea1Neurons = sUseNeuron(indArea1Neurons);
		
		%% get spike times
		cellSpikeTimes1 = {sArea1Neurons.SpikeTimes};
		matSpikeCountsArea1 = getSpikeCounts(cellSpikeTimes1,vecStimOnTime,vecStimOffTime);
		
		%% split data
		%prep data
		matData = matSpikeCountsArea1;
		vecTrialStimType = vecOrientation;
		%remove nulls
		matNullNeurons = doDimRemoveNulls(matData,vecTrialStimType);
		indRem = any(matNullNeurons,2);
		matData(indRem,:) = [];
		
		%prep params
		vecUseStimTypes = unique(vecOrientation);
		sParams = struct;
		sParams.intSizeX = 25;
		sParams.intSizeY = 10;
		sParams.vecUseStimTypes = vecUseStimTypes;
		sParams.intWithinArea = 1;
		sParams.vecCellArea = ones(1,size(matData,1));
		%split
		[cellMatX,	cellNeuronsX,cellMatY,cellNeuronsY,cellTrials] = doDimDataSplits(matData,vecTrialStimType,sParams);
		%dblLambda = 1;
		%matPredictionsR2 = doDimPredFull(cellMatX,cellMatY,dblLambda);
		%mean(matPredictionsR2(:))
		%matPredDimDepR2 = doDimPredRR(cellMatX,cellMatY)
		
		%% find optimal regularisation parameter for ridge regression
		%folds
		intFoldK = 10;
		
		%line search
		vecLambda = 2.^[0:0.2:20];%[0:(ceil(log2(dblLambda))+2)];
		vecLambdaPerf = nan(size(vecLambda));
		for intL=1:numel(vecLambda)
			dblL = vecLambda(intL);
			[vecR2_Ridge_CV, sRidge] = doRidge_CV(cellMatX{1,1}, cellMatY{1,1}, intFoldK, dblL);
			vecLambdaPerf(intL) = mean(vecR2_Ridge_CV);
		end
		intW = 10;
		vecF = [nan([1 intW]) vecLambdaPerf nan([1 intW])];
		vecLPF = filter((1/intW)*ones(1,intW),1,vecF);
		vecLPF = vecLPF((intW+1):(end-intW));
		[dummy,intP] = max(vecLPF);
		dblLambda = vecLambda(intP);
		
		%plot
		%scatter(vecLambda,vecLambdaPerf);
		%set(gca,'xscale','log')
		%% decode
		matDecData = cell2mat(cellMatX(1,:)');
		vecDecStimType = cell2vec(cellfun(@(x,y) y*ones(1,x),vec2cell(cellfun(@(x) size(x,1),cellMatX(1,:))),vec2cell(1:size(cellMatX,2)),'UniformOutput',false));
		
		
		[dblPerformance1,vecDecodedIndexCV,matPosteriorProbability,matWeights,dblMeanErrorDegs,matConfusion1] = ...
			doCrossValidatedDecodingLR(matData,vecTrialStimType,2,dblLambda);
		
		[dblPerformance2,vecDecodedIndexCV,matPosteriorProbability,matWeights,dblMeanErrorDegs,matConfusion2] = ...
			doCrossValidatedDecodingLR(matDecData,vecDecStimType,2,dblLambda);
		
		close all;
		figure
		subplot(2,3,1)
		imagesc(matConfusion1);colorbar;
		title(sprintf('%s %s, full set, n=%d, diag=%.3f',strRecID,strArea1,size(matData,1),dblPerformance1),'interpreter','none')
		xlabel('Real stim');
		ylabel('Decoded stim #');
		fixfig;
		grid off
		
		subplot(2,3,2)
		imagesc(matConfusion2);colorbar;
		title(sprintf('split set, n=%d, diag=%.3f',size(matDecData,2),dblPerformance2),'interpreter','none')
		xlabel('Real stim');
		ylabel('Decoded stim #');
		fixfig;
		grid off
		
		%% run self-prediction
		intFoldK = 10;
		%dblLambda = 1/3000;
		vecRank = 1:sParams.intSizeY;
		tic
		dblLambda1 =1/100;
		dblLambda2 =100;
		warning('off','MATLAB:sqrtm:SingularMatrix');
		matPredictionsR2_1 = doDimPredRRR_CV(cellMatX,cellMatY,vecRank,intFoldK,dblLambda1);
		matPredictionsR2_2 = doDimPredRRR_CV(cellMatX,cellMatY,vecRank,intFoldK,dblLambda2);
		warning('on','MATLAB:sqrtm:SingularMatrix');
		toc
		intReps = size(cellMatX{1},1);
		
		subplot(2,3,4)
		plot(squeeze(mean(matPredictionsR2_1,1))')
		title(sprintf('%s, CV RRRR, %s=%.0f, reps=%d',strArea1,getGreek('lambda'),dblLambda1,intReps),'interpreter','none')
		xlabel('Dimensionality');
		ylabel('Noise prediction, R^2');
		fixfig;
		
		subplot(2,3,5)
		plot(squeeze(mean(matPredictionsR2_2,1))')
		title(sprintf('CV RRRR, %s=%.0f, reps=%d',getGreek('lambda'),dblLambda2,intReps),'interpreter','none')
		xlabel('Dimensionality');
		ylabel('Noise prediction, R^2');
		fixfig;
		
		strFigFile = sprintf('SelfPred_%s_%s_%s',strRecID,strArea1,getDate);
		maxfig;
		drawnow;
		export_fig([strFigDir strFigFile '.tif']);
		export_fig([strFigDir strFigFile '.pdf']);
		%vecR2 = getRegInSpace(matX,matSubspace,matY,dblLambda)
		
		%% select area 2
		for intArea2=1:intAreas
			strArea2 = cellUseAreas{intArea2};
			indArea2Neurons = contains({sUseNeuron.Area},strArea2,'IgnoreCase',true);
			if sum(indArea2Neurons) == 0 || intArea1 == intArea2, continue;end
			
			%% get orientation responses & single-trial population noise
			sArea2Neurons = sUseNeuron(indArea2Neurons);
			
			%% get spike times
			cellSpikeTimes2 = {sArea2Neurons.SpikeTimes};
			matSpikeCountsArea2 = getSpikeCounts(cellSpikeTimes2,vecStimOnTime,vecStimOffTime);
			
			%remove nulls
			matData2 = matSpikeCountsArea2;
			matNullNeurons = doDimRemoveNulls(matData2,vecTrialStimType);
			indRem = any(matNullNeurons,2);
			matData2(indRem,:) = [];
			if size(matData2,1) < (sParams.intSizeY + 2),continue;end
			
			%split data
			sParams.vecCellArea = cat(1,ones(size(matData,1),1),2*ones(size(matData2,1),1));
			sParams.intWithinArea = 0;
			matDataCombined = cat(1,matData,matData2);
			[cellMatX,	cellNeuronsX,cellMatY,cellNeuronsY,cellTrials] = doDimDataSplits(matDataCombined,vecTrialStimType,sParams);
			
			%% run across prediction
			intFoldK = 10;
			%dblLambda = 1/3000;
			vecRank = 1:sParams.intSizeY;
			tic
			dblLambda1 =1/100;
			dblLambda2 =100;
			warning('off','MATLAB:sqrtm:SingularMatrix');
			matPredictionsR2_1 = doDimPredRRR_CV(cellMatX,cellMatY,vecRank,intFoldK,dblLambda1);
			matPredictionsR2_2 = doDimPredRRR_CV(cellMatX,cellMatY,vecRank,intFoldK,dblLambda2);
			warning('on','MATLAB:sqrtm:SingularMatrix');
			toc
			intReps = size(cellMatX{1},1);
			
			figure;
			subplot(2,3,1)
			plot(squeeze(mean(matPredictionsR2_1,1))')
			title(sprintf('%s,%s => %s, %s=%.0f',strRecID,strArea1,strArea2,getGreek('lambda'),dblLambda1),'interpreter','none')
			xlabel('Dimensionality');
			ylabel('Noise prediction, R^2');
			fixfig;
			
			subplot(2,3,2)
			plot(squeeze(mean(matPredictionsR2_2,1))')
			title(sprintf('CV RRRR, %s=%.0f, reps=%d',getGreek('lambda'),dblLambda2,intReps),'interpreter','none')
			xlabel('Dimensionality');
			ylabel('Noise prediction, R^2');
			fixfig;
			
			strFigFile = sprintf('xPred_%s_%s_%s_%s',strRecID,strArea1,strArea2,getDate);
			maxfig;
			drawnow;
			export_fig([strFigDir strFigFile '.tif']);
			export_fig([strFigDir strFigFile '.pdf']);
		end
	end
end
