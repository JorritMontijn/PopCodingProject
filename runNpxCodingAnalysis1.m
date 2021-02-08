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
cellRepStr = {...
	'RunDriftingGratings','-DG';...
	'RunNaturalMovie','-NM';...
	'lateral geniculate','LGN';...
	'Primary visual','V1';...
	'Lateral posterior nucleus','LP';...
	'Anterior pretectal nucleus','APN';...
	'Posterior complex of the thalamus','PCT';...
	'Nucleus of the optic tract','NOT';...
	'Superior colliculus','SC';...
	'Anteromedial visual','AM';...
	'posteromedial visual','PM';...
	'Anterolateral visual','AL';...
	'Lateral visual','L';...
	'Rostrolateral area','RL';...
	'Anterior area','A';...
	'Subiculum','H-Sub';...
	'Field CA1','H-CA1';...
	'Field CA2','H-CA2';...
	'Field CA3','H-CA3';...
	'Dentate gyrus','H-DG';...
	'Retrosplenial','RetSpl';...
	};
strFigDir = 'F:\Data\Results\NpxDims\';


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
		if sum(indArea1Neurons) < 20, continue;end
		%% get orientation responses & single-trial population noise
		sArea1Neurons = sUseNeuron(indArea1Neurons);
		
		%% get spike times
		cellSpikeTimes1 = {sArea1Neurons.SpikeTimes};
		matSpikeCountsArea1 = getSpikeCounts(cellSpikeTimes1,vecStimOnTime,vecStimOffTime);
		
		%% MD area 1
		%real[vecNoiseParallel,vecNoiseOrthogonal,vecNoiseTotal]
		[vecNoiseParallel1,vecNoiseOrthogonal1,vecNoiseTotal1] = getNoiseParaOrtho(matSpikeCountsArea1,vecOrientation);
		% shuffled
		boolRandomFprime = true;
		matSpikeCountsShuffled1 = matSpikeCountsArea1(:,randperm(size(matSpikeCountsArea1,2)));
		[vecNoiseParallel1_S,vecNoiseOrthogonal1_S,vecNoiseTotal1_S] = getNoiseParaOrtho(matSpikeCountsArea1,vecOrientation,boolRandomFprime);
		
		%% select area 2
		for intArea2=(intArea1):intAreas
			strArea2 = cellUseAreas{intArea2}
			indArea2Neurons = contains({sUseNeuron.Area},strArea2,'IgnoreCase',true);
			if sum(indArea2Neurons) < 10, continue;end
			
			%% get orientation responses & single-trial population noise
			sArea2Neurons = sUseNeuron(indArea2Neurons);
			
			%% get spike times
			cellSpikeTimes2 = {sArea2Neurons.SpikeTimes};
			if intArea2 == intArea1
				%% area 1
				cellSpikeTimes1 = {sArea1Neurons.SpikeTimes};
				matSpikeCountsArea1 = getSpikeCounts(cellSpikeTimes1,vecStimOnTime,vecStimOffTime);
				intUseN = floor(size(matSpikeCountsArea1,1)/2);
				matSpikeCountsArea2 = matSpikeCountsArea1((intUseN+1):(intUseN*2),:);
				matSpikeCountsArea1 = matSpikeCountsArea1(1:intUseN,:);
				
				%real[vecNoiseParallel,vecNoiseOrthogonal,vecNoiseTotal]
				[vecNoiseParallel1,vecNoiseOrthogonal1,vecNoiseTotal1] = getNoiseParaOrtho(matSpikeCountsArea1,vecOrientation);
				% shuffled
				boolRandomFprime = true;
				matSpikeCountsShuffled1 = matSpikeCountsArea1(:,randperm(size(matSpikeCountsArea1,2)));
				[vecNoiseParallel1_S,vecNoiseOrthogonal1_S,vecNoiseTotal1_S] = getNoiseParaOrtho(matSpikeCountsArea1,vecOrientation,boolRandomFprime);
				
				% shuffle
				matSpikeCountsShuffled2 = matSpikeCountsArea2(:,randperm(size(matSpikeCountsArea2,2)));
			elseif intArea2 == (intArea1+1)
				%% area 1
				cellSpikeTimes1 = {sArea1Neurons.SpikeTimes};
				matSpikeCountsArea1 = getSpikeCounts(cellSpikeTimes1,vecStimOnTime,vecStimOffTime);
				
				%real[vecNoiseParallel,vecNoiseOrthogonal,vecNoiseTotal]
				[vecNoiseParallel1,vecNoiseOrthogonal1,vecNoiseTotal1] = getNoiseParaOrtho(matSpikeCountsArea1,vecOrientation);
				% shuffled
				boolRandomFprime = true;
				matSpikeCountsShuffled1 = matSpikeCountsArea1(:,randperm(size(matSpikeCountsArea1,2)));
				[vecNoiseParallel1_S,vecNoiseOrthogonal1_S,vecNoiseTotal1_S] = getNoiseParaOrtho(matSpikeCountsArea1,vecOrientation,boolRandomFprime);
				
				%area 2
				matSpikeCountsArea2 = getSpikeCounts(cellSpikeTimes2,vecStimOnTime,vecStimOffTime);
				% shuffle
				matSpikeCountsShuffled2 = matSpikeCountsArea2(:,randperm(size(matSpikeCountsArea2,2)));
			else
				%real
				matSpikeCountsArea2 = getSpikeCounts(cellSpikeTimes2,vecStimOnTime,vecStimOffTime);
				% shuffle
				matSpikeCountsShuffled2 = matSpikeCountsArea2(:,randperm(size(matSpikeCountsArea2,2)));
			end
			%% area 2
			%real
			[vecNoiseParallel2,vecNoiseOrthogonal2,vecNoiseTotal2] = getNoiseParaOrtho(matSpikeCountsArea2,vecOrientation);
			
			%matSpikeCountsShuffled2 = reshape(matSpikeCountsArea2(randperm(numel(matSpikeCountsArea2))),size(matSpikeCountsArea2));
			[vecNoiseParallel2_S,vecNoiseOrthogonal2_S,vecNoiseTotal2_S] = getNoiseParaOrtho(matSpikeCountsArea2,vecOrientation,boolRandomFprime);
			
			%% compare
			%real
			dblR_PP = corr(vecNoiseParallel1(:),vecNoiseParallel2(:));
			dblR_OP = corr(vecNoiseOrthogonal1(:),vecNoiseParallel2(:));
			dblR_OO = corr(vecNoiseOrthogonal1(:),vecNoiseOrthogonal2(:));
			dblR_PO = corr(vecNoiseParallel1(:),vecNoiseOrthogonal2(:));
			vecR = [dblR_PP dblR_OP dblR_OO dblR_PO];
			
			%shuffled
			dblR_PP_S = corr(vecNoiseParallel1_S(:),vecNoiseParallel2_S(:));
			dblR_OP_S = corr(vecNoiseOrthogonal1_S(:),vecNoiseParallel2_S(:));
			dblR_OO_S = corr(vecNoiseOrthogonal1_S(:),vecNoiseOrthogonal2_S(:));
			dblR_PO_S = corr(vecNoiseParallel1_S(:),vecNoiseOrthogonal2_S(:));
			vecR_S = [dblR_PP_S dblR_OP_S dblR_OO_S dblR_PO_S];
			
			%add to matrices
			matR_PP_All(intArea1,intArea2,intRec,1) = vecR(1);
			matR_OP_All(intArea1,intArea2,intRec,1) = vecR(2);
			matR_OO_All(intArea1,intArea2,intRec,1) = vecR(3);
			matR_PO_All(intArea1,intArea2,intRec,1) = vecR(4);
			
			matR_PP_All(intArea2,intArea1,intRec,1) = vecR(1);
			matR_OP_All(intArea2,intArea1,intRec,1) = vecR(2);
			matR_OO_All(intArea2,intArea1,intRec,1) = vecR(3);
			matR_PO_All(intArea2,intArea1,intRec,1) = vecR(4);
			
			%add to matrices
			matR_PP_All(intArea1,intArea2,intRec,2) = vecR_S(1);
			matR_OP_All(intArea1,intArea2,intRec,2) = vecR_S(2);
			matR_OO_All(intArea1,intArea2,intRec,2) = vecR_S(3);
			matR_PO_All(intArea1,intArea2,intRec,2) = vecR_S(4);
			
			matR_PP_All(intArea2,intArea1,intRec,2) = vecR_S(1);
			matR_OP_All(intArea2,intArea1,intRec,2) = vecR_S(2);
			matR_OO_All(intArea2,intArea1,intRec,2) = vecR_S(3);
			matR_PO_All(intArea2,intArea1,intRec,2) = vecR_S(4);
			
			
			%% pairwise cross-correlations
			%normal
			matNoiseCrossCorrs_All = getPairwiseNoiseCrossCorrelations(matSpikeCountsArea1,matSpikeCountsArea2,vecOrientation);
			matNoiseCrossCorrs = nanmean(matNoiseCrossCorrs_All,3);
			mat_xR_All(intArea2,intArea1,intRec,1) = nanmean(matNoiseCrossCorrs(:));
			
			%shuffled
			matNoiseCrossCorrs_All_S = getPairwiseNoiseCrossCorrelations(matSpikeCountsShuffled1,matSpikeCountsShuffled2,vecOrientation);
			matNoiseCrossCorrs_S = nanmean(matNoiseCrossCorrs_All_S,3);
			mat_xR_All(intArea2,intArea1,intRec,2) = nanmean(matNoiseCrossCorrs_S(:));
			
		end
	end
end
toc
%% average & remove empty rows
%shorten names
cellAreaAbbr = cell(size(cellUseAreas));
for intArea=1:numel(cellAreaAbbr)
	cellAreaAbbr{intArea} = cellRepStr{cellfun(@contains,cellfill(cellUseAreas{intArea},size(cellRepStr(:,1))),cellRepStr(:,1)),2};
end

%proc data
matR_PP = nanmean(matR_PP_All(:,:,:,1),3);
matR_OP = nanmean(matR_OP_All(:,:,:,1),3);
matR_PO = nanmean(matR_PO_All(:,:,:,1),3);
matR_OO = nanmean(matR_OO_All(:,:,:,1),3);
matR_PP_S = nanmean(matR_PP_All(:,:,:,2),3);
matR_OP_S = nanmean(matR_OP_All(:,:,:,2),3);
matR_PO_S = nanmean(matR_PO_All(:,:,:,2),3);
matR_OO_S = nanmean(matR_OO_All(:,:,:,2),3);

vecRemR = all(isnan(matR_PP),1);
vecRemC = all(isnan(matR_PP),2);

matR_PP(vecRemR,:) = [];
matR_PP(:,vecRemC) = [];
matR_OP(vecRemR,:) = [];
matR_OP(:,vecRemC) = [];
matR_OO(vecRemR,:) = [];
matR_OO(:,vecRemC) = [];
matR_PO(vecRemR,:) = [];
matR_PO(:,vecRemC) = [];

matR_PP_S(vecRemR,:) = [];
matR_PP_S(:,vecRemC) = [];
matR_OP_S(vecRemR,:) = [];
matR_OP_S(:,vecRemC) = [];
matR_OO_S(vecRemR,:) = [];
matR_OO_S(:,vecRemC) = [];
matR_PO_S(vecRemR,:) = [];
matR_PO_S(:,vecRemC) = [];

vecLimC = [-1 1];%max(flat(abs([matR_OP matR_PO])))*[-1 1];
figure
subplot(3,4,1)
imagesc(matR_PP,[-1 1]);
colormap(redbluepurple);
hcb = nancolorbar(matR_PP,[-1 1],redblue);
title(sprintf('Corr Parallel-parallel var'))
set(gca,'ytick',1:numel(cellAreaAbbr),'yticklabel',cellAreaAbbr)
set(gca,'xtick',1:numel(cellAreaAbbr),'xticklabel',cellAreaAbbr)
xtickangle(gca,45);

subplot(3,4,2)
imagesc(matR_OP,vecLimC);
colormap(redbluepurple);
hcb = nancolorbar(matR_OP,vecLimC,redblue);
title(sprintf('Corr Orthogonal-parallel var'))
set(gca,'ytick',1:numel(cellAreaAbbr),'yticklabel',cellAreaAbbr)
set(gca,'xtick',1:numel(cellAreaAbbr),'xticklabel',cellAreaAbbr)
xtickangle(gca,45);

subplot(3,4,3)
imagesc(matR_OO,[-1 1]);
colormap(redbluepurple);
hcb = nancolorbar(matR_OO,[-1 1],redblue);
title(sprintf('Corr Orthogonal-Orthogonal var'))
set(gca,'ytick',1:numel(cellAreaAbbr),'yticklabel',cellAreaAbbr)
set(gca,'xtick',1:numel(cellAreaAbbr),'xticklabel',cellAreaAbbr)
xtickangle(gca,45);

subplot(3,4,4)
imagesc(matR_PO,vecLimC);
colormap(redbluepurple);
hcb = nancolorbar(matR_PO,vecLimC,redblue);
title(sprintf('Corr Parallel-Orthogonal var'))
set(gca,'ytick',1:numel(cellAreaAbbr),'yticklabel',cellAreaAbbr)
set(gca,'xtick',1:numel(cellAreaAbbr),'xticklabel',cellAreaAbbr)
xtickangle(gca,45);

subplot(3,4,5)
imagesc(matR_PP_S,[-1 1]);
colormap(redbluepurple);
hcb = nancolorbar(matR_PP_S,[-1 1],redblue);
title(sprintf('Shuffled, Corr Para-para var'))
set(gca,'ytick',1:numel(cellAreaAbbr),'yticklabel',cellAreaAbbr)
set(gca,'xtick',1:numel(cellAreaAbbr),'xticklabel',cellAreaAbbr)
xtickangle(gca,45);

subplot(3,4,6)
imagesc(matR_OP_S,vecLimC);
colormap(redbluepurple);
hcb = nancolorbar(matR_OP_S,vecLimC,redblue);
title(sprintf('Shuffled, Corr orth-para var'))
set(gca,'ytick',1:numel(cellAreaAbbr),'yticklabel',cellAreaAbbr)
set(gca,'xtick',1:numel(cellAreaAbbr),'xticklabel',cellAreaAbbr)
xtickangle(gca,45);

subplot(3,4,7)
imagesc(matR_OO_S,[-1 1]);
colormap(redbluepurple);
hcb = nancolorbar(matR_OO_S,[-1 1],redblue);
title(sprintf('Shuffled, Corr orth-orth var'))
set(gca,'ytick',1:numel(cellAreaAbbr),'yticklabel',cellAreaAbbr)
set(gca,'xtick',1:numel(cellAreaAbbr),'xticklabel',cellAreaAbbr)
xtickangle(gca,45);

subplot(3,4,8)
imagesc(matR_PO_S,vecLimC);
colormap(redbluepurple);
hcb = nancolorbar(matR_PO_S,vecLimC,redblue);
title(sprintf('Shuffled, Corr Para-orth var'))
set(gca,'ytick',1:numel(cellAreaAbbr),'yticklabel',cellAreaAbbr)
set(gca,'xtick',1:numel(cellAreaAbbr),'xticklabel',cellAreaAbbr)
xtickangle(gca,45);

subplot(3,4,9)
imagesc(matR_PP - matR_PP_S,[-1 1]);
colormap(redbluepurple);
hcb = nancolorbar(matR_PP - matR_PP_S,[-1 1],redblue);
title(sprintf('Diff, Corr Parallel-parallel var'))
set(gca,'ytick',1:numel(cellAreaAbbr),'yticklabel',cellAreaAbbr)
set(gca,'xtick',1:numel(cellAreaAbbr),'xticklabel',cellAreaAbbr)
xtickangle(gca,45);

subplot(3,4,10)
imagesc(matR_OP - matR_OP_S,vecLimC);
colormap(redbluepurple);
hcb = nancolorbar(matR_OP - matR_OP_S,vecLimC,redblue);
title(sprintf('Diff, Corr Orthogonal-parallel var'))
set(gca,'ytick',1:numel(cellAreaAbbr),'yticklabel',cellAreaAbbr)
set(gca,'xtick',1:numel(cellAreaAbbr),'xticklabel',cellAreaAbbr)
xtickangle(gca,45);

subplot(3,4,11)
imagesc(matR_OO - matR_OO_S,[-1 1]);
colormap(redbluepurple);
hcb = nancolorbar(matR_OO - matR_OO_S,[-1 1],redblue);
title(sprintf('Diff, Corr Orthogonal-Orthogonal var'))
set(gca,'ytick',1:numel(cellAreaAbbr),'yticklabel',cellAreaAbbr)
set(gca,'xtick',1:numel(cellAreaAbbr),'xticklabel',cellAreaAbbr)
xtickangle(gca,45);

subplot(3,4,12)
imagesc(matR_PO - matR_PO_S,vecLimC);
colormap(redbluepurple);
hcb = nancolorbar(matR_PO - matR_PO_S,vecLimC,redblue);
title(sprintf('Diff, Corr Parallel-Orthogonal var'))
set(gca,'ytick',1:numel(cellAreaAbbr),'yticklabel',cellAreaAbbr)
set(gca,'xtick',1:numel(cellAreaAbbr),'xticklabel',cellAreaAbbr)
xtickangle(gca,45);

strFigFile = sprintf('CorrOrthPara_%s',getDate);
			maxfig;
			drawnow;
			export_fig([strFigDir strFigFile '.tif']);
			export_fig([strFigDir strFigFile '.pdf']);
			
%% cross corrs
matR_xR = nanmean(mat_xR_All(:,:,:,1),3);
matR_xR_S = nanmean(mat_xR_All(:,:,:,2),3);

figure
subplot(2,3,1)
imagesc(matR_xR,[-1 1]);
colormap(redbluepurple);
hcb = nancolorbar(matR_xR,[-1 1],redblue);
title(sprintf('xCorr'))

subplot(2,3,2)
imagesc(matR_xR_S,[-1 1]);
colormap(redbluepurple);
hcb = nancolorbar(matR_xR_S,[-1 1],redblue);
title(sprintf('xCorr Shuffled'))

%% means
intRecs = size(matR_PP_All,3);
matTril3 = repmat(tril(true(size(matR_PP_All,[1 2]))),[1 1 intRecs]);
matDiag3 = repmat(logical(eye(size(matR_PP_All,[1 2]))),[1 1 intRecs]);
matR_PP_D = matR_PP_All(:,:,:,1) - matR_PP_All(:,:,:,2);
matR_OP_D = matR_OP_All(:,:,:,1) - matR_OP_All(:,:,:,2);
matR_PO_D = matR_PO_All(:,:,:,1) - matR_PO_All(:,:,:,2);
matR_OO_D = matR_OO_All(:,:,:,1) - matR_OO_All(:,:,:,2);

vecSelfD_PP = matR_PP_D(matDiag3);
vecSelfD_OP = matR_OP_D(matDiag3);
vecSelfD_PO = matR_PO_D(matDiag3);
vecSelfD_OO = matR_OO_D(matDiag3);
vecSelfD_PP = vecSelfD_PP(~isnan(vecSelfD_PP));
vecSelfD_OP = vecSelfD_OP(~isnan(vecSelfD_OP));
vecSelfD_PO = vecSelfD_PO(~isnan(vecSelfD_PO));
vecSelfD_OO = vecSelfD_OO(~isnan(vecSelfD_OO));

vecCrossAreaD_PP = matR_PP_D(matTril3);
vecCrossAreaD_OP = matR_OP_D(matTril3);
vecCrossAreaD_PO = matR_PO_D(matTril3);
vecCrossAreaD_OO = matR_OO_D(matTril3);
vecCrossAreaD_PP = vecCrossAreaD_PP(~isnan(vecCrossAreaD_PP));
vecCrossAreaD_OP = vecCrossAreaD_OP(~isnan(vecCrossAreaD_OP));
vecCrossAreaD_PO = vecCrossAreaD_PO(~isnan(vecCrossAreaD_PO));
vecCrossAreaD_OO = vecCrossAreaD_OO(~isnan(vecCrossAreaD_OO));

%t-tests
[h,pSPP0] = ttest(vecSelfD_PP);
[h,pSOO0] = ttest(vecSelfD_OO);
[h,pSPO] = ttest(vecSelfD_PP,vecSelfD_OO);

[h,pCPP0] = ttest(vecCrossAreaD_PP);
[h,pCOO0] = ttest(vecCrossAreaD_OO);
[h,pCPO] = ttest(vecCrossAreaD_PP,vecCrossAreaD_OO);

%plot
figure
subplot(2,3,1)
hold on
%plot([0 5],[0 0],'color',0*[0.5 0.5 0.5]);
bplot(vecSelfD_PP,1,'outliers','linewidth',10);
bplot(vecSelfD_OO,2,'outliers','linewidth',10);
ylabel('dCorr var magnitude over shuff')
hold off
title(sprintf('r(Self);PP-0,p=%.3f;OO-0,p=%.3f;PP-OO,p=%.3f',pSPP0,pSOO0,pSPO))
set(gca,'xtick',1:2,'xticklabel',{'Para-para','Orth-orth'});
xtickangle(gca,45)
fixfig;
%grid off

subplot(2,3,2)
hold on
%plot([0 5],[0 0],'k--');
bplot(vecCrossAreaD_PP,1,'outliers','linewidth',10);
bplot(vecCrossAreaD_OO,2,'outliers','linewidth',10);
ylabel('dCorr var magnitude over shuff')
hold off
title(sprintf('r(Across);PP-0,p=%.3f;OO-0,p=%.3f;PP-OO,p=%.3f',pCPP0,pCOO0,pCPO))
set(gca,'xtick',1:2,'xticklabel',{'Para-para','Orth-orth'});
xtickangle(gca,45)
fixfig;
maxfig;

strFigFile = sprintf('MeanCorrOrthParaOverShuff_%s',getDate);
maxfig;
drawnow;
export_fig([strFigDir strFigFile '.tif']);
export_fig([strFigDir strFigFile '.pdf']);
		