function vecR = getAnalysisNpxMD(matSpikeCountsArea1,matSpikeCountsArea2,vecOrientation)
	%UNTITLED Summary of this function goes here
	%   Detailed explanation goes here
	%
	%"distinct subnetworks in mouse visual system for visual and non-visual signals"
	%
	%At single trial basis, split which info information limiting noise
	%(projected unto f') and non-limiting noise (orthogonal directions). Then
	%look at correlation of those noise types between areas. E.g., lgn diff
	%corrs correlate with v1, but SC non-diff corrs correlate with v1
	
	%% area 1
	%real
	[vecNoiseParallel1,vecNoiseOrthogonal1] = getNoiseParaOrtho(matSpikeCountsArea1,vecOrientation);
	
	%% area 2
	%real
	[vecNoiseParallel2,vecNoiseOrthogonal2] = getNoiseParaOrtho(matSpikeCountsArea2,vecOrientation);
	
	%% compare
	dblR_PP = corr(vecNoiseParallel1(:),vecNoiseParallel2(:));
	dblR_OP = corr(vecNoiseOrthogonal1(:),vecNoiseParallel2(:));
	dblR_OO = corr(vecNoiseOrthogonal1(:),vecNoiseOrthogonal2(:));
	dblR_PO = corr(vecNoiseParallel1(:),vecNoiseOrthogonal2(:));
	vecR = [dblR_PP dblR_OP dblR_OO dblR_PO];
	
end
%{
%% plot
	subplot(2,3,1)
	histx(vecNoiseParallel)
	mean(vecNoiseParallel)
	
	subplot(2,3,2)
	histx(vecNoiseOrthogonal)
	mean(vecNoiseOrthogonal)
	
	subplot(2,3,4)
	histx(vecNoiseParallelShuff)
	mean(vecNoiseParallelShuff)
	
	subplot(2,3,5)
	histx(vecNoiseOrthogonalShuff)
	mean(vecNoiseOrthogonalShuff)

	return
	%%
	vecDims = [3 4];
	vecRef2 = vecRef*10;
	cla;
	hold on
	plot3([0 vecRef2(vecDims(1))],[0 vecRef2(vecDims(2))],[0 vecRef2(vecDims(3))]);
	scatter3(vecPoint(vecDims(1)),vecPoint(vecDims(2)),vecPoint(vecDims(3)))
	scatter3(vecParallel(vecDims(1)),vecParallel(vecDims(2)),vecParallel(vecDims(3)))
	hold off
	
	%%
	vecDims = [3 4];
	
	%stim1
		vecTrialsOri1 = find(vecOriIdx==intOriIdx);
		intRepNum = numel(vecTrialsOri1);
		matAct1 = matSpikeCountsArea1(vecDims,vecTrialsOri1);
		vecMu1 = mean(matAct1,2);
		
		%stim2
		intOriIdx2 = modx(intOriIdx+1,intOriNum);
		matAct2 = matSpikeCountsArea1(vecDims,vecOriIdx==intOriIdx2);
		vecMu2 = mean(matAct2,2);
		
		%get ref & project
		vecRef = vecMu2 - vecMu1;
		vecRef = vecRef / norm(vecRef);
		matRef = matSpikeCountsArea1(vecDims,:)  - vecMu1;
		
		vecNoiseParallel = nan(1,intRepNum);
		vecNoiseOrthogonal = nan(1,intRepNum);
		
		intTrial = vecTrialsOri1(intRep);
		
			vecPoint = matRef(:,intTrial);
			
			%project vecPoint onto vecRef and decompose in parallel & orthogonal
			vecParallel = ((vecRef'*vecPoint)/(norm(vecRef).^2)).*vecRef;
			vecOrtho = vecPoint - vecParallel;
			
			
			vecRef2 = vecRef*10;
	cla;
	hold on
	plot([0 vecRef2(1)],[0 vecRef2(2)]);
	scatter(vecPoint(1),vecPoint(2),'ok')
	scatter(vecParallel(1),vecParallel(2),'rx')
	scatter(vecOrtho(1),vecOrtho(2),'bx')
	hold off
	xlim([-2 10])
	ylim([-2 10])
%}

