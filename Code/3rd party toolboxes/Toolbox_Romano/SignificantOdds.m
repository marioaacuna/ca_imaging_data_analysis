function [mapOfOdds] = SignificantOdds(deltaFoF, sigma, densityData, densityNoise, xev)
        
    pCutOff=(100 - 95) / 100;
    mapOfOdds = densityNoise <= pCutOff * densityData;
    
%     if plotFlag
%         transfDataMatrix = deltaFoF  ./ sigma;
%         points = [transfDataMatrix(1:end-1,:), transfDataMatrix(2:end,:)];
%         points(any(isnan(points).'), :) = [];
%         % We plot the result
%         figure; plot(points(:,1),points(:,2),'k.'); axis equal; hold on
%         h=imagesc(xev(1,:),xev(1,:), mapOfOdds); axis xy; axis tight
%         alpha(h,0.6);
%         hxlab=xlabel('z-transformed value @ sample i'); hylab=ylabel('z-transformed value @ sample i+1'); set(gcf,'color','w');
%         set(gca,'FontSize',14); set([hxlab hylab],'FontSize',14)
%     end
end
