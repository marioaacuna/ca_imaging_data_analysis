function [ movie, ResultsOut ] = turboregCellMaps( p,movie )

    options = p.turboreg.options;

%% Get registration coordinates with turboreg
    
    refPic = single(movie(:,:,1));
    mask = zeros(size(refPic),'single'); 

    % run turboreg
    for i = 1:size(movie,3)
        pic = single(movie(:,:,i));
        [~,ResultsOut{i}] = turboreg(refPic,pic,mask,mask,options);
    end
    
    
%% Perform registration with transfturboreg
   
    mask = ones(size(movie(:,:,1)),'single');
    for i = 1:size(movie,3)
        pic = single(movie(:,:,i));
        movie(:,:,i) = transfturboreg(pic,mask,ResultsOut{i});
    end    

end

