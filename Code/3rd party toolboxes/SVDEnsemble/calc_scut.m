function [scut] = calc_scut(data, num_shuff, p)
% calculate threshold of similarity values from shuffled data

% num_shuff = 20;
% p = 0.98;
dims = size(data);

% calculate similarity matrix
warning('off')
S = zeros(dims(2),dims(2),num_shuff);
for n = 1:num_shuff
    data_shuff = shuffle(data,'time');
    S(:,:,n) = 1-pdist2(data_shuff',data_shuff','cosine');
end
warning('on')

% determine threshold
bins = 0:0.01:1;
cd = histc(S(:), bins);
cd = cumsum(cd / sum(cd));
scut = bins(find(cd > p, 1));

end
