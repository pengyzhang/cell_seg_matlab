%distance map
function [d randomMat] = distance_map(Img, dict, transform)
d = [];
[dim,nSample] = size(dict);
[nRow,nCol] = size(Img);
for iSample = 1:nSample
    oddIndex = 1:2:dim;
    evenIndex = 2:2:dim;
    X = dict(oddIndex,iSample);
    Y = dict(evenIndex,iSample);
    xyCordinate = [X, Y];
    xyCordinate = transform.b*xyCordinate*transform.T+transform.c;
    X = xyCordinate(:,1);
    Y = xyCordinate(:,2);
    BW = roipoly(Img, X, Y);
    edgeBW = edge(BW);
    dm = -bwdist(edgeBW);
    dm(BW==1) = -1*dm(BW==1);
    dm = dm/norm(dm,'fro');
    d{iSample} = dm;
%     figure,imshow(dm);
%     dm = reshape(dm,nRow*nCol,1);
%     dm = dm/norm(dm);
%     d = [d dm];
end

downDim = 50;
randomMat = normrnd(0,1,[downDim nRow*nCol]);
for i = 1:downDim
    randomMat(i,:) = randomMat(i,:)/norm(randomMat(1,:));
end
% d = randomMat*d;