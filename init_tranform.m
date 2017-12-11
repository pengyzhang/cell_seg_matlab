function tranform = init_tranform(u, dict)
nPhi = size(u,2);
for iPhi = 1:nPhi
    BW = u{iPhi}>0;
    BWbound = bwboundaries(u{iPhi}>0);
    BWbound = BWbound{1};
    corShapeX = BWbound(:,2);
    corShapeY = BWbound(:,1);
    nPoint = length(corShapeX);
    STATS = regionprops(BW,'Orientation','Centroid','MajorAxisLength','Perimeter');
    majorAxisAngle = STATS.Orientation*pi/180;
    startpointColumn = STATS.Centroid(1)+STATS.MajorAxisLength*cos(majorAxisAngle)/2; % Approximation of the intersection point between contour and the major axis, which is the startpoint of the shape
    startpointRow = STATS.Centroid(2)-STATS.MajorAxisLength*sin(majorAxisAngle)/2;
    
    while ~BW(floor(startpointRow), floor(startpointColumn)) % Refine the intersection point
        startpointColumn = startpointColumn - cos(majorAxisAngle)/2;
        startpointRow = startpointRow + sin(majorAxisAngle)/2;
            
    end
    distFromStartpoint = sqrt((corShapeX-startpointColumn*ones(nPoint,1)).^2+(corShapeY-startpointRow*ones(nPoint,1)).^2); % Sort shape points according to Euclidean dist to startpoint
    [unused,minPoint] = min(distFromStartpoint);
    corShapeX = [corShapeX;corShapeX]; % Duplicate corShapeX
    corShapeY = [corShapeY;corShapeY];
    corShapeX = corShapeX(minPoint:minPoint+nPoint-1); % Rearrange point on ellipse
    corShapeY = corShapeY(minPoint:minPoint+nPoint-1);
    segLength = 0;
    
    for kPoint = 1:nPoint

        if kPoint == 1
            segLength(kPoint) = sqrt((corShapeX(kPoint)-startpointColumn)^2+(corShapeY(kPoint)-startpointRow)^2)+epsilon;
        elseif kPoint <= nPoint
            segLength(kPoint) = sqrt((corShapeX(kPoint)-corShapeX(kPoint-1))^2+(corShapeY(kPoint)-corShapeY(kPoint-1))^2)+epsilon;
        end
    end

    segLength(nPoint+1) = sqrt((corShapeX(nPoint)-startpointColumn)^2+(corShapeY(nPoint)-startpointRow)^2)+epsilon;
    perimeter = sum(segLength);
    arcLength = cumsum(segLength,2)./perimeter; % Uniform arclength in [0,1]
    arcLength = [0,arcLength];
    arcGrid = 0:0.01:1;
    interpColumn = interp1(arcLength',[startpointColumn;corShapeX;startpointColumn],arcGrid','spline'); % Interpolation with 101 points on each shape
    interpRow = interp1(arcLength',[startpointRow;corShapeY;startpointRow],arcGrid','spline');
    interpColumn = interpColumn - STATS.Centroid(1); % Relative position
    interpRow = interpRow - STATS.Centroid(2);
        
    y = [interpColumn,interpRow];
    y = reshape(y,2*length(interpColumn),1);
    
end
end