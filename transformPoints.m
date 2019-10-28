function [outData] = transformPoints(inpData)
    
    ccwArea = 0;
    
    cwSum = 0;
    
    scaleFactor = 5;

    centroid = [0, 0];
    sArea = 0;
    
    tXPos = -1e9;
    tIndex = 0;
    
    for i = 1:size(inpData, 1)
        
        target = i+1;
        if target > size(inpData, 1)
            target = 1;
        end
        
        sPoint = inpData(i, :);
        tPoint = inpData(target, :);
        
        
        cwSum = cwSum + (tPoint(1) - sPoint(1))*(tPoint(2) + sPoint(2));
        
        
        s = (sPoint(1) * tPoint(2)) - (sPoint(2) * tPoint(1));
        
        if abs(s) < 1e-2
            s = 1e-2;
        end
        
        sArea = sArea + s;
        
        sx = sPoint(1) + tPoint(1);
        sy = sPoint(2) + tPoint(2);
        
        centroid = centroid + ( [sx, sy]*s );       
        
        if sPoint(1) > tXPos
            tXPos = sPoint(1);
            tIndex = i;
        end
        
    end
    
    if cwSum < 0
        inpData = flipud(inpData);
    end
    
    sArea = sArea * 0.5;
    centroid = centroid / (6*sArea);
    
    inpData = circshift(inpData, -(tIndex-1), 1);
    
    for i = 1:size(inpData, 1)
        inpData(i, 1) = inpData(i, 1) - centroid(1);
        inpData(i, 2) = inpData(i, 2) - centroid(2);
    end
    
    mag = zeros(size(inpData, 1), 1);
    for i = 1:size(inpData, 1)
        mag(i, 1) = norm(inpData(i, :));
    end
    inpData = (inpData/max(mag))*scaleFactor;
    
    outData = inpData;    
end