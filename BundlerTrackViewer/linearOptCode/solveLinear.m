function result = solveLinear(noImages, alignmentsFile, cycleFile, newAlignmentsFile)

    fid = fopen(alignmentsFile, 'r');

    %the assumptions in the next two lines are made just to ease the memory allocation process
    maxPlanesPerImage = 2; %we assume there are at most 2 facade planes in each image
    maxGridAlignmentsPerPair = 5; %we assume there are at most 5 different 3D grids

    noAlignments = 0;
    noAlignments = fscanf(fid, '%d', 1);

    alignmentIds = zeros(noImages, noImages);
    alignmentImages = [];
    rowShifts = zeros(noImages*maxGridAlignmentsPerPair, noImages);
    colShifts = zeros(noImages*maxGridAlignmentsPerPair, noImages);
    weights = [];
    noGridAlignmentsPerPair = zeros(noImages, noImages);
    gridIdsInFirstImage = zeros(noImages*maxGridAlignmentsPerPair, noImages);
    gridIdsInSecondImage = zeros(noImages*maxGridAlignmentsPerPair, noImages);
    rotationMatrices = zeros(noImages*9, noImages);
    inconsistentCycles = [];
    inconsistencyMeasure = [];
    consistentCycles = [];
    noInconsistentCyclesInvolved = zeros(noAlignments, 1);

    imgIndex1 = 0;
    imgIndex2 = 0;
    imgIndex3 = 0;

    for i=1:noAlignments
        
        x = fscanf(fid, '%d', 2);
        imgIndex1 = x(1);
        imgIndex2 = x(2);
        display(['alg:' num2str(i) ' img1:' num2str(imgIndex1) 'img2:' num2str(imgIndex2)]);
        weight = fscanf(fid, '%f', 1);
        gridMatch = fscanf(fid, '%d', 1);
        alignmentIds(imgIndex1+1, imgIndex2+1) = i;
        alignmentIds(imgIndex2+1, imgIndex1+1) = i;
        alignmentImages = [alignmentImages; imgIndex1+1 imgIndex2+1];
        weights = [weights; weight];
        
        x = fscanf(fid, '%f', 3);
        x = fscanf(fid, '%d', 1);
        noShifts = x(1);
        noGridAlignmentsPerPair(imgIndex1+1, imgIndex2+1) =  noShifts;
        noGridAlignmentsPerPair(imgIndex2+1, imgIndex1+1) =  noShifts;
        for j=1:noShifts
            x = fscanf(fid, '%d', 7);
            display(['read' num2str(x(1)) ' ' num2str(x(2)) ' ' num2str(x(3)) ' ' num2str(x(4)) ' ' num2str(x(5)) ' ' num2str(x(6)) ' ' num2str(x(7))]);
            templateId = x(1);
            plane1 = x(2); plane2 = x(3);
            grid1 = x(4); grid2 = x(5);
            
            gridIdsInFirstImage(imgIndex1*5+j, imgIndex2+1) = templateId*maxGridAlignmentsPerPair+grid1;
            gridIdsInSecondImage(imgIndex1*5+j, imgIndex2+1) = templateId*maxGridAlignmentsPerPair+grid2;
            
            gridIdsInFirstImage(imgIndex2*5+j, imgIndex1+1) = templateId*maxGridAlignmentsPerPair+grid2;
            gridIdsInSecondImage(imgIndex2*5+j, imgIndex1+1) = templateId*maxGridAlignmentsPerPair+grid1;
            
            colShifts(imgIndex1*5+j, imgIndex2+1) = x(6);
            rowShifts(imgIndex1*5+j, imgIndex2+1) = x(7);
            
            colShifts(imgIndex2*5+j, imgIndex1+1) = -x(6);
            rowShifts(imgIndex2*5+j, imgIndex1+1) = -x(7);
        end
        x = fscanf(fid, '%f', 9);
        rotationMatrices(imgIndex1*9+1, imgIndex2+1) = x(1);
        rotationMatrices(imgIndex1*9+2, imgIndex2+1) = x(2);
        rotationMatrices(imgIndex1*9+3, imgIndex2+1) = x(3);
        rotationMatrices(imgIndex1*9+4, imgIndex2+1) = x(4);
        rotationMatrices(imgIndex1*9+5, imgIndex2+1) = x(5);
        rotationMatrices(imgIndex1*9+6, imgIndex2+1) = x(6);
        rotationMatrices(imgIndex1*9+7, imgIndex2+1) = x(7);
        rotationMatrices(imgIndex1*9+8, imgIndex2+1) = x(8);
        rotationMatrices(imgIndex1*9+9, imgIndex2+1) = x(9);
        R = [x(1) x(2) x(3); x(4) x(5) x(6); x(7) x(8) x(9)];
        R = R';
        rotationMatrices(imgIndex2*9+1, imgIndex1+1) = R(1,1);
        rotationMatrices(imgIndex2*9+2, imgIndex1+1) = R(1,2);
        rotationMatrices(imgIndex2*9+3, imgIndex1+1) = R(1,3);
        rotationMatrices(imgIndex2*9+4, imgIndex1+1) = R(2,1);
        rotationMatrices(imgIndex2*9+5, imgIndex1+1) = R(2,2);
        rotationMatrices(imgIndex2*9+6, imgIndex1+1) = R(2,3);
        rotationMatrices(imgIndex2*9+7, imgIndex1+1) = R(3,1);
        rotationMatrices(imgIndex2*9+8, imgIndex1+1) = R(3,2);
        rotationMatrices(imgIndex2*9+9, imgIndex1+1) = R(3,3);
    end
    
    cycles = [];
    fid = fopen(cycleFile, 'r');
    while(~feof(fid))
        x = fscanf(fid, '%d %d %d ', 3);
        y = sort(x);
        imgIndex1 = y(1);
        imgIndex2 = y(2);
        imgIndex3 = y(3);
        cycles = [cycles; imgIndex1 imgIndex2 imgIndex3];
    end
    
    [noCycles a] = size(cycles);

    fixedEdges = [];
    for i=1:noAlignments
        fixedEdges = [fixedEdges; i];
    end
    noFixedEdges = 0;
    tmpWeights = weights;

    tmpInconsistentCycles = inconsistentCycles;
    noChanged = 1;
    
    for iter=1:30
        for i=1:noAlignments
            if(noInconsistentCyclesInvolved(i) == 0)
                tmpWeights(i) = 1;
            else
                tmpWeights(i) = 1.0 / noInconsistentCyclesInvolved(i);
            end
        end

        %weights = weights/max(weights);
        
        noChosenEdges = (noAlignments-noFixedEdges) / 20 + 1;
        
        cvx_begin
            variable x(noAlignments-noFixedEdges);
            l = ones(noAlignments-noFixedEdges,1);
            l = l*0.1;
            u = ones(noAlignments-noFixedEdges,1);
    
            minimize(weights'*x);

            subject to
                l <= x <= u;
        
                %%%check for inconsistent cycles
                for i=1:noCycles
                     img1 = cycles(i, 1);
                     img2 = cycles(i, 2);
                     img3 = cycles(i, 3);
                     
                     gridIds = [];
                     cols = [];
                     rows = [];
                     
                     %if at least one pair does not have grid alignments, check for rotation consistency
                     if(noGridAlignmentsPerPair(img1+1, img2+1) == 0 || noGridAlignmentsPerPair(img2+1, img3+1) == 0 || noGridAlignmentsPerPair(img3+1, img1+1) == 0)
                         R1 = [rotationMatrices(img1*9+1,img2+1) rotationMatrices(img1*9+2,img2+1) rotationMatrices(img1*9+3,img2+1); rotationMatrices(img1*9+4,img2+1) rotationMatrices(img1*9+5,img2+1) rotationMatrices(img1*9+6,img2+1); rotationMatrices(img1*9+7,img2+1) rotationMatrices(img1*9+8,img2+1) rotationMatrices(img1*9+9,img2+1)];
                         R2 = [rotationMatrices(img2*9+1,img3+1) rotationMatrices(img2*9+2,img3+1) rotationMatrices(img2*9+3,img3+1); rotationMatrices(img2*9+4,img3+1) rotationMatrices(img2*9+5,img3+1) rotationMatrices(img2*9+6,img3+1); rotationMatrices(img2*9+7,img3+1) rotationMatrices(img2*9+8,img3+1) rotationMatrices(img2*9+9,img3+1)];
                         R3 = [rotationMatrices(img3*9+1,img1+1) rotationMatrices(img3*9+2,img1+1) rotationMatrices(img3*9+3,img1+1); rotationMatrices(img3*9+4,img1+1) rotationMatrices(img3*9+5,img1+1) rotationMatrices(img3*9+6,img1+1); rotationMatrices(img3*9+7,img1+1) rotationMatrices(img3*9+8,img1+1) rotationMatrices(img3*9+9,img1+1)];
                         R = R1*R2*R3;
                         cos_theta = 0.5 * (R(1,1) + R(2,2) + R(3,3) - 1.0);
                         cos_theta = max(-1.0, min(1.0, cos_theta));
                         angle = acos(cos_theta);
                         angleDeg = angle*180.0/pi;
                         if(angleDeg > 5)
                             display(['******incons:' num2str(img1) ',' num2str(img2) ',' num2str(img3), 'angle:' num2str(angleDeg)]);
                            x(alignmentIds(img1+1,img2+1)) + x(alignmentIds(img2+1,img3+1)) + x(alignmentIds(img3+1,img1+1)) >= 1.0
                         end
                     else
                        for a=1:noGridAlignmentsPerPair(img1+1, img2+1)
                            gridIds = [gridIds; gridIdsInSecondImage((img1)*5+a, img2+1)];
                            cols = [cols; colShifts((img1)*5+a, img2+1)];
                            rows = [rows; rowShifts((img1)*5+a, img2+1)];
                            %display(['grid ' num2str(gridIds(length(gridIds))) ':(' num2str(cols(length(cols))) ',' num2str(rows(length(rows))) '),']);
                        end
            
                        tmpGridIds = [];
                        tmpCols = [];
                        tmpRows = [];
                        for a=1:noGridAlignmentsPerPair(img2+1, img3+1)
                            grid = gridIdsInFirstImage((img2)*5+a, img3+1);
                            found = 0;
                            for b=1:length(gridIds)
                                if(grid == gridIds(b))
                                    found = b;
                                break;
                                end
                            end
                            if(found > 0)
                                tmpGridIds = [tmpGridIds; gridIdsInSecondImage((img2)*5+a, img3+1)];
                                tmpCols = [tmpCols; cols(found) + colShifts((img2)*5+a, img3+1)];
                                tmpRows = [tmpRows; rows(found) + rowShifts((img2)*5+a, img3+1)];
                                %display(['grid ' num2str(tmpGridIds(length(tmpGridIds))) ':(' num2str(tmpCols(length(tmpCols))) ',' num2str(tmpRows(length(tmpRows))) '),']);
                            end
                        end
                        cols = tmpCols;
                        rows = tmpRows;
                        gridIds = tmpGridIds;
                
                        tmpGridIds = [];
                        tmpCols = [];
                        tmpRows = [];
                        for a=1:noGridAlignmentsPerPair(img3+1, img1+1)
                            grid = gridIdsInFirstImage((img3)*5+a, img1+1);
                            found = 0;
                            for b=1:length(gridIds)
                                if(grid == gridIds(b))
                                    found = b;
                                    break;
                                end
                            end
                            if(found > 0)
                                tmpGridIds = [tmpGridIds; gridIdsInSecondImage((img3)*5+a, img1+1)];
                                tmpCols = [tmpCols; cols(found) + colShifts((img3)*5+a, img1+1)];
                                tmpRows = [tmpRows; rows(found) + rowShifts((img3)*5+a, img1+1)];
                                %display(['grid ' num2str(tmpGridIds(length(tmpGridIds))) ':(' num2str(tmpCols(length(tmpCols))) ',' num2str(tmpRows(length(tmpRows))) '),']);
                            end
                        end
                        cols = tmpCols;
                        rows = tmpRows;
                        gridIds = tmpGridIds;
            
                        for a=1:length(gridIds)
                            if(cols(a) ~= 0 || rows(a) ~= 0)
                                display(['******incons:' num2str(img1) ',' num2str(img2) ',' num2str(img3)]);
                                x(alignmentIds(img1+1,img2+1)) + x(alignmentIds(img2+1,img3+1)) + x(alignmentIds(img3+1,img1+1)) >= 1.0
                                break;
                            end
                        end
                     end
               end
        cvx_end
        
        for i=0:noImages-1
            for j=0:noImages-1
                if(alignmentIds(i+1, j+1) ~= 0)
                    display(['alignment from ' num2str(i) ' to ' num2str(j) ' with var: ' num2str(x(alignmentIds(i+1, j+1)))]); 
                    for a=1:noGridAlignmentsPerPair(i+1, j+1)
                        display(['grid ' num2str(gridIdsInFirstImage(i*5+a, j+1)) ' to ' num2str(gridIdsInSecondImage(i*5+a, j+1)) ':(' num2str(colShifts(i*5+a, j+1)) ',' num2str(rowShifts(i*5+a, j+1)) ')']);
                    end
                end
            end
        end
        
        %after the optimization, for every two node, find the shortest path and compare the cost of the path with
        %the cost of the edge connecting the two nodes. If necessary update the alignment of the edge with the composite
        %alignment along the path
        noChanged = 0;
        newColShifts = colShifts;
        newRowShifts = rowShifts;
        newWeights = weights;
        newNoGridAlignmentsPerPair = noGridAlignmentsPerPair;
        newGridIdsInFirstImage = gridIdsInFirstImage;
        newGridIdsInSecondImage = gridIdsInSecondImage;
        newRotationMatrices = rotationMatrices;
        for i=1:noImages
            distances = ones(noImages,1);
            distances = distances * (inf);
            processed = zeros(noImages,1);
            predecessors = zeros(noImages, 1);
            distances(i, 1) = 0.0;
            count = 0;
            while(count < noImages)
                %pick next node
                nextNode = 0;
                minDistance = inf;
                for j=1:noImages
                    if(processed(j,1) == 1)
                        continue;
                    end
                    if(minDistance >= distances(j, 1))
                        minDistance = distances(j, 1);
                        nextNode = j;
                    end
                end
                
                processed(nextNode, 1) = 1;
                count = count + 1;
                for j=1:noImages
                    if(processed(j,1) == 0 && alignmentIds(nextNode, j) ~= 0)
                        if(distances(j,1) > distances(nextNode,1) + x(alignmentIds(nextNode, j))/weights(alignmentIds(nextNode, j)))
                            distances(j, 1) = distances(nextNode,1) + x(alignmentIds(nextNode, j))/weights(alignmentIds(nextNode, j));
                            predecessors(j) = nextNode;
                        end
                    end
                end
            end
            
            for j=i+1:noImages
                if(alignmentIds(i, j) == 0)
                    continue;
                end
                
                currentNode = j;
                valid = 1;
                while(predecessors(currentNode) ~= 0)
                    if(x(alignmentIds(currentNode, predecessors(currentNode))) > 0.2)
                        valid = 0;
                        break;
                    end
                    currentNode = predecessors(currentNode);
                end
                    
                if(valid == 0)
                    continue;
                end
                    
                if(distances(j) < x(alignmentIds(i, j))/weights(alignmentIds(i, j)))
                    
                    
                    display(['changing alignment from ' num2str(i-1) ' to ' num2str(j-1)]);
                    
                    display(['old weight:' num2str(weights(alignmentIds(i,j))) ' old var:' num2str(x(alignmentIds(i,j)))]);
                    
                    %no grid alignment, just change the rotation matrices
                    if(noGridAlignmentsPerPair(i,j) == 0)
                        currentNode = j;
                        weightSum = 1.0;
                        R = [1 0 0; 0 1 0; 0 0 1];
                        while(predecessors(currentNode) ~= 0)
                            img1 = currentNode-1;
                            img2 = predecessors(currentNode);
                            tmpR = [rotationMatrices(img1*9+1,img2) rotationMatrices(img1*9+2,img2) rotationMatrices(img1*9+3,img2); rotationMatrices(img1*9+4,img2) rotationMatrices(img1*9+5,img2) rotationMatrices(img1*9+6,img2); rotationMatrices(img1*9+7,img2) rotationMatrices(img1*9+8,img2) rotationMatrices(img1*9+9,img2)];
                            R = R*tmpR;
                            weightSum = weightSum * weights(alignmentIds(currentNode, predecessors(currentNode)));
                            currentNode = predecessors(currentNode);
                        end
                        newRotationMatrices((j-1)*9+1,i) = R(1,1); newRotationMatrices((j-1)*9+2,i) = R(1,2); newRotationMatrices((j-1)*9+3,i) = R(1,3);
                        newRotationMatrices((j-1)*9+4,i) = R(2,1); newRotationMatrices((j-1)*9+5,i) = R(2,2); newRotationMatrices((j-1)*9+6,i) = R(2,3);
                        newRotationMatrices((j-1)*9+7,i) = R(3,1); newRotationMatrices((j-1)*9+8,i) = R(3,2); newRotationMatrices((j-1)*9+9,i) = R(3,3);
                        R = R';
                        newRotationMatrices((i-1)*9+1,j) = R(1,1); newRotationMatrices((i-1)*9+2,j) = R(1,2); newRotationMatrices((i-1)*9+3,j) = R(1,3);
                        newRotationMatrices((i-1)*9+4,j) = R(2,1); newRotationMatrices((i-1)*9+5,j) = R(2,2); newRotationMatrices((i-1)*9+6,j) = R(2,3);
                        newRotationMatrices((i-1)*9+7,j) = R(3,1); newRotationMatrices((i-1)*9+8,j) = R(3,2); newRotationMatrices((i-1)*9+9,j) = R(3,3);
                        if(weightSum > newWeights(alignmentIds(i,j)))
                            newWeights(alignmentIds(i,j)) = weightSum;% + newWeights(alignmentIds(i,j));
                        end
                    else    
                    	%change also grid alignments
                        for a=1:noGridAlignmentsPerPair(i,j)
                            display(['old alignment:(' num2str(colShifts((i-1)*5+a,j)) ',' num2str(rowShifts((i-1)*5+a,j)) ')']);
                        end
                    
                        currentNode = j;
                        weightSum = 1.0;
                        varSum = 0;
                        gridIds = [];
                        cols = [];
                        rows = [];
                        index = 0;
                        while(predecessors(currentNode) ~= 0)
                            display(['path from ' num2str(currentNode-1) ' to ' num2str(predecessors(currentNode)-1) ' weight: ' num2str(weights(alignmentIds(currentNode, predecessors(currentNode)))) ' var:' num2str(x(alignmentIds(currentNode, predecessors(currentNode)))) ]);
                            if(index == 0)
                                index = 1;
                                for a=1:noGridAlignmentsPerPair(currentNode, predecessors(currentNode))
                                    gridIds = [gridIds; gridIdsInSecondImage((currentNode-1)*5+a, predecessors(currentNode))];
                                    cols = [cols; -colShifts((currentNode-1)*5+a, predecessors(currentNode))];
                                    rows = [rows; -rowShifts((currentNode-1)*5+a, predecessors(currentNode))];
                                end
                            else
                                tmpGridIds = [];
                                tmpCols = [];
                                tmpRows = [];
                                for a=1:noGridAlignmentsPerPair(currentNode, predecessors(currentNode))
                                    grid = gridIdsInFirstImage((currentNode-1)*5+a, predecessors(currentNode));
                                    found = 0;
                                    for b=1:length(gridIds)
                                        if(grid == gridIds(b))
                                            found = b;
                                            break;
                                        end
                                    end
                                    if(found > 0)
                                        tmpGridIds = [tmpGridIds; gridIdsInSecondImage((currentNode-1)*5+a, predecessors(currentNode))];
                                        tmpCols = [tmpCols; cols(found) - colShifts((currentNode-1)*5+a, predecessors(currentNode))];
                                        tmpRows = [tmpRows; rows(found) - rowShifts((currentNode-1)*5+a, predecessors(currentNode))];
                                    end
                                end
                                cols = tmpCols;
                                rows = tmpRows;
                                gridIds = tmpGridIds;
                            end
                            weightSum = weightSum * weights(alignmentIds(currentNode, predecessors(currentNode)));% / x(alignmentIds(currentNode, predecessors(currentNode)));
                            varSum = varSum + 1.0 / x(alignmentIds(currentNode, predecessors(currentNode)));
                            currentNode = predecessors(currentNode);
                        end
                    
                        count = 0;
                        changed = 0;
                        for a=1:noGridAlignmentsPerPair(i, j)
                            grid = gridIdsInFirstImage((i-1)*5+a,j);
                            found = 0;
                            for b=1:length(gridIds)
                                if(grid == gridIds(b))
                                    found = b;
                                    break;
                                end
                            end
                            if(found > 0)
                                count = count + 1;
                                newGridIdsInFirstImage((i-1)*5+count, j) =  gridIdsInFirstImage((i-1)*5+a,j);
                                newGridIdsInSecondImage((i-1)*5+count, j) = gridIdsInSecondImage((i-1)*5+a,j);
                          
                                newGridIdsInFirstImage((j-1)*5+count, i) =  newGridIdsInSecondImage((i-1)*5+count, j);
                                newGridIdsInSecondImage((j-1)*5+count, i) = newGridIdsInFirstImage((i-1)*5+count, j);
                          
                                newColShifts((i-1)*5+count, j) = cols(found);
                                newRowShifts((i-1)*5+count, j) = rows(found);
                                newColShifts((j-1)*5+count, i) = -cols(found);
                                newRowShifts((j-1)*5+count, i) = -rows(found);
                          
                                display(['changing alignment from ' num2str(newGridIdsInFirstImage((i-1)*5+count, j)) ' to ' num2str(newGridIdsInSecondImage((i-1)*5+count, j)) ':(' num2str(newColShifts((i-1)*5+count, j)) ' ' num2str(newRowShifts((i-1)*5+count, j)) ')']);
                                if(newColShifts((i-1)*5+count, j) ~= colShifts((i-1)*5+a, j) || newRowShifts((i-1)*5+count, j) ~= rowShifts((i-1)*5+a, j))
                                    changed = 1;
                                end;
                            end
                        end
                        newNoGridAlignmentsPerPair(i,j) = count;
                        newNoGridAlignmentsPerPair(j,i) = count;
                    
                        if(changed)
                            newWeights(alignmentIds(i,j)) = weightSum;
                        else
                            if(weightSum > newWeights(alignmentIds(i,j)))
                                newWeights(alignmentIds(i,j)) = weightSum;% + newWeights(alignmentIds(i,j));
                            end
                        end
                        display(['new weight:' num2str(newWeights(alignmentIds(i,j)))]);
                        if(weightSum > 1)
                            int debug = 1;
                        end
                        if(newWeights(alignmentIds(i,j)) ~= weights(alignmentIds(i,j)) || changed == 1)
                            noChanged = noChanged + 1;
                        end
                    end
                end
            end
        end
        
        weights = newWeights;
        colShifts = newColShifts;
        rowShifts = newRowShifts;
        noGridAlignmentsPerPair = newNoGridAlignmentsPerPair;
        gridIdsInFirstImage = newGridIdsInFirstImage;
        gridIdsInSecondImage = newGridIdsInSecondImage;
        rotationMatrices = newRotationMatrices;
        
        %find inconsistent cycles
        noInconsistentCyclesInvolved = zeros(noAlignments, 1);
        for i=1:noCycles
             img1 = cycles(i, 1);
             img2 = cycles(i, 2);
             img3 = cycles(i, 3);

             gridIds = [];
             cols = [];
             rows = [];

             if(noGridAlignmentsPerPair(img1+1, img2+1) == 0 || noGridAlignmentsPerPair(img2+1, img3+1) == 0 || noGridAlignmentsPerPair(img3+1, img1+1) == 0)
                 R1 = [rotationMatrices(img1*9+1,img2+1) rotationMatrices(img1*9+2,img2+1) rotationMatrices(img1*9+3,img2+1); rotationMatrices(img1*9+4,img2+1) rotationMatrices(img1*9+5,img2+1) rotationMatrices(img1*9+6,img2+1); rotationMatrices(img1*9+7,img2+1) rotationMatrices(img1*9+8,img2+1) rotationMatrices(img1*9+9,img2+1)];
                 R2 = [rotationMatrices(img2*9+1,img3+1) rotationMatrices(img2*9+2,img3+1) rotationMatrices(img2*9+3,img3+1); rotationMatrices(img2*9+4,img3+1) rotationMatrices(img2*9+5,img3+1) rotationMatrices(img2*9+6,img3+1); rotationMatrices(img2*9+7,img3+1) rotationMatrices(img2*9+8,img3+1) rotationMatrices(img2*9+9,img3+1)];
                 R3 = [rotationMatrices(img3*9+1,img1+1) rotationMatrices(img3*9+2,img1+1) rotationMatrices(img3*9+3,img1+1); rotationMatrices(img3*9+4,img1+1) rotationMatrices(img3*9+5,img1+1) rotationMatrices(img3*9+6,img1+1); rotationMatrices(img3*9+7,img1+1) rotationMatrices(img3*9+8,img1+1) rotationMatrices(img3*9+9,img1+1)];
                 R = R1*R2*R3;
                 cos_theta = 0.5 * (R(1,1) + R(2,2) + R(3,3) - 1.0);
                 cos_theta = max(-1.0, min(1.0, cos_theta));
                 angle = acos(cos_theta);
                 angleDeg = angle*180.0/pi;
                 if(angleDeg > 5)
                     noInconsistentCyclesInvolved(alignmentIds(img1+1,img2+1)) = noInconsistentCyclesInvolved(alignmentIds(img1+1,img2+1)) + 1;
                     noInconsistentCyclesInvolved(alignmentIds(img2+1,img3+1)) = noInconsistentCyclesInvolved(alignmentIds(img2+1,img3+1)) + 1;
                     noInconsistentCyclesInvolved(alignmentIds(img3+1,img1+1)) = noInconsistentCyclesInvolved(alignmentIds(img3+1,img1+1)) + 1;
                 end
             else
                for a=1:noGridAlignmentsPerPair(img1+1, img2+1)
                    gridIds = [gridIds; gridIdsInSecondImage((img1)*5+a, img2+1)];
                    cols = [cols; colShifts((img1)*5+a, img2+1)];
                    rows = [rows; rowShifts((img1)*5+a, img2+1)];
                    %display(['grid ' num2str(gridIds(length(gridIds))) ':(' num2str(cols(length(cols))) ',' num2str(rows(length(rows))) '),']);
                end

                tmpGridIds = [];
                tmpCols = [];
                tmpRows = [];
                for a=1:noGridAlignmentsPerPair(img2+1, img3+1)
                    grid = gridIdsInFirstImage((img2)*5+a, img3+1);
                    found = 0;
                    for b=1:length(gridIds)
                        if(grid == gridIds(b))
                            found = b;
                        break;
                        end
                    end
                    if(found > 0)
                        tmpGridIds = [tmpGridIds; gridIdsInSecondImage((img2)*5+a, img3+1)];
                        tmpCols = [tmpCols; cols(found) + colShifts((img2)*5+a, img3+1)];
                        tmpRows = [tmpRows; rows(found) + rowShifts((img2)*5+a, img3+1)];
                        %display(['grid ' num2str(tmpGridIds(length(tmpGridIds))) ':(' num2str(tmpCols(length(tmpCols))) ',' num2str(tmpRows(length(tmpRows))) '),']);
                    end
                end
                cols = tmpCols;
                rows = tmpRows;
                gridIds = tmpGridIds;

                tmpGridIds = [];
                tmpCols = [];
                tmpRows = [];
                for a=1:noGridAlignmentsPerPair(img3+1, img1+1)
                    grid = gridIdsInFirstImage((img3)*5+a, img1+1);
                    found = 0;
                    for b=1:length(gridIds)
                        if(grid == gridIds(b))
                            found = b;
                            break;
                        end
                    end
                    if(found > 0)
                        tmpGridIds = [tmpGridIds; gridIdsInSecondImage((img3)*5+a, img1+1)];
                        tmpCols = [tmpCols; cols(found) + colShifts((img3)*5+a, img1+1)];
                        tmpRows = [tmpRows; rows(found) + rowShifts((img3)*5+a, img1+1)];
                        %display(['grid ' num2str(tmpGridIds(length(tmpGridIds))) ':(' num2str(tmpCols(length(tmpCols))) ',' num2str(tmpRows(length(tmpRows))) '),']);
                    end
                end
                cols = tmpCols;
                rows = tmpRows;
                gridIds = tmpGridIds;

                for a=1:length(gridIds)
                    if(cols(a) ~= 0 || rows(a) ~= 0)
                        noInconsistentCyclesInvolved(alignmentIds(img1+1,img2+1)) = noInconsistentCyclesInvolved(alignmentIds(img1+1,img2+1)) + 1;
                        noInconsistentCyclesInvolved(alignmentIds(img2+1,img3+1)) = noInconsistentCyclesInvolved(alignmentIds(img2+1,img3+1)) + 1;
                        noInconsistentCyclesInvolved(alignmentIds(img3+1,img1+1)) = noInconsistentCyclesInvolved(alignmentIds(img3+1,img1+1)) + 1;
                        break;
                    end
                end
             end
       end
        

        iter = iter+1;

    end

    fid = fopen(newAlignmentsFile, 'w');
    fprintf(fid, '%d\n', noAlignments);
    for i=1:noAlignments
        img1 = alignmentImages(i, 1);
        img2 = alignmentImages(i, 2);
        fprintf(fid, '%d %d %f %f\n', img1-1, img2-1, x(i), weights(i));
        fprintf(fid, '%d\n', noGridAlignmentsPerPair(img1, img2));
       
        for a=1:noGridAlignmentsPerPair(img1, img2)
            if(floor(gridIdsInFirstImage((img1-1)*5+a, img2)/maxGridAlignmentsPerPair) ~= floor(gridIdsInSecondImage((img1-1)*5+a, img2)/maxGridAlignmentsPerPair))
                     display('Error: The grids in alignment do not come from the same template!');
            end
            g1 = gridIdsInFirstImage((img1-1)*5+a, img2) - floor(gridIdsInFirstImage((img1-1)*5+a, img2)/maxGridAlignmentsPerPair)*maxGridAlignmentsPerPair;
            g2 = gridIdsInSecondImage((img1-1)*5+a, img2) - floor(gridIdsInSecondImage((img1-1)*5+a, img2)/maxGridAlignmentsPerPair)*maxGridAlignmentsPerPair;
            fprintf(fid, '%d %d %d %d %d %d %d\n', floor(gridIdsInFirstImage((img1-1)*5+a, img2)/maxGridAlignmentsPerPair), 0, 0, g1, g2, colShifts((img1-1)*5+a, img2), rowShifts((img1-1)*5+a, img2));
        end
        
        fprintf(fid, '%f %f %f\n', rotationMatrices((img1-1)*9+1, img2), rotationMatrices((img1-1)*9+2, img2), rotationMatrices((img1-1)*9+3, img2));
        fprintf(fid, '%f %f %f\n', rotationMatrices((img1-1)*9+4, img2), rotationMatrices((img1-1)*9+5, img2), rotationMatrices((img1-1)*9+6, img2));
        fprintf(fid, '%f %f %f\n', rotationMatrices((img1-1)*9+7, img2), rotationMatrices((img1-1)*9+8, img2), rotationMatrices((img1-1)*9+9, img2));
    end
    
    fclose(fid);
    
    resultWeights = [];
    for i=1:noAlignments
        resultWeights = [resultWeights; x(i)/weights(i)];
    end
    
    [w, indices] = sort(resultWeights);
    for i=1:noAlignments
        display(['alignment from ' num2str(alignmentImages(indices(i), 1)-1) ' to ' num2str(alignmentImages(indices(i), 2)-1) ' w:' num2str(w(i)) ' var:' num2str(x(indices(i))) ]);
        img1 = alignmentImages(indices(i), 1);
        img2 = alignmentImages(indices(i), 2);
       
        for a=1:noGridAlignmentsPerPair(img1, img2)
            display(['from grid ' num2str(gridIdsInFirstImage((img1-1)*5+a, img2)) ' to grid ' num2str(gridIdsInSecondImage((img1-1)*5+a,img2)) ':(' num2str(colShifts((img1-1)*5+a, img2)) ',' num2str(rowShifts((img1-1)*5+a, img2)) ')' ]);
        end
        
        if(noGridAlignmentsPerPair(img1, img2) == 0)
            display([num2str(rotationMatrices((img1-1)*9+1, img2)) ' ' num2str(rotationMatrices((img1-1)*9+2, img2)) ' ' num2str(rotationMatrices((img1-1)*9+3, img2))]);
            display([num2str(rotationMatrices((img1-1)*9+4, img2)) ' ' num2str(rotationMatrices((img1-1)*9+5, img2)) ' ' num2str(rotationMatrices((img1-1)*9+6, img2))]);
            display([num2str(rotationMatrices((img1-1)*9+7, img2)) ' ' num2str(rotationMatrices((img1-1)*9+8, img2)) ' ' num2str(rotationMatrices((img1-1)*9+9, img2))]);
        end
    end

end