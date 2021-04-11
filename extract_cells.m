%% --------------------------- Binarization -------------------------------
I = imread('C:\Videos\Images_make_up_length\Cycle 2\6h\2\2.tif');
if numel(size(I)) == 3
    I = rgb2gray(I);
end
figure,subplot(1,2,1);
imshow(I);
Iblur = imgaussfilt(I,1); % Blur
%I = Iblur;

%BW = zeros(size(I)); BW(I >= 205 & I <= 255) = 1;
BW = imbinarize(Iblur, 'adaptive','Sensitivity',0.35);
BW = bwmorph(BW,'open');
subplot(1,2,2);
imshow(BW);hold on;

%% ------------------------ Get Cells Properties --------------------------
BW = imclearborder(BW); % Remove cells on boarders
CC = bwconncomp(BW); % Find connected objects (cells)
stats = regionprops(CC,'Centroid', 'Area','MajorAxisLength','MinorAxisLength', 'ConvexImage',...
                     'solidity','Image'); 

for ii=1:length(stats)
   cp = regionprops(+stats(ii).ConvexImage,'Area','Perimeter');
   stats(ii).Holes = CheckHoles(stats(ii).Image);
end
showBetterCells(I, BW, stats, 'k');


%% ------------------------------------------------------------------------
%---------------------------- MEASUREMENTS --------------------------------
%--------------------------------------------------------------------------
%% ------------------------ Extract Small Cells ---------------------------
BW2idx = find([stats.Area] >15 & [stats.Area] <= 350 & [stats.MajorAxisLength] >= 1.3...
                        * [stats.MinorAxisLength] & [stats.Solidity] > 0.9);
BW2 = ismember(labelmatrix(CC),BW2idx);
% figure,subplot(1,2,1);title('Small Cells Extracted');
% imshow(BW);x
% subplot(1,2,2);
% imshow(BW2); 
showBetterCells(I,BW2,stats(BW2idx), 'g');

%% ---------------------- Measure Short Cells -----------------------------
cc = bwconncomp(BW2);
LStruct = regionprops(cc, 'MajorAxisLength');
SmallBacLength = [LStruct.MajorAxisLength] ./ 5.9034; 
fprintf('Short Cell Measurements: DONE!\n');

%% ------------------ Extract Larger Straight Cells -----------------------
BW4idx = find([stats.Area] < 3500 & [stats.Area] > 350 & [stats.MajorAxisLength] >= 2 * [stats.MinorAxisLength]...
                & [stats.Solidity] >= 0.892);% & [stats.Solidity] <= 0.9);
BW4 = ismember(labelmatrix(CC),BW4idx);

% figure,subplot(1,2,1);
% imshow(BW);
% subplot(1,2,2);
% imshow(BW4);
showBetterCells(I,BW4,stats(BW4idx), 'c');


%% ------------------- Measure Larger Straight Cells ----------------------
cc = bwconncomp(BW4);
LStruct = regionprops(cc, 'MajorAxisLength');
LagStrBacLength = [LStruct.MajorAxisLength] ./ 5.9034; 
fprintf('Larger Straight Cell Measurements: DONE!\n');

%% --------------------- Extract Long Curvy Cells -------------------------
BW3idx = find([stats.Area] > 2000 & [stats.Solidity] < 0.785);
BW3idx = GetRidOfBadCells(BW3idx, stats, 200);
BW3 = ismember(labelmatrix(CC),BW3idx);

% figure,subplot(1,2,1);
% imshow(BW);
% subplot(1,2,2);
% imshow(BW3);
showBetterCells(I,BW3,stats(BW3idx), 'y');

%% --------------------- Extract Short Curvy Cells ------------------------
BW5idx = find([stats.Area] > 300 & [stats.Area] <= 2000 & [stats.MajorAxisLength] >= 2 * [stats.MinorAxisLength] ...
                & [stats.Solidity] >= 0.49 &  [stats.Solidity] < 0.892);
BW3idx = GetRidOfBadCells(BW3idx, stats, 100);
BW5 = ismember(labelmatrix(CC),BW5idx);

% figure,subplot(1,2,1);
% imshow(BW);
% subplot(1,2,2);
% imshow(BW5);
showBetterCells(I,BW5,stats(BW5idx), 'y');
%% --------------------- Show All Extracted Cells -------------------------
BWTotalidx = [BW2idx BW3idx BW4idx BW5idx];
BWTotal = ismember(labelmatrix(CC), BWTotalidx);
showBetterCells(I, BWTotal, stats(BWTotalidx), 'k');

%% ----------------------------- Skeleton ---------------------------------
skel = bwskel(logical(BW3+BW5));
skelFilled = imfill(skel, 'hole');
skel = bwskel(skelFilled);
figure, imshow(skel);

skelcc = bwconncomp(skel); % Find connected objects (cells)
%% --------------- Measure Long and Short Curved Cells --------------------
indivCell = regionprops(skelcc, 'Image');

CurvBacContourLength = zeros(1,length(indivCell));
for ii = 1:length(indivCell)
    ima = indivCell(ii).Image;
    figure, imshow(ima);
    %figure, imshow(bwmorph(ima, 'spur'));
    skelRem = SpurRemover(ima);
    E = bwmorph(logical(skelRem), 'endpoints');
    Ecorrect = EndptsCorrect(E, logical(skelRem));
    [i,j] = find(Ecorrect);
    counter = 0;
    while length(i) > 2
        skelRem = SpurRemover(logical(skelRem));
        Eloop = bwmorph(logical(skelRem), 'endpoints');
        EloopCorrect = EndptsCorrect(Eloop, logical(skelRem));
        [i, j] = find(EloopCorrect);
        
        counter = counter + 1;
        if counter == 20
            fprintf('Error in repetitive remover at ii = %d \n', ii);
            break;
        end
    end
    figure,imshow(logical(skelRem));
    hold on;
    plot(j,i,'wo');
    if length(i) == 2
        geodisMap = bwdistgeodesic(logical(skelRem), j(1), i(1), 'quasi-euclidean');
        CurvBacContourLength(ii) = geodisMap(i(end), j(end)) / 5.9034; % Print contour length
    else
        fprintf('ERROR at ii = %d \n', ii);
        continue
    end
end
display([SmallBacLength, LagStrBacLength, CurvBacContourLength]);

%% --------------------------- Write Excel --------------------------------
AllLength = [SmallBacLength LagStrBacLength CurvBacContourLength];

destination = 'C:\Videos\Images_make_up_length\makeup_data.xlsx';
xlswrite(destination, {'C2 H6'},'sheet1', 'D1');
xlswrite(destination, AllLength', 'sheet1', 'D140');
%% Close all windows
if input('Close all windows?','s') =='y' || input('Close all windows?','s') =='Y'
    close all;
end
%% ------------------------------------------------------------------------
%
%%                          END OF PROGRAM HERE
%
%% ------------------------------------------------------------------------
%% Functions
function skelRem = SpurRemover(skel)
    B = bwmorph(skel, 'branchpoints');
    E = bwmorph(skel, 'endpoints');
    Ecorr = EndptsCorrect(E, skel);
    [Ey,Ex] = find(Ecorr);
    B_loc = B;
    
    Dmask = false(size(skel));
    for k = 1:numel(Ex)
        D = bwdistgeodesic(skel,Ex(k),Ey(k));
        distanceToBranchPt = min(D(B_loc));
        if distanceToBranchPt < 15
            Dmask(D < distanceToBranchPt) = true;
        end
    end
    skelRem = skel - Dmask;
    %figure,imshow(logical(skelRem));
    %hold all;
    %[By,Bx] = find(B);% plot(Bx,By,'ro');
    %plot(Ex,Ey,'wo');
end

function endCorr = EndptsCorrect(E, image)
    endCorr = E;
    [Ey, Ex] = find(E);
    [row,col] = size(image);
    for ii = 1:length(Ex)
        if(Ey(ii) > 1 && Ex(ii) > 1)
            neighbors(1) = image(Ey(ii)-1,Ex(ii)-1); % Upper left. 
        end
        if(Ey(ii) > 1)
            neighbors(2) = image(Ey(ii)-1,Ex(ii)); % Upper middle.   
        end
        if(Ey(ii) > 1 && Ex(ii) < col)
            neighbors(3) = image(Ey(ii)-1,Ex(ii)+1); % Upper right. 
        end
        if(Ex(ii) > 1)  
            neighbors(4) = image(Ey(ii),Ex(ii)-1); % left.  
        end
        if(Ey(ii) < row && Ex(ii) > 1)
            neighbors(8) = image(Ey(ii)+1,Ex(ii)-1); % Lower left. 
        end
        if(Ex(ii) < col)
            neighbors(5) = image(Ey(ii),Ex(ii)+1); % right. 
        end
        if(Ey(ii) < row && Ex(ii) < col)
            neighbors(6) = image(Ey(ii)+1,Ex(ii)+1); % Lower left.  
        end
        if(Ey(ii) < row)
            neighbors(7) = image(Ey(ii)+1,Ex(ii)); % lower middle. 
        end

        if length(find(neighbors)) > 2
            endCorr(Ey(ii),Ex(ii)) = false; 
        end
    end            
end

function holeProp = CheckHoles(image)
    filledImage = imfill(image, 'holes');
    holes = filledImage - image;
    holeProp = false;
    if sum(holes(:)) > 50
        holeProp = true;
    end
end

function showBetterCells(I, bw, stat, color)
    [B,~] = bwboundaries(bw,'noholes');
    figure,imshow(I); hold on;
    for ii=1:length(B)
      boundary = B{ii};
      plot(boundary(:,2), boundary(:,1),...
          color,'LineWidth',1);
      h = text(stat(ii).Centroid(1),stat(ii).Centroid(2)+1, num2str(ii));
      set(h,'Color','r','FontSize',8);
    end
end

function idx = GetRidOfBadCells(idx, stats, area) 
    for ii = idx
       holesBW = bwconvhull(stats(ii).Image) - stats(ii).Image;
       cc = bwconncomp(holesBW);
       stat = regionprops(cc, 'Area');
       if length(find([stat.Area] >= area)) > 3 && stats(ii).Area <= 5000
           idx(idx == ii) = [];
       end
    end
end

function showCells(bw)
    [B,L,~,~] = bwboundaries(bw,'noholes');

    figure,imshow(bw);
    hold on;
    for k=1:length(B)
      boundary = B{k};
      plot(boundary(:,2), boundary(:,1),...
          'r','LineWidth',2);

      %randomize text position for better visibility
      rndRow = ceil(length(boundary)/(mod(rand*k,7)+1));
      col = boundary(rndRow,2); row = boundary(rndRow,1);
      h = text(col+1, row-1, num2str(L(row,col)));
      set(h,'Color','w','FontSize',8);
    end
end