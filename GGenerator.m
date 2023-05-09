% This function will return a model as a sparse matrix.
% Additionally, if it returns a 2D model, then it also returns its dual
% representation.
% Subfunctions are also returned as handles for debugging.
function generator = GGenerator(type)
type = string(type);
if type == 'square'
    generator = @Square;
elseif type=='honeycomb'
    generator = @Honeycomb;
elseif type=='triangle'
    generator = @Triangle;
elseif type=='triangle2'
    generator = @Triangle2;
elseif type=='cubic'
    generator = @Cubic;
elseif type=='bcc'
    generator = @BCC;
elseif type=='rrg';
    generator = @RRG;
elseif type=='rpg'
    generator = @RPG;
elseif type=='squarehpbc'
    generator = @SquareHPBC;
elseif type=='squarepbc'
    generator = @SquarePBC;
elseif type=="rg"
    generator = @RG;
else
    error 'There is no such type!'
end
end


function [JJJ index1 index2 nodes center1 center2] = Square(hsize, wsize, seed)
nsites = (hsize+1)*(wsize+1);
if seed~=0
    rng(seed);
else
    rng('shuffle');
end
lefts = zeros(4*nsites,1);
rights = zeros(4*nsites,1);
values = zeros(4*nsites,1);
scount = 0;
JJJ = spalloc(nsites,nsites,4*nsites);
for i = 1:hsize+1
    for j = 1:wsize+1
        index = (i-1)*(wsize+1)+j;
        neibs = zeros(4,1);
        count = 1;
        if j > 1
            neibs(count) = index-1;
            count = count+1;
        end
        if j<wsize+1
            neibs(count) = index+1;
            count = count+1;
        end
        if i>1
            neibs(count) = index-wsize-1;
            count = count+1;
        end
        if i<hsize+1
            neibs(count) = index+wsize+1;
        end
        neibs = neibs(neibs>index);
        for ineib = 1:length(neibs)
            neib = neibs(ineib);
            scount = scount+1;
            lefts(scount) = index;
            rights(scount) = neib;
            values(scount) = randn();
        end
    end
end
lefts = lefts(lefts>0);
rights = rights(lefts>0);
values = values(lefts>0);
JJJ = sparse(lefts,rights,values,nsites,nsites);
JJJ = JJJ+JJJ';
JJJ = round(1e6*JJJ+sign(JJJ));
index1 = nsites/2;
index2 = index1+1;
nodes = -1;
center1 = -1;
center2 = -1;

end




function [JJJ index1 index2 nodes center1 center2] = SquareBackUp(hsize, wsize, seed)
if mod(hsize,2)~=0 || mod(wsize,2)==0
    error('Input sizes must both be even in height and odd in width!');
end
if seed~=0
    rng(seed);
end
ncenters = hsize*wsize+1;   
nodes = [];
for i = 1:ncenters
    nodes(i).neib = [];
    nodes(i).weight = [];
    nodes(i).left = [];
    nodes(i).right = [];
end
nsites = (hsize+1)*(wsize+1);
adjmatrix = spalloc(nsites,nsites,nsites*4);
for i = 1:hsize
    for j = 1:wsize
        index = (i-1)*wsize+j;
        up = index - wsize;
        if up<1
            up = ncenters;
        end
        down = index + wsize;
        if down>hsize*wsize
            down = ncenters;
        end
        left = 0;
        right = 0;
        if mod(index,wsize)==1
            left = ncenters;
        else
            left = index-1;
        end
        if mod(index,wsize)==0
            right = ncenters;
        else
            right = index+1;
        end
        anglecount = 0;
        for neibindex = [left right up down]
            anglecount = anglecount+1;
            if neibindex>index
                rn = randn();
                randomweight = round(rn*1e6+sign(rn));
                nodes(index).neib(end+1) = neibindex;
                nodes(index).weight(end+1) = randomweight;
                nodes(neibindex).neib(end+1) = index;
                nodes(neibindex).weight(end+1) = randomweight;
                NE = (i-1)*(wsize+1)+j+1;
                NW = (i-1)*(wsize+1)+j;
                SE = i*(wsize+1)+j+1;
                SW = i*(wsize+1)+j;
                if anglecount == 1
                    adjmatrix(NW,SW) = randomweight;
                    adjmatrix(SW,NW) = randomweight;
                    nodes(index).left(end+1) = NW;
                    nodes(index).right(end+1) = SW;
                    nodes(neibindex).left(end+1) = NW;
                    nodes(neibindex).right(end+1) = SW;
                elseif anglecount == 2
                    adjmatrix(NE,SE) = randomweight;
                    adjmatrix(SE,NE) = randomweight;
                    nodes(index).left(end+1) = NE;
                    nodes(index).right(end+1) = SE;
                    nodes(neibindex).left(end+1) = NE;
                    nodes(neibindex).right(end+1) = SE;
                elseif anglecount == 3
                    adjmatrix(NW,NE) = randomweight;
                    adjmatrix(NE,NW) = randomweight;
                    nodes(index).left(end+1) = NW;
                    nodes(index).right(end+1) = NE;
                    nodes(neibindex).left(end+1) = NW;
                    nodes(neibindex).right(end+1) = NE;
                elseif anglecount == 4
                    adjmatrix(SW,SE) = randomweight;
                    adjmatrix(SE,SW) = randomweight;
                    nodes(index).left(end+1) = SW;
                    nodes(index).right(end+1) = SE;
                    nodes(neibindex).left(end+1) = SW;
                    nodes(neibindex).right(end+1) = SE;
                end
            end
        end
    end
end
JJJ = adjmatrix;
index1 = (hsize/2)*(wsize+1)+(wsize-1)/2+1;
index2 = index1+1;
center1 = (hsize/2-1)*wsize+(wsize+1)/2;
center2 = center1+wsize;
end

function [JJJ index1 index2 nodes center1 center2 cords0 bd_indices] = Honeycomb(hsize,wsize,seed)
if mod(hsize,2)~=0 || mod(wsize,2)==0
    error('Input sizes must both be even in height and odd in width!');
end
if seed~=0
    rng(seed);
end
cords = zeros(6*hsize*wsize,1);
iii = sqrt(-1);
count = 0;
for i = 1:hsize
    for j = 1:wsize
        xx = 5+(j-1)*3;
        yy = 5+(i-1)*2-mod(j+1,2);
        centercord = xx + yy*iii;
        count = count+1;
        sixnums = [centercord+1-iii centercord+2 centercord+1+iii centercord-1+iii ...,
            centercord-2 centercord-1-iii];
        cords(count:count+5) = sixnums;
        count = count+5;
    end
end
cords = unique(cords);
nsites = length(cords);
modelold = spalloc(nsites,nsites,12*nsites);
ncenters = hsize*wsize+1; 
nodes = [];
for i = 1:ncenters
    nodes(i).neib = [];
    nodes(i).weight = [];
    nodes(i).left = [];
    nodes(i).right = [];
end

for i = 1:hsize
    for j = 1:wsize
        index = (i-1)*wsize+j;
        up = index - wsize;
        if up<1
            up = ncenters;
        end
        down = index + wsize;
        if down>hsize*wsize
            down = ncenters;
        end
        ni = 0;
        nj = 0;
        leftup = 0;
        if mod(j,2)==1
            ni = i;
        else
            ni = i-1;
        end
        nj = j-1;
        if ni<1 || nj<1
            leftup = ncenters;
        else
            leftup = (ni-1)*wsize+nj;
        end
        ni = 0;
        nj = 0;
        rightup = 0;
        if mod(j,2)==1
            ni = i;
        else
            ni = i-1;
        end
        nj = j+1;
        if ni<1 || nj>wsize
            rightup = ncenters;
        else
            rightup = (ni-1)*wsize+nj;
        end
        ni = 0;
        nj = 0;
        leftdown = 0;
        if mod(j,2)==1
            ni = i+1;
        else
            ni = i;
        end
        nj = j-1;
        if ni>hsize || nj<1
            leftdown = ncenters;
        else
            leftdown = (ni-1)*wsize+nj;
        end
        ni = 0;
        nj = 0;
        rightdown = 0;
        if mod(j,2)==1
            ni = i+1;
        else
            ni = i;
        end
        nj = j+1;
        if ni>hsize || nj>wsize
            rightdown = ncenters;
        else
            rightdown = (ni-1)*wsize+nj;
        end
        
        
        neibcount = 0;
        if index == 2
            a = 3;
        end
        xx = 5+(j-1)*3;
        yy = 5+(i-1)*2-mod(j+1,2);
        centercord = xx + yy*iii;
        sixnums = [centercord+1-iii centercord+2 centercord+1+iii centercord-1+iii ...,
            centercord-2 centercord-1-iii];
        for neibindex = [rightup rightdown down leftdown leftup up]
            neibcount = neibcount+1;
            if neibindex>index
                rn = randn();
                randomweight = round(rn*1e6+sign(rn));
                nodes(index).neib(end+1) = neibindex;
                nodes(index).weight(end+1) = randomweight;
                nodes(neibindex).neib(end+1) = index;
                nodes(neibindex).weight(end+1) = randomweight;
                leftsite = find(cords==sixnums(neibcount));
                rightsite = find(cords==sixnums(mod(neibcount,6)+1));
                modelold(leftsite,rightsite) = randomweight;
                nodes(index).left(end+1) = leftsite;
                nodes(index).right(end+1) = rightsite;
                nodes(neibindex).left(end+1) = leftsite;
                nodes(neibindex).right(end+1) = rightsite;
            end
        end
    end
end
JJJ = modelold+modelold';
      

centerloc = 5+floor(wsize/2)*3+5*iii+(hsize/2-1)*2*iii;
center1 = centerloc+1+iii;
center2 = centerloc-1+iii;
index1 = find(cords==center1);
index2 = find(cords==center2);
cords0 = zeros(length(cords),2);
bd_indices = [];
for i = 1:length(cords)
    xx = real(cords(i));
    yy = imag(cords(i));
    text(xx,yy,num2str(i));
    cords0(i,:) = [xx yy];
end
xmax = min(maxk(unique(cords0(:,1)),2));
xmin = max(mink(unique(cords0(:,1)),2));
ymax = min(maxk(unique(cords0(:,2)),2));
ymin =  max(mink(unique(cords0(:,2)),2));
for i = 1:length(cords0)
    if cords0(i,1)>=xmax || cords0(i,1)<=xmin || cords0(i,2)>=ymax || cords0(i,2)<=ymin
        bd_indices(end+1) = i;
    end
end
if_plot = 1;
if if_plot
    plot(cords,'.')
    hold on
    set(gca,'YDir','reverse')
    [cols rows] = find(modelold~=0);
    locs = zeros(length(cols),2);
    for i = 1:length(cols)
        left = cols(i);
        right = rows(i);
        plot(cords([left right]),'blue')
    end
end
center1 = (hsize/2-1)*wsize+(wsize+1)/2;
center2 = center1+wsize;
end

function [JJJ index1 index2 nodes center1 center2 cords bd_indices] = Triangle(hsize,wsize,seed)

if mod(hsize,4)~=0 || mod(wsize,2)==0
    error('Input sizes must both be like 4m (m is a integer) AND an odd number!');
end
if seed~=0
    rng(seed);
end
ncenters = hsize*wsize+1;   
layersites = floor(wsize/2)+2;
nsites = (hsize/2)*(layersites*2-1)+layersites;
modelold = spalloc(nsites,nsites,nsites*6);
nodes = [];
for i = 1:ncenters
    nodes(i).neib = [];
    nodes(i).weight = [];
    nodes(i).left = [];
    nodes(i).right = [];
end

for i = 1:hsize
    for j = 1:wsize
        lsite = 0; rsite = 0; udsite = 0;
        if mod(i,2)==1
            if mod(j,2)==1
                lsite = floor(i/2)*(layersites*2-1)+floor(j/2)+1;
                rsite = lsite+1;
                udsite = lsite+layersites;
            else
                udsite = floor(i/2)*(layersites*2-1)+floor(j/2)+1;
                lsite = udsite+layersites-1;
                rsite = lsite+1;
            end
        else
            if mod(j,2)==1
                udsite = (i/2-1)*(layersites*2-1)+floor(j/2)+1+layersites;
                lsite = udsite+layersites-1;
                rsite = lsite+1;
            else
                lsite = (i/2-1)*(layersites*2-1)+floor(j/2)+1-1+layersites;
                rsite = lsite+1;
                udsite = lsite+layersites;
            end
        end
        %[i j lsite rsite udsite]
        index = (i-1)*wsize+j;
        left = index-1;
        if j==1
            left = ncenters;
        end
        right = index+1;
        if j==wsize
            right = ncenters;
        end
        updown = 0;
        if mod(i+j,2) == 0
            updown = index-wsize;
        else
            updown = index+wsize;
        end
        if updown<=0
            updown = ncenters;
        end
        if updown>ncenters
            updown = ncenters;
        end
        %[left right updown]
        neibcount = 0;
        for neibindex = [left right updown]
            neibcount = neibcount+1;
            if neibindex>index
                rn = randn();
                randomweight = round(rn*1e6+sign(rn));
                if neibcount == 1
                    modelold(lsite,udsite) = randomweight;
                    nodes(index).left(end+1) = lsite;
                    nodes(index).right(end+1) = udsite;
                    nodes(neibindex).left(end+1) = lsite;
                    nodes(neibindex).right(end+1) = udsite;
                elseif neibcount == 2
                    modelold(rsite,udsite) = randomweight;
                    nodes(index).left(end+1) = rsite;
                    nodes(index).right(end+1) = udsite;
                    nodes(neibindex).left(end+1) = rsite;
                    nodes(neibindex).right(end+1) = udsite;
                elseif neibcount == 3
                    modelold(lsite,rsite) = randomweight;
                    nodes(index).left(end+1) = lsite;
                    nodes(index).right(end+1) = rsite;
                    nodes(neibindex).left(end+1) = lsite;
                    nodes(neibindex).right(end+1) = rsite;
                end
                nodes(index).neib(end+1) = neibindex;
                nodes(index).weight(end+1) = randomweight;
                nodes(neibindex).neib(end+1) = index;
                nodes(neibindex).weight(end+1) = randomweight;
            end
        end
    end
end
        
center1 = (hsize/2-1)*wsize+(wsize+1)/2;
center2 = center1+wsize;
JJJ = modelold+modelold';
index1 = nsites/2;
index2 = nsites/2+1;
xs = zeros(nsites,1);
ys = zeros(nsites,1);
bd_indices = [];
for i = 1:(hsize/2+1)
    starting = (i-1)*(2*layersites-1)+1;
    end1 = starting+layersites-1;
    end2 = i*(2*layersites-1);
    xs(starting:end1) = 0:layersites-1;
    ys(starting:end1) = sqrt(3)*i;
    if i<(hsize/2+1)
        xs(end1+1:end2) = 0.5+(0:layersites-2);
        ys(end1+1:end2) = sqrt(3)*(i+0.5); 
    end
    if i == 1
        bd_indices = [bd_indices starting:end1];
        bd_indices = [bd_indices end1+1 end2];
    elseif i==(hsize/2+1)
        bd_indices = [bd_indices starting:end1];
    else
        bd_indices = [bd_indices starting end1 end1+1 end2];
    end
end
cords = [xs ys];
if_plot = 1;
if if_plot
    figure();
    hold on

    scatter(xs,ys,'filled')
    hold on
    for k = 1:length(xs)
        %text(xs(k),ys(k),num2str(k));
    end
    [rows cols] = find(JJJ~=0);
    for i = 1:length(cols)
        site1 = cols(i);
        site2 = rows(i);
        if site1<site2
            plot([xs(site1) xs(site2)],[ys(site1) ys(site2)],'blue')
        end
        daspect([1 1 1])
    end
end
end


function [JJJ index1 index2 nodes center1 center2 cords bd_indices] = Triangle2(hsize,wsize,seed)

if mod(hsize,4)~=0 || mod(wsize,2)==0
    error('Input sizes must both be like 4m (m is a integer) AND an odd number!');
end
if seed~=0
    rng(seed);
end
wsize = hsize*2+1;
ncenters = hsize*wsize+1;   
layersites = floor(wsize/2)+2;
nsites = (hsize/2)*(layersites*2-1)+layersites;
modelold = spalloc(nsites,nsites,nsites*6);
nodes = [];
for i = 1:ncenters
    nodes(i).neib = [];
    nodes(i).weight = [];
    nodes(i).left = [];
    nodes(i).right = [];
end

for i = 1:hsize
    for j = 1:wsize
        lsite = 0; rsite = 0; udsite = 0;
        if mod(i,2)==1
            if mod(j,2)==1
                lsite = floor(i/2)*(layersites*2-1)+floor(j/2)+1;
                rsite = lsite+1;
                udsite = lsite+layersites;
            else
                udsite = floor(i/2)*(layersites*2-1)+floor(j/2)+1;
                lsite = udsite+layersites-1;
                rsite = lsite+1;
            end
        else
            if mod(j,2)==1
                udsite = (i/2-1)*(layersites*2-1)+floor(j/2)+1+layersites;
                lsite = udsite+layersites-1;
                rsite = lsite+1;
            else
                lsite = (i/2-1)*(layersites*2-1)+floor(j/2)+1-1+layersites;
                rsite = lsite+1;
                udsite = lsite+layersites;
            end
        end
        %[i j lsite rsite udsite]
        index = (i-1)*wsize+j;
        left = index-1;
        if j==1
            left = ncenters;
        end
        right = index+1;
        if j==wsize
            right = ncenters;
        end
        updown = 0;
        if mod(i+j,2) == 0
            updown = index-wsize;
        else
            updown = index+wsize;
        end
        if updown<=0
            updown = ncenters;
        end
        if updown>ncenters
            updown = ncenters;
        end
        %[left right updown]
        neibcount = 0;
        for neibindex = [left right updown]
            neibcount = neibcount+1;
            if neibindex>index
                rn = randn();
                randomweight = round(rn*1e6+sign(rn));
                if neibcount == 1
                    modelold(lsite,udsite) = randomweight;
                    nodes(index).left(end+1) = lsite;
                    nodes(index).right(end+1) = udsite;
                    nodes(neibindex).left(end+1) = lsite;
                    nodes(neibindex).right(end+1) = udsite;
                elseif neibcount == 2
                    modelold(rsite,udsite) = randomweight;
                    nodes(index).left(end+1) = rsite;
                    nodes(index).right(end+1) = udsite;
                    nodes(neibindex).left(end+1) = rsite;
                    nodes(neibindex).right(end+1) = udsite;
                elseif neibcount == 3
                    modelold(lsite,rsite) = randomweight;
                    nodes(index).left(end+1) = lsite;
                    nodes(index).right(end+1) = rsite;
                    nodes(neibindex).left(end+1) = lsite;
                    nodes(neibindex).right(end+1) = rsite;
                end
                nodes(index).neib(end+1) = neibindex;
                nodes(index).weight(end+1) = randomweight;
                nodes(neibindex).neib(end+1) = index;
                nodes(neibindex).weight(end+1) = randomweight;
            end
        end
    end
end
        
center1 = (hsize/2-1)*wsize+(wsize+1)/2;
center2 = center1+wsize;
JJJ = modelold+modelold';
index1 = nsites/2;
index2 = nsites/2+1;
xs = zeros(nsites,1);
ys = zeros(nsites,1);
bd_indices = [];
for i = 1:(hsize/2+1)
    starting = (i-1)*(2*layersites-1)+1;
    end1 = starting+layersites-1;
    end2 = i*(2*layersites-1);
    xs(starting:end1) = 0:layersites-1;
    ys(starting:end1) = sqrt(3)*i;
    if i<(hsize/2+1)
        xs(end1+1:end2) = 0.5+(0:layersites-2);
        ys(end1+1:end2) = sqrt(3)*(i+0.5); 
    end
    if i == 1
        bd_indices = [bd_indices starting:end1];
        bd_indices = [bd_indices end1+1 end2];
    elseif i==(hsize/2+1)
        bd_indices = [bd_indices starting:end1];
    else
        bd_indices = [bd_indices starting end1 end1+1 end2];
    end
end
cords = [xs ys];
if_plot = 1;
if if_plot
    figure();
    hold on

    scatter(xs,ys,'filled')
    hold on
    for k = 1:length(xs)
        text(xs(k),ys(k),num2str(k));
    end
    [rows cols] = find(JJJ~=0);
    for i = 1:length(cols)
        site1 = cols(i);
        site2 = rows(i);
        if site1<site2
            plot([xs(site1) xs(site2)],[ys(site1) ys(site2)],'blue')
        end
        daspect([1 1 1])
    end
end
end

function [JJJ index1 index2 nodes center1 center2] = Cubic(SIZE,seed)
if seed~=0
    rng(seed);
end
nodes = -1;
center1 = -1;
center2 = -1;
if mod(SIZE,2)==0
	size1 = SIZE+1;
	size2 = SIZE+1;
	size3 = SIZE;
elseif mod(SIZE,2)==1
	size1 = SIZE;
	size2 = SIZE;
	size3 = SIZE+1;
end
NNN = size1*size2*size3;
JJJ = zeros(NNN,NNN);
if seed == 0 
    rng('shuffle');
else
    rng(seed);
end
indices = [];
coords = zeros(NNN,3);
for i = 1:size1
    for j=1:size2
        for k=1:size3
            index = 1 + (i-1)*size2*size3 + (j-1)*size3 + (k-1);
            coords(index,:) =[i,j,k];
            indices(end+1) = index;
            if k<size3
                JJJ(index,index+1)=1;
            end
            if j<size2
                JJJ(index,index+size3)=1;
            end
            if i<size1
                JJJ(index,index+size3*size2)=1;
            end
        end
    end
end
for i = 1:NNN
    for j = i+1:NNN
        if JJJ(i,j)~=0
			rn = randn();
            JJJ(i,j) = round(rn*1e6+sign(rn));
            JJJ(j,i) = JJJ(i,j);
        end
    end
end
if_plot = 1;
if if_plot
scatter3(coords(:,1),coords(:,2),coords(:,3))
hold on
[rows cols] = find(JJJ~=0);
for i = 1:length(rows)
    left = rows(i);
    right = cols(i);
    x1 = coords(left,:);
    x2 = coords(right,:);
    if left<=NNN && right<=NNN
        plot3([x1(1) x2(1)],[x1(2) x2(2)],[x1(3) x2(3)],'blue')
    elseif left>NNN && right>NNN
        plot3([x1(1) x2(1)],[x1(2) x2(2)],[x1(3) x2(3)],'red')
    else
        plot3([x1(1) x2(1)],[x1(2) x2(2)],[x1(3) x2(3)],'green')
    end
end
end
index1= NNN/2;
index2 = NNN/2+1;
end

function [JJJ index1 index2 nodes center1 center2] = BCC(SIZE,seed)
nodes = -1;
center1 = -1;
center2 = -1;
if mod(SIZE,2)==0
	size1 = SIZE+1;
	size2 = SIZE+1;
	size3 = SIZE;
elseif mod(SIZE,2)==1
	size1 = SIZE;
	size2 = SIZE;
	size3 = SIZE+1;
end
NNN = size1*size2*size3;
MMM = (size1-1)*(size2-1)*(size3-1);
JJJ = zeros(NNN+MMM,NNN+MMM);
coords = zeros(NNN+MMM,3);
if seed == 0 
    rng('shuffle');
else
    rng(seed);
end
indices = [];
for i = 1:size1
    for j=1:size2
        for k=1:size3
            index = 1 + (i-1)*size2*size3 + (j-1)*size3 + (k-1);
            coords(index,:) = [i j k];
            if k<size3
                JJJ(index,index+1)=1;
            end
            if j<size2
                JJJ(index,index+size3)=1;
            end
            if i<size1
                JJJ(index,index+size3*size2)=1;
            end
        end
    end
end
for i = 1:size1-1
    for j = 1:size2-1
        for k = 1:size3-1
            center_index = NNN + 1 + (i-1)*(size2-1)*(size3-1) + (j-1)*(size3-1) + (k-1);
%             if k<size3-1
%                 JJJ(center_index,center_index+1)=1;
%             end
%             if j<size2-1
%                 JJJ(center_index,center_index+size3-1)=1;
%             end
%             if i<size1-1
%                 JJJ(center_index,center_index+(size3-1)*(size2-1))=1;
%             end
            coords(center_index,:) = [i+0.5 j+0.5 k+0.5];
            site1 = 1 + (i-1)*size2*size3 + (j-1)*size3 + (k-1);
            site2 = site1+1;
            site3 = site1+size3;
            site4 = site1+size3+1;
            site5 = site1+size2*size3;
            site6 = site5+1;
            site7 = site5+size3;
            site8 = site5+size3+1;
            sites = [site1 site2 site3 site4 site5 site6 site7 site8];
            JJJ(center_index,sites) = 1;
            JJJ(sites,center_index) = 1;
        end
    end
end
for i = 1:NNN+MMM
    for j = i+1:NNN+MMM
        if JJJ(i,j)~=0
			rn = randn();
            JJJ(i,j) = round(rn*1e6+sign(rn));
            JJJ(j,i) = JJJ(i,j);
        end
    end
end
if_plot = 1;
if if_plot
    scatter3(coords(:,1),coords(:,2),coords(:,3))
    hold on
    [rows cols] = find(JJJ~=0);
    for i = 1:length(rows)
       left = rows(i);
       right = cols(i);
       x1 = coords(left,:);
       x2 = coords(right,:);
       if left<=NNN && right<=NNN
           plot3([x1(1) x2(1)],[x1(2) x2(2)],[x1(3) x2(3)],'blue')
       elseif left>NNN && right>NNN
           plot3([x1(1) x2(1)],[x1(2) x2(2)],[x1(3) x2(3)],'red')
       else
           plot3([x1(1) x2(1)],[x1(2) x2(2)],[x1(3) x2(3)],'green')
       end
    end
end
index1 = NNN/2;
index2 = index1+1;
end

%golan pundak (2023). Random Regular generator (https://www.mathworks.com/matlabcentral/fileexchange/29786-random-regular-generator), MATLAB Central File Exchange. Retrieved May 9, 2023. 
function [JJJ index1 index2 nodes center1 center2] = RRG(vertNum, deg)
nodes = -1;
center1 = -1;
center2 = -1;
% createRegularGraph - creates a simple d-regular undirected graph
% simple = without loops or double edges
% d-reglar = each vertex is adjecent to d edges
%
% input arguments :
%   vertNum - number of vertices
%   deg - the degree of each vertex
%
% output arguments :
%   A - A sparse matrix representation of the graph
%
% algorithm :
% "The pairing model" : create n*d 'half edges'.
% repeat as long as possible: pick a pair of half edges 
%   and if it's legal (doesn't creat a loop nor a double edge)
%   add it to the graph
% reference: http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.67.7957&rep=rep1&type=pdf

if_fix = false;
if mod(deg,1)~=0
    if_fix = true;
    rng(123);
    deg = floor(deg);
end
n = vertNum;
d = deg;
matIter = 10;

%check parameters
if mod(n*d,2)==1   
    disp('createRandRegGraph input err: n*d must be even!');
    A=[];
    return;
end

%a list of open half-edges
U = repmat(1:n,1,d);

%the graphs adajency matrix
A=sparse(n,n);

edgesTested=0; 
repetition=1;

%continue until a proper graph is formed
while ~isempty(U) && repetition < matIter
    
    edgesTested = edgesTested + 1;

    %print progress
    if mod(edgesTested, 5000)==0 
        fprintf('createRandRegGraph() progress: edges=%d/%d\n', edgesTested, n*d);    
    end

    %chose at random 2 half edges
    i1 = ceil(rand*length(U));
    i2 = ceil(rand*length(U));
    v1 = U(i1);
    v2 = U(i2);

    %check that there are no loops nor parallel edges
    if (v1 == v2) || (A(v1,v2) == 1)
        
        %restart process if needed
        if (edgesTested == n*d)           
            repetition=repetition+1;            
            edgesTested = 0;
            U = repmat(1:n,1,d);
            A = sparse(n,n);
        end
    else
        %add edge to graph
        A(v1, v2)=1;
        A(v2, v1)=1;
        
        %remove used half-edges
        v = sort([i1,i2]);
        U = [U(1:v(1)-1), U(v(1)+1:v(2)-1), U(v(2)+1:end)];
    end
end

%check that A is indeed simple regular graph
msg=isRegularGraph(A);
if ~isempty(msg)    
    disp(msg);
end
JJJ = A;
rng('shuffle');
randmat = randn(vertNum,vertNum);
randmat = sign(randmat)+round(randmat*1e6);
JJJ = JJJ.*randmat;
JJJ = triu(JJJ) + triu(JJJ)';
[rows cols] = find(JJJ~=0);
rnindex = randperm(length(rows),1);
if if_fix
    rnindex = 1;
end
index1 = rows(rnindex);
index2 = cols(rnindex);
end

%-------------------------------------------------

function msg=isRegularGraph(G)
%is G a simple d-regular graph the function returns []
%otherwise it returns a message describing the problem in G

msg=[];

%check symmetry
if (norm(G-G','fro')>0)
    msg=[msg,' is not symmetric, '];
end

%check parallel edged
if (max(G(:))>1)
    msg=[msg,sprintf(' has %d parallel edges, ',length(find(G(:)>1)) )];
end

%check that d is d-regular
d_vec=sum(G);
if min(d_vec)<d_vec(1) || max(d_vec)>d_vec(1)
    msg=[msg,' not d-regular, '];
end

%check that g doesn't contain any loops
if (norm(diag(G))>0)
    msg=[msg,sprintf(' has %d self loops, ',length(find(diag(G)>0)) )];
end
end
   



function [JJJ i1 i2 nodes c1 c2] = RPG(nsites,~,seed)
if seed~=0
    rng(seed);
end
cords = rand(nsites,2)-0.5;
dt = delaunayTriangulation(cords);
PTS = dt.Points;
CL = dt.ConnectivityList;
belongs = Convert_CL_to_DUAL(CL);
JJJ = sparse(CL, CL(:, [2:end 1]), 1);
JJJ = sign(JJJ+JJJ');
randmat = randn(nsites,nsites);
randmat = triu(randmat)+triu(randmat)';
randmat = round(1e6*randmat)+sign(randmat);
JJJ = JJJ.* randmat;
g = graph(JJJ);
ncenters = size(CL,1)+1;
nodes = [];
for i = 1:ncenters
    nodes(i).neib = [];
    nodes(i).weight = [];
    nodes(i).left = [];
    nodes(i).right = [];
end
[rows,cols] = find(JJJ~=0);
for iedge = 1:length(rows)
    left = rows(iedge);
    right = cols(iedge);
    if left>right
        continue;
    end
    node1 = belongs(left,right,1);
    node2 = belongs(right,left,2);
    jij = JJJ(left,right);
    nodes(node1).neib(end+1) = node2;
    nodes(node2).neib(end+1) = node1;
    nodes(node1).weight(end+1) = jij;
    nodes(node2).weight(end+1) = jij;
    nodes(node1).left(end+1) = left;
    nodes(node2).left(end+1) = left;
    nodes(node1).right(end+1) = right;
    nodes(node2).right(end+1) = right;
end

dists = cords(:,1).^2+cords(:,2).^2;
[dists_sorts idx] = sort(dists);
i1 = idx(1);
i2 = idx(2);
c1 = belongs(i1,i2,1);
c2 = belongs(i1,i2,2);
%plot(g, 'XData', dt.Points(:, 1), 'YData', dt.Points(:, 2));
end

function [JJJ, index1, index2, nodes, center1, center2] = SquarePBC(lsize,seed)
if seed==0
    rng('shuffle');
else
    rng(seed);
end
nodes = -1;
center1 = -1;
center2 = -1;
nsites = lsize^2;
JJJ = zeros(nsites,nsites);
index1 = nsites/2;
index2 = index1+1;
for i = 1:lsize
    for j =1:lsize
        index = (i-1)*lsize+j;
        right = index+1;
        if j==lsize
            right = right-lsize;
        end
        up = index+lsize;
        if i==lsize
            up = up-lsize^2;
        end
        rn = randn();
        rn = sign(rn)+round(rn*1e6);
        JJJ(index,right) = rn;
        rn = randn();
        rn = sign(rn)+round(rn*1e6);
        JJJ(index,up) = rn;
    end
end
JJJ = JJJ+JJJ';
end

function [JJJ index1 index2 nodes center1 center2] = SquareHPBC(hsize, wsize, seed)
% if mod(hsize,2)~=0 || mod(wsize,2)==0
%     error('Input sizes must both be even in height and odd in width!');
% end
wsize = hsize+2;
ncenters = hsize*wsize+2;   
nodes = [];
for i = 1:ncenters
    nodes(i).neib = [];
    nodes(i).weight = [];
    nodes(i).left = [];
    nodes(i).right = [];
end
nsites = (hsize+1)*(wsize);
adjmatrix = zeros(nsites,nsites);
for i = 1:hsize
    for j = 1:wsize
        index = (i-1)*wsize+j;
        up = index - wsize;
        if up<1
            up = ncenters-1;
        end
        down = index + wsize;
        if down>hsize*wsize
            down = ncenters;
        end
        left = 0;
        right = 0;
        if mod(index,wsize)==1
            left = index+wsize-1;
        else
            left = index-1;
        end
        if mod(index,wsize)==0
            right = index-wsize+1;
        else
            right = index+1;
        end
        anglecount = 0;
        for neibindex = [left right up down]
            anglecount = anglecount+1;
            if neibindex>index
                rn = randn();
                randomweight = round(rn*1e6+sign(rn));
                nodes(index).neib(end+1) = neibindex;
                nodes(index).weight(end+1) = randomweight;
                nodes(neibindex).neib(end+1) = index;
                nodes(neibindex).weight(end+1) = randomweight;
                NE = (i-1)*(wsize)+j+1;
                NW = (i-1)*(wsize)+j;
                SE = i*(wsize)+j+1;
                SW = i*(wsize)+j;
                if j == wsize
                    NE = (i-1)*wsize+1;
                    SE = i*wsize+1;
                end
                if anglecount == 1
                    i1 = NW;
                    i2 = SW;
                elseif anglecount == 2
                    i1 = NE;
                    i2 = SE;
                elseif anglecount == 3
                    i1 = NW;
                    i2 = NE;
                elseif anglecount == 4
                    i1 = SW;
                    i2 = SE;
                end
                adjmatrix(i1,i2) = randomweight;
                adjmatrix(i2,i1) = randomweight;
                nodes(index).left(end+1) = i1;
                nodes(index).right(end+1) = i2;
                nodes(neibindex).left(end+1) = i1;
                nodes(neibindex).right(end+1) = i2;
            end
        end
    end
end
JJJ = adjmatrix;
index1 = floor(nsites/2);
index2 = index1+1;
center1 = (ncenters-2)/2-wsize/2;
center2 = center1+wsize;
end

function [JJJ index1 index2 nodes c1 c2] = RG(lsize,deg)
nodes = -1; c1 = -1; c2 = -1;
rng('shuffle');
JJJ = zeros(lsize,lsize);
idxs = zeros((lsize)*(lsize-1)/2,2);
count = 0;
for i = 1:lsize-1
    for j = i+1:lsize
        count = count+1;
        idxs(count,:) = [i j];
    end
end
ridxs = randperm(lsize*(lsize-1)/2,lsize*deg/2);
for i = 1:length(ridxs)
    ridx = ridxs(i);
    left = idxs(ridx,1);
    right = idxs(ridx,2);
    if i == 1
        index1 = left;
        index2 = right;
    end
    randomweight = randn();
    randomweight = round(randomweight*1e6)+sign(randomweight);
    JJJ(left,right) = randomweight;
end
JJJ = sparse(JJJ+JJJ');
JJJ = randomize(JJJ);
end
