function solver = SSolver(type)
if strcmp(type,'square') || strcmp(type,'honeycomb') || strcmp(type,'square') || strcmp(type,'triangle')  || strcmp(type,'triangle2')|| strcmp(type,'rpg')|| strcmp(type,'squarehpbc')
    solver = @GS_nodes;
elseif strcmp(type,'cubic') || strcmp(type,'bcc') || strcmp(type,'rrg') || strcmp(type,'rg') || strcmp(type,'squarepbc')
    solver = @GurobiComputeGround;
end
if type=="square"
    solver = @k4Solver;
end
end

%% Getting GS via MWPM, only deal with free boundary conditions.
function state = GS_nodes(nodes)
ncenters = length(nodes);
strangecell = zeros(ncenters,1); 
for i = 1:ncenters
    iweights = nodes(i).weight;
    if mod(sum(iweights<0),2)~=0
        strangecell(i) = i;
    end
end
strangecell = strangecell(strangecell~=0);
nstrangecell = length(strangecell);
ggg = spalloc(ncenters,ncenters,20*ncenters);
linkproducts = spalloc(ncenters,ncenters,20*ncenters);
for i = 1:ncenters
    iweights = nodes(i).weight;
    ineibs = nodes(i).neib;
    nneibs = length(ineibs);
    for j = 1:nneibs
        if ggg(i,ineibs(j))==0 
            ggg(i,ineibs(j)) = iweights(j);
            linkproducts(i,ineibs(j)) = sign(iweights(j));
        else
            if abs(iweights(j)) < abs(ggg(i, ineibs(j)))
                ggg(i,ineibs(j)) = iweights(j);
                linkproducts(i,ineibs(j)) = sign(iweights(j));
            end
        end
    end
end
dualmodel = abs(ggg);
ggg = graph(dualmodel);

dists = distances(ggg,strangecell,strangecell);
noutputsites = nstrangecell;
noutputedges = (nstrangecell-1)*nstrangecell/2;
distarray = triu(dists,1);
distarray = distarray';
distarray = distarray(distarray~=0);
[matchl matchr time] = blossom(noutputsites,noutputedges,distarray);
nmatching = noutputsites/2;
unsat_bonds = zeros(3*ncenters,2);
count = 0;
for i = 1:nmatching
    source = strangecell(matchl(i)+1);
    dest = strangecell(matchr(i)+1);
    spath = shortestpath(ggg,source,dest);
    for j = 1:length(spath)-1
        count = count+1;
        unsat_bonds(count,1) = spath(j);
        unsat_bonds(count,2) = spath(j+1);
    end
end
unsat_bonds = unsat_bonds(1:count,:);

for i=1:size(unsat_bonds,1)
    linkproducts(unsat_bonds(i,1),unsat_bonds(i,2)) = ...,
        linkproducts(unsat_bonds(i,1),unsat_bonds(i,2))*(-1);
    linkproducts(unsat_bonds(i,2),unsat_bonds(i,1)) =...,
        linkproducts(unsat_bonds(i,2),unsat_bonds(i,1))*(-1);
end


nsites = max([nodes.left]);
nsites = max(nsites,max([nodes.right]));
edgestate = spalloc(nsites,nsites,12*nsites);
state = zeros(nsites,1);
for index = 1:ncenters
    ineibs = nodes(index).neib;
    ilefts = nodes(index).left;
    irights = nodes(index).right;
    iweights = nodes(index).weight;
    boundmatch = find(ineibs==ncenters);
    if length(boundmatch)>=2
        iweights = nodes(index).weight;
        [minim minid] = min(abs(iweights(boundmatch)));
        smallindex = boundmatch(minid);
    end
    for neibi = 1:length(ineibs)
        left = ilefts(neibi);
        right = irights(neibi);
        if edgestate(left,right)==0

        if ineibs(neibi) == ncenters && neibi ~= smallindex && smallindex~=0
            edgestate(left,right) = sign(iweights(neibi));
            edgestate(right,left) = sign(iweights(neibi));
        else
            edgestate(left,right) = linkproducts(index,ineibs(neibi));
            edgestate(right,left) = linkproducts(index,ineibs(neibi));
        end
        
        end
    end
end
state(1) = 1;
ggg = graph(edgestate);
T = dfsearch(ggg,1,'allevents');
indices = T.Event=='edgetonew';
findings = T(indices,:);
for i = 1:height(findings)
    site1 = findings(i,:).Edge(1);
    site2 = findings(i,:).Edge(2);
    state(site2) = state(site1)*edgestate(site1,site2);
end
end

function result = GurobiComputeGround(JJJ)
model.Q = -sparse(2*JJJ);
model.obj = (2*sum(JJJ));
model.objcon = -0.5*sum(sum(JJJ));
model.vtype = 'B';
model.A = sparse(zeros(length(JJJ),length(JJJ)));
model.modelsense = 'min';
params.OutputFlag = 0;
result = gurobi(model,params);
result = result.x;
result = 2*result-1;
end

function state = k4Solver(JJJ)
nsites = size(JJJ,1);
hsize = floor(sqrt(nsites))-1;
if nsites~=(hsize+1)*(hsize+2) && nargin<2
    error('Can not solve such a system size!');
end
wsize = hsize+1;
ncenters = (hsize+2)*(wsize+2);
k4JJJs = zeros(ncenters*4*4,3); 
k4JJJcount = 0;
connections = k4JJJs; 
connectioncount = 0;
k4site1s = k4JJJs;  
k4site2s = k4JJJs;
for i = 1:hsize+1
    for j = 0:wsize+1
        left = (i-1)*(wsize+1)+j;
        right = left+1;
        up = (i-1)*(wsize+2)+j+1+ncenters;
        down = i*(wsize+2)+j+1;
        if j==0 || j==wsize+1
            connectioncount = connectioncount+1;
            connections(connectioncount,:) = [up down 1];
            continue;
        end
        k4JJJcount = k4JJJcount+1;
        k4JJJs(k4JJJcount,:) = [up down JJJ(left,right)];
        k4site1s(k4JJJcount,:) = [up down left];
        k4site2s(k4JJJcount,:) = [up down right];
    end
end
for i = 0:hsize+1
    for j = 1:wsize+1
        up = (i-1)*(wsize+1)+j;
        down = up+(wsize+1);
        left = i*(wsize+2)+j+ncenters*3;
        right = left+1-ncenters;
        if i==0 || i==hsize+1
            connectioncount = connectioncount+1;
            connections(connectioncount,:) = [left right 1];
            continue;
        end
        k4JJJcount = k4JJJcount+1;
        k4JJJs(k4JJJcount,:) = [left right JJJ(up,down)];
        k4site1s(k4JJJcount,:) = [left right up];
        k4site2s(k4JJJcount,:) = [right left down];
    end
end
for i = 1:ncenters*4
    for neib = i+ncenters:ncenters:ncenters*4
        connectioncount = connectioncount+1;
        connections(connectioncount,:) = [i neib 1];
    end
end
k4JJJs = k4JJJs(1:k4JJJcount,:);
k4site1s = k4site1s(1:k4JJJcount,:);
k4site2s = k4site2s(1:k4JJJcount,:);
connections = connections(1:connectioncount,:);
k4JJJ = sparse(k4JJJs(:,1),k4JJJs(:,2),k4JJJs(:,3),ncenters*4,ncenters*4);
k4JJJ = k4JJJ+k4JJJ';
connection = sparse(connections(:,1),connections(:,2),connections(:,3),ncenters*4,ncenters*4);
connection = connection + connection';
k4site1 = sparse(k4site1s(:,1),k4site1s(:,2),k4site1s(:,3),ncenters*4,ncenters*4);
k4site2 = sparse(k4site2s(:,1),k4site2s(:,2),k4site2s(:,3),ncenters*4,ncenters*4);
connection = connection+sign(abs(k4JJJ));
k4site1 = k4site1+k4site1';
k4site2 = k4site2+k4site2'; 
[rows, cols] = find(connection~=0);
nedges = sum(sum(connection~=0))/2;
weights = zeros(nedges,1);
edges = zeros(nedges,2);
count = 0;
for i = 1:length(rows)
    site1 = rows(i);
    site2 = cols(i);
    if(site2>site1)
        count = count+1;
        weights(count) = k4JJJ(site1,site2);
        edges(count,:) = [site1 site2];
    end
end
nedges = full(nedges);
[matchl, matchr] = blossom2(ncenters*4,nedges,weights,edges-1);
edgestate = sign(abs(JJJ));
for i = 1:length(matchl)
    match1 = matchl(i)+1;
    match2 = matchr(i)+1;
    site1 = k4site1(match1,match2);
    site2 = k4site2(match1,match2);
    if site1~=0
        edgestate(site1,site2) = -1;
        edgestate(site2,site1) = -1;
    end
end
state = zeros(nsites,1);
state(1) = 1;
for i = 1:hsize+1
    for j = 1:wsize+1
        index = (i-1)*(wsize+1)+j;
        if index==1
            continue;
        end
        index0 = 0;
        if i>1
            index0 = index-wsize-1;
        end
        if j>1
            index0 = index-1;
        end
        state(index) = state(index0)*edgestate(index,index0);
    end
end
end