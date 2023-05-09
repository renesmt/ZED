function CritDropletFunc(type,lsize,nsamples,seed,para1)
if strcmp(type,"honeycomb") || strcmp(type,"triangle") || strcmp(type,"triangle2")|| strcmp(type,"rpg")|| strcmp(type,"squarehpbc")
    general_type = "2D";
elseif strcmp(type,"bcc") || strcmp(type,"cubic") || strcmp(type,"squarepbc")
    general_type = "3D";
elseif strcmp(type,"rrg") || strcmp(type,"rg") 
    general_type = type;
elseif strcmp(type,"square")
    general_type = "square";
end
dirname = sprintf('./data/%s/%d/',type,lsize)
if type == "rrg"
    dirname = sprintf('data/%s/%d/%d/',type,para1,lsize);
end
if type == "rg"
    dirname = sprintf('data/%s_%d/%d/',type,para1,lsize);
end
if ~isdir(dirname)
    disp('Dir Not Exist!')
    mkdir(dirname);
end
generator = GGenerator(string(type));
solver = SSolver(string(type));
general_type
parfor sample = 1:nsamples
    sample
    filename = [dirname RandStr(30) '.mat'];
    if general_type == "rrg" || general_type == "rg"
        [JJJ i1 i2 nodes c1 c2] = generator(lsize,para1);
    elseif general_type == "3D"
        [JJJ i1 i2 nodes c1 c2] = generator(lsize,seed);
    elseif general_type == "2D" || general_type == "square"
        [JJJ i1 i2 nodes c1 c2] = generator(lsize,lsize+1,seed);
    end
    if general_type=="2D"
        stateold = solver(nodes);
    elseif general_type=="square"
        stateold = solver(JJJ);
    else
        stateold = solver(full(JJJ));
    end
    lpcenter = stateold(i1)*stateold(i2);
    if lpcenter>0
        newJ = -2e8;
    else
        newJ = 2e8;
    end
    if general_type=="2D"
        changeindex = find(nodes(c1).neib==c2);
        nodes(c1).weight(changeindex) = newJ;
        changeindex = find(nodes(c2).neib==c1);
        nodes(c2).weight(changeindex) = newJ;
        statenew = solver(nodes);
    else
        JJJnew = JJJ;
        JJJnew(i1,i2) = newJ;
        JJJnew(i2,i1) = newJ;
        if general_type=="square"
            statenew = solver(JJJnew);
        else
            statenew = solver(full(JJJnew));
        end
    end
    tic
    JJJ= sparse(JJJ);
    conn = abs(sign(JJJ));
    lpnew = Gen_LP(conn,statenew);
    lpold = Gen_LP(conn,stateold);
    diff = (sum(sum(abs(lpnew-lpold))))/4;
    statenew = stateold(1)/statenew(1)*statenew;
    vol = sum(stateold~=statenew);
    egynew = -0.5*statenew'*JJJ*statenew;
    egyold = -0.5*stateold'*JJJ*stateold;
    dpegy = egynew-egyold;
    cproduct = stateold(i1)*stateold(i2);
    jc = JJJ(i1,i2)-cproduct*dpegy/2;
    parsave(filename,i1,i2,c1,c2,JJJ,stateold,statenew,newJ,diff, ...,
        vol,egynew,egyold,dpegy,jc)
    toc
end
end


function lp = Gen_LP(conn,state)
[rows cols] = find(conn~=0);
nedges = length(rows);
nsites = length(state);
lps = zeros(nedges,3);
for i = 1:nedges
    left = rows(i);
    right = cols(i);
    lps(i,:) = [left right state(left)*state(right)];
end
lp = sparse(lps(:,1),lps(:,2),lps(:,3),nsites,nsites);
end
% function lp = LinkProduct(JJJ,state)
% if size(state,2)~=1
%     state = state';
% end
% statemat = repmat(state,1,length(state));
% lp = statemat.*JJJ.*statemat';
% end

function parsave(filename,i1,i2,c1,c2,JJJ,stateold,statenew,newJ,diff,vol,egynew,egyold,dpegy,jc)
    save(filename,"i1","i2","c1","c2","JJJ","stateold","statenew","newJ","diff", ...,
        "vol","egynew","egyold","dpegy","jc")
end



