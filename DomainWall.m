for lsize = 128:-4:16
    DomainWallFunc('squarehpbc',lsize,1e3,0,0);
end
function DomainWallFunc(type,lsize,nsamples,seed,para1)
dirname = sprintf('data/%s_DomainWall/%d/',type,lsize);
if ~isdir(dirname)
    disp('Dir Not Exist!')
    mkdir(dirname);
end
addpath(dirname);
generator = GGenerator(string(type));
solver = SSolver(string(type));
randstrs = {};
for i = 1:nsamples
    randstrs{end+1} = RandStr(30);
end
parfor sample = 1:nsamples
    sample
    filename = [dirname randstrs{sample} '.mat'];
    [JJJ i1 i2 nodes c1 c2] = generator(lsize,lsize+2,seed);
    JJJnew = JJJ;
    stateold = solver(nodes);
    nodesnew = nodes;
    wsize = c2-c1;
    for i = 1:length(nodesnew)
        nodei = nodesnew(i);
        for jj = 1:length(nodei.neib)
            j = nodei.neib(jj);
            ileft = min(nodei.left(jj),nodei.right(jj));
            iright = max(nodei.left(jj),nodei.right(jj));
            if mod(ileft,wsize)==wsize/2 && iright==ileft+1
                nodesnew(i).weight(jj) = -nodesnew(i).weight(jj);
                if j>i
                    JJJnew(ileft,iright) = -JJJnew(ileft,iright);
                    JJJnew(iright,ileft) = -JJJnew(iright,ileft);
                end
            end
        end
    end
    statenew = solver(nodesnew);
    mathandle = matfile(filename);
    mathandle.Properties.Writable = true;
    mathandle.i1 = i1;
    mathandle.i2 = i2;
    mathandle.c1 = c1;
    mathandle.c2 = c2;
    mathandle.JJJ = sparse(JJJ);
    mathandle.stateold = stateold;
    mathandle.statenew = statenew;
    mathandle.JJJnew = sparse(JJJnew);
    conn = abs(sign(JJJ));
    lpnew = conn.*(statenew*statenew');
    lpold = conn.*(stateold*stateold');
    diff = (sum(sum(abs(lpnew-lpold))))/4;
    mathandle.diff = full(diff);
    statenew = stateold(1)/statenew(1)*statenew;
    mathandle.vol = sum(stateold~=statenew);
    egynew = -0.5*statenew'*JJJ*statenew;
    egyold = -0.5*stateold'*JJJ*stateold;
    mathandle.egynew = egynew;
    mathandle.egyold = egyold;
    dpegy = egynew-egyold;
    mathandle.dpegy = dpegy;
end
end

% function lp = LinkProduct(JJJ,state)
% if size(state,2)~=1
%     state = state';
% end
% statemat = repmat(state,1,length(state));
% lp = statemat.*JJJ.*statemat';
% end


