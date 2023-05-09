% for lsize = 6
%     parent_dir = ['./data/rrg/' num2str(lsize) '/'];
%     Summarizer(parent_dir);
% end
Summarizer('./data/squarehpbc/',"squarehpbc",["96"])
%Summarizer('./data/rpg/')
% Summarizer('./data/squarehpbc_DomainWall/',"squarehpbc_domainwall")
%  Summarizer('./data/cubic/',"square",["5","6"])
%Summarizer('./data/honeycomb/',"honeycomb",["96"])
%Summarizer('./data/triangle/',"triangle")
function Summarizer(parent_dir,type,dirs)
    if parent_dir(end) ~= '/'
        parent_dir(end+1) = '/';
    end
    subdirs = dir(parent_dir);
    for idir = 1:length(subdirs)
        dirname = subdirs(idir).name
        dirchar = char(dirname);
        if dirname=="." || dirname==".." || dirchar(end)=='t'
            continue;
        end
        if nargin>2 && length(find(dirs==string(dirname)))==0
            continue;
        end
        lsize = str2num(dirname);
        if type == "honeycomb"
            generator =  GGenerator("honeycomb");
            [~,~,~,~,~,~,cords,bd_indices] = generator(lsize,lsize+1,0);
            save('./data/auxdata.mat','cords','bd_indices');
        end
        if type == "triangle"
            generator =  GGenerator("triangle");
            [~,~,~,~,~,~,cords,bd_indices] = generator(lsize,lsize+1,0);
            save('./data/auxdata.mat','cords','bd_indices');
        end
        disp("Getting Infos!")
        infos = Get_Infos([parent_dir dirname '/'],type);
        filename = [parent_dir dirname '.mat'];
        save(filename,"infos");
    end
end

function infos = Get_Infos(dirname,type)
    infos = struct;
    files = dir([dirname '*.mat']);
    nfiles = length(files);
    ifile = 1;
    filename = [dirname files(ifile).name];
    info = Get_Info(filename,type);
    quantities = fieldnames(info);
    for i = 1:length(quantities)
        eval(sprintf('infos.%ss=zeros(%d,1);',quantities{i},nfiles));
    end
    for ifile = 1:nfiles
        ifile
        filename = [dirname files(ifile).name]
        try
            info = Get_Info(filename,type);
        catch
            continue;
        end
        quantities = fieldnames(info);
        for i = 1:length(quantities)
            eval(sprintf('infos.%ss(%d,1)=%f;',quantities{i},ifile,full(info.(quantities{i}))));
        end
    end
end

function info = Get_Info(filename,type)
    info = struct;
    matload = load(filename);
    if type~= 'squarehpbc_domainwall'
    info.jc = matload.jc;
    end
    info.dpegy = matload.dpegy;
    info.vol = matload.vol;
    info.diff = matload.diff;
    [rad bdt] = Droplet_Geometry(matload.JJJ,matload.statenew,matload.stateold);
    info.rad = rad;
    info.bdt = bdt;
    [rad2 bdt2] = Droplet_Geometry2(matload.JJJ,matload.statenew,matload.stateold,type);
    info.rad2 = rad2;
    info.bdt2 = bdt2;
%     stateold = matload.stateold;
%     statenew = matload.statenew;
%     JJJ = matload.JJJ;
%     satis = sign(JJJ).*(stateold*stateold');
%     info.sat = full(sum(sum(satis>0))/2);
%     info.unsat = full(sum(sum(satis<0))/2);
%     lpdiff = (stateold*stateold'-statenew*statenew').*abs(sign(JJJ));
%     lpdiff = abs(sign(lpdiff));
%     satis_dp = sign(lpdiff).*(stateold*stateold');
%     info.sat_dp = full(sum(sum(satis_dp>0))/2);
%     info.unsat_dp = full(sum(sum(satis_dp<0))/2);
end

% This function is to remove the trivial folders, namely '.' and '..'.
function mod_folders = Remove_Trivial_Folder(folders)
    indices = 1:length(folders);
    rm_indices = [];
    for i = 1:length(folders)
        if strcmp(folders(i).name,'.') || strcmp(folders(i).name,'..')
            rm_indices(end+1) = i;
        end
    end
    indices(rm_indices) = [];
    mod_folders =folders(indices);
end

% function Summarizer(parent_dir)
%     workdir = sprintf('./data/%s/');
%     filename = sprintf('./data/%s/%s_summary.mat');
%     if ~isfile(filename)
%         delete(filename);
%     end
%     mathandle = matfile(filename);
%     mathandle.Properties.Writable = true;
%     folders = dir(workdir);
%     quantities = {'dpegy','jc','diff','vol'};
%     for i = 1:length(quantities)
%         eval([quantities{i} 'ss={};']);
%     end
%     for ifo = 1:length(folders)
%         lsize = str2num(folders(ifo).name);
%         if length(lsize)==0
%             continue;
%         end
%         dirname = [workdir num2str(lsize) '/'];
%         infos = Summ_OneSize(dirname);
%         for i = 1:length(quantities)
%             eval(sprintf('%sss{end}=infos.%ss;'), ...
%                 quantities{i},quantities{i});
%         end
%     end
%     for i = 1:length(quantities)
%         eval(sprintf('mathandle.%sss=%sss;'),quantities{i},quantities{i});
%     end
% end
