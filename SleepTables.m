addpath('P:\3013077.01\Scripts_301307701\Shervin\ST\sleeptrip-master')
st_defaults
%%
files=dir('*.txt');
% for file = files
%     varname = regexp(file.name, '^.\w+', 'match');
%     varname = genvarname(varname{:});
%     data.(varname) = load(file.name);
% end
files=struct2cell(files)
load('P:\3013077.01\Scripts_301307701\Shervin\All Hyps\HypFiles.mat')
%%
for N=1:numel(files)
cfg = [];
cfg.scoringfile   = files{N}
cfg.scoringformat = 'spisop';
cfg.standard      = 'aasm'; % 'aasm' or 'rk'

scoring{N} = st_read_scoring(cfg);
end
%%
for N=1:numel(scoring)
% create a sleep table with the scoring parameters
cfg = [];
res_sleepdescriptive{N} = st_scoringdescriptives(cfg, scoring{N});

% lets take a look what we got in there
res_sleepdescriptive{N}.table
end

%%
subject.name = 'sub#'
% export the results
cfg = [];
cfg.prefix = 'example';
cfg.infix  = subject.name;
cfg.posfix = '';
filelist_res_sleepdescriptives = st_write_res(cfg, res_sleepdescriptive); 

% note: you can write mutliple results with st_write_res(cfg, res_sleep1, res_sleep2, ...)
%       or if they are in a cell array with st_write_res(cfg, res_all{:})
for N=1:numel(res_sleepdescriptive)
BigTable(N,:)=res_sleepdescriptive{N}.table;
end

writetable(BigTable, 'SleepParameters.xlsx')

    