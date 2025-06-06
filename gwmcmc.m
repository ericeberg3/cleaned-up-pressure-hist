function [models,logP,dhat,RunTime,acceptrate]=gwmcmc(minit,logPfuns,mccount,varargin)
%% Cascaded affine invariant ensemble MCMC sampler. "The MCMC hammer"
%
% GWMCMC is an implementation of the Goodman and Weare 2010 Affine
% invariant ensemble Markov Chain Monte Carlo (MCMC) sampler. MCMC sampling
% enables bayesian inference. The problem with many traditional MCMC samplers
% is that they can have slow convergence for badly scaled problems, and that
% it is difficult to optimize the random walk for high-dimensional problems.
% This is where the GW-algorithm really excels as it is affine invariant. It
% can achieve much better convergence on badly scaled problems. It is much
% simpler to get to work straight out of the box, and for that reason it
% truly deserves to be called the MCMC hammer.
%
% (This code uses a cascaded variant of the Goodman and Weare algorithm).
%
% USAGE:
%  [models,logP]=gwmcmc(minit,logPfuns,mccount, Parameter,Value,Parameter,Value);
%
% INPUTS:
%     minit: an MxW matrix of initial values for each of the walkers in the
%            ensemble. (M:number of model params. W: number of walkers). W
%            should be atleast 2xM. (see e.g. mvnrnd).
%  logPfuns: a cell of function handles returning the log probality of a
%            proposed set of model parameters. Typically this cell will
%            contain two function handles: one to the logprior and another
%            to the loglikelihood. E.g. {@(m)logprior(m) @(m)loglike(m)}
%   mccount: What is the desired total number of monte carlo proposals.
%            This is the total number, -NOT the number per chain.
%
% Named Parameter-Value pairs:
%   'StepSize': unit-less stepsize (default=2.5).
%   'ThinChain': Thin all the chains by only storing every N'th step (default=10)
%   'ProgressBar': Show a text progress bar (default=true)
%   'Parallel': Run in ensemble of walkers in parallel. (default=false)
%   'BurnIn': fraction of the chain that should be removed. (default=0)
%
%    added YQW Nov 22, 2019
%   'FileName': save file of models, logP (default='')
%   'OutputData': whether to output data into dhat (default=0)
%   'NCores': number of cores to use if running in parallel (default=2)
%
% OUTPUTS:
%    models: A MxWxT matrix with the thinned markov chains (with T samples
%            per walker). T=~mccount/p.ThinChain/W.
%      logP: A PxWxT matrix of log probabilities for each model in the
%            models. here P is the number of functions in logPfuns.
%      dhat: (Added YQW, Nov 22, 2019)
%            A NdxWxT matrix of the predicted data (for convenience it is
%            just the 2nd of output of the two probability functions.
%            priors should output nan, likelihood will output dhat)
%
%
% Note on cascaded evaluation of log probabilities:
% The logPfuns-argument can be specifed as a cell-array to allow a cascaded
% evaluation of the probabilities. The computationally cheapest function should be
% placed first in the cell (this will typically the prior). This allows the
% routine to avoid calculating the likelihood, if the proposed model can be
% rejected based on the prior alone.
% logPfuns={logprior loglike} is faster but equivalent to
% logPfuns={@(m)logprior(m)+loglike(m)}
%
% TIP: if you aim to analyze the entire set of ensemble members as a single
% sample from the distribution then you may collapse output models-matrix
% thus: models=models(:,:); This will reshape the MxWxT matrix into a
% Mx(W*T)-matrix while preserving the order.
%
%
% EXAMPLE: Here we sample a multivariate normal distribution.
%
% %define problem:
% mu = [5;-3;6];
% C = [.5 -.4 0;-.4 .5 0; 0 0 1];
% iC=pinv(C);
% logPfuns={@(m)-0.5*sum((m-mu)'*iC*(m-mu))}
%
% %make a set of starting points for the entire ensemble of walkers
% minit=randn(length(mu),length(mu)*2);
%
% %Apply the MCMC hammer
% [models,logP]=gwmcmc(minit,logPfuns,100000);
% models(:,:,1:floor(size(models,3)*.2))=[]; %remove 20% as burn-in
% models=models(:,:)'; %reshape matrix to collapse the ensemble member dimension
% scatter(models(:,1),models(:,2))
% prctile(models,[5 50 95])
%
%
% References:
% Goodman & Weare (2010), Ensemble Samplers With Affine Invariance, Comm. App. Math. Comp. Sci., Vol. 5, No. 1, 65�80
% Foreman-Mackey, Hogg, Lang, Goodman (2013), emcee: The MCMC Hammer, arXiv:1202.3665
%
% WebPage: https://github.com/grinsted/gwmcmc
%
% -Aslak Grinsted 2015


persistent isoctave;
if isempty(isoctave)
    isoctave = (exist ('OCTAVE_VERSION', 'builtin') > 0);
end

if nargin<3
    error('GWMCMC:toofewinputs','GWMCMC requires atleast 3 inputs.')
end
M=size(minit,1);
if size(minit,2)==1
    minit=bsxfun(@plus,minit,randn(M,M*5));
end


p=inputParser;
if isoctave
    p=p.addParamValue('StepSize',2,@isnumeric); %addParamValue is chosen for compatibility with octave. Still Untested.
    p=p.addParamValue('ThinChain',10,@isnumeric);
    p=p.addParamValue('ProgressBar',true,@islogical);
    p=p.addParamValue('Parallel',false,@islogical);
    p=p.addParamValue('BurnIn',0,@(x)(x>=0)&&(x<1));
    
    % adding filename, OutputData, Ncores. YQW Nov 22, 2019
    p=p.addParamValue('FileName','',@ischar);
    p=p.addParamValue('OutputData',0,@isnumeric);
    p=p.addParamValue('Ncores',2,@isnumeric);
    p=p.parse(varargin{:});
else
    p.addParameter('StepSize',2,@isnumeric); %addParamValue is chose for compatibility with octave. Still Untested.
    p.addParameter('ThinChain',10,@isnumeric);
    p.addParameter('ProgressBar',true,@islogical);
    p.addParameter('Parallel',false,@islogical);
    p.addParameter('BurnIn',0,@(x)(x>=0)&&(x<1));
    
    % adding filename, OutputData, Ncores. YQW Nov 22, 2019
    p.addParameter('FileName','',@ischar);
    p.addParameter('OutputData',0,@isnumeric);
    p.addParameter('Ncores',2,@isnumeric);
    p.parse(varargin{:});
end
p=p.Results;

p.save = 0;
if ~isempty(p.FileName), p.save = 1; end % adding filename, YQW Nov 22, 2019

Nwalkers=size(minit,2);

if size(minit,1)*2>size(minit,2)
    warning('GWMCMC:minitdimensions','Check minit dimensions.\nIt is recommended that there be atleast twice as many walkers in the ensemble as there are model dimension.')
end

if p.ProgressBar
    progress=@textprogress;
else
    progress=@noaction;
end

if (p.Parallel)
    ParpoolObj = parpool(min([p.Ncores,mccount]));
    p.Ncores = ParpoolObj.NumWorkers;
else
    p.Ncores = 0;
end



Nkeep=ceil(mccount/p.ThinChain/Nwalkers); %number of samples drawn from each walker
mccount=(Nkeep-1)*p.ThinChain+1;

models=nan(M,Nwalkers,Nkeep); %pre-allocate output matrix
RunTime = nan(Nwalkers,Nkeep);

models(:,:,1)=minit;

if ~iscell(logPfuns)
    logPfuns={logPfuns};
end

NPfun=numel(logPfuns);
logP=nan(NPfun,Nwalkers,Nkeep);

% run one set to get sizes of 2nd output (needed to store dhat)
% added YQW Nov 22, 2019
dhatTmp = cell(NPfun,1);
tw = tic;
for fix=1:NPfun
    if p.OutputData
        [v,dhatTmp{fix}]=logPfuns{fix}(minit(:,1));
    else
        v=logPfuns{fix}(minit(:,1));
    end
    logP(fix,1,1)=v;
end
RunTime(1,1) = toc(tw);

dhat = cell(length([dhatTmp{:}]'),Nwalkers,Nkeep);
dhat(:,1,1) = num2cell([dhatTmp{:}]');
%calculate logP state initial pos of walkers
parfor (wix=2:Nwalkers, p.Ncores)
    
    tw = tic;
    logpinit = nan(NPfun,1);
    dhatinit = cell(NPfun,1);
    
    for fix=1:NPfun
        if p.OutputData
            [logpinit(fix),dhatinit{fix}]=logPfuns{fix}(minit(:,wix));
        else
            logpinit(fix)=logPfuns{fix}(minit(:,wix));
        end
        
    end
    
    RunTime(wix,1) = toc(tw);
    
    logP(:,wix,1) = logpinit;
    dhat(:,wix,1) = num2cell([dhatinit{:}]');
end

if ~all(all(isfinite(logP(:,:,1))))
    error('Starting points for all walkers must have finite logP')
end


reject=zeros(Nwalkers,1);


curm=models(:,:,1);
curlogP=logP(:,:,1);
curdhat=dhat(:,:,1);

progress(0,0,0)
totcount=Nwalkers;
for row=1:Nkeep
    for jj=1:p.ThinChain
        %generate proposals for all walkers
        %(done outside walker loop, in order to be compatible with parfor - some penalty for memory):
        %-Note it appears to give a slight performance boost for non-parallel.
        rix=mod((1:Nwalkers)+floor(rand*(Nwalkers-1)),Nwalkers)+1; %pick a random partner
        zz=((p.StepSize - 1)*rand(1,Nwalkers) + 1).^2/p.StepSize;
        proposedm=curm(:,rix) - bsxfun(@times,(curm(:,rix)-curm),zz);
        logrand=log(rand(NPfun+1,Nwalkers)); %moved outside because rand is slow inside parfor
        
        RunTimeTmp = zeros(Nwalkers,1);
        %parallel/non-parallel code is currently mirrored in
        %order to enable experimentation with separate optimization
        %techniques for each branch. Parallel is not really great yet.
        %TODO: use SPMD instead of parfor.
        %           %YQW, Nov 22, 2019: added specified number of cores
        
        parfor (wix=1:Nwalkers, p.Ncores)
            cp=curlogP(:,wix);
            lr=logrand(:,wix);
            acceptfullstep=true;
            proposedlogP=nan(NPfun,1);
            proposeddhat=cell(NPfun,1);
            if lr(1)<(numel(proposedm(:,wix))-1)*log(zz(wix))
                tw = tic;
                for fix=1:NPfun
                    if p.OutputData
                        [proposedlogP(fix), proposeddhat{fix}]=logPfuns{fix}(proposedm(:,wix)); %have tested workerobjwrapper but that is slower.
                    else
                        proposedlogP(fix)=logPfuns{fix}(proposedm(:,wix)); %have tested workerobjwrapper but that is slower.
                    end
                    if lr(fix+1)>proposedlogP(fix)-cp(fix) || ~isreal(proposedlogP(fix)) || isnan( proposedlogP(fix) )
                        %if ~(lr(fix+1)<proposedlogP(fix)-cp(fix))
                        acceptfullstep=false;
                        break
                    end
                end
                RunTimeTmp(wix) = toc(tw); 
            else
                acceptfullstep=false;
            end
            if acceptfullstep
                curm(:,wix)=proposedm(:,wix);
                curlogP(:,wix)=proposedlogP;
                curdhat(:,wix) = num2cell([proposeddhat{:}]');
            else
                reject(wix)=reject(wix)+1;
            end
        end
        
        
        totcount=totcount+Nwalkers;
        progress((row-1+jj/p.ThinChain)/Nkeep,curm,sum(reject)/totcount)
    end
    models(:,:,row)=curm;
    logP(:,:,row)=curlogP;
    dhat(:,:,row)=curdhat;
    RunTime(:,row) = RunTimeTmp;
    
    if (p.save)
        modelsSave = models(:,:,1:row);
        logpSave   = logP(:,:,1:row);
        dhatSave   = dhat(:,:,1:row);
        RunTimeSave= RunTime(:,1:row);
        save(p.FileName, 'modelsSave', 'logpSave', 'dhatSave','RunTimeSave');
        fprintf('Nkeep = %d. File saved.\n', row);
    end
    
    %progress bar
    
end
acceptrate = 1 - sum(reject)/totcount;
progress(1,0,0);
if p.BurnIn>0
    crop=ceil(Nkeep*p.BurnIn);
    models(:,:,1:crop)=[]; %TODO: never allocate space for them ?
    logP(:,:,1:crop)=[];
    dhat(:,:,1:crop)=[];
end


% TODO: make standard diagnostics to give warnings...
% TODO: make some diagnostic plots if nargout==0;



function textprogress(pct,curm,rejectpct)
persistent lastNchar lasttime starttime
if isempty(lastNchar)||pct==0
    lasttime=cputime-10;starttime=cputime;lastNchar=0;
    pct=1e-16;
end
if pct==1
    fprintf('%s',repmat(char(8),1,lastNchar));lastNchar=0;
    return
end
if (cputime-lasttime>0.1)
    
    ETA=datestr((cputime-starttime)*(1-pct)/(pct*60*60*24),13);
    progressmsg=[183-uint8((1:40)<=(pct*40)).*(183-'*') ''];
    %progressmsg=['-'-uint8((1:40)<=(pct*40)).*('-'-'�') ''];
    %progressmsg=[uint8((1:40)<=(pct*40)).*'#' ''];
    curmtxt=sprintf('% 9.3g\n',curm(1:min(end,20),1));
    %curmtxt=mat2str(curm);
    progressmsg=sprintf('\nGWMCMC %5.1f%% [%s] %s\n%3.0f%% rejected\n%s\n',pct*100,progressmsg,ETA,rejectpct*100,curmtxt);
    
    fprintf('%s%s',repmat(char(8),1,lastNchar),progressmsg);
    drawnow;lasttime=cputime;
    lastNchar=length(progressmsg);
end

function noaction(varargin)

% Acknowledgements: I became aware of the GW algorithm via a student report
% which was using emcee for python. Great stuff.
