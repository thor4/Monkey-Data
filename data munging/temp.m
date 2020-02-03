% monkey: [ "betty", "clark" ]
% day:(betty) [   '090615', '090616', '090617', '090618', '090622', '090625', 
%                 '090626', '090629', '090701', '090702', '090706', '090708', 
%                 '090709', '090901', '090903', '090916', '090917', '090921', 
%                 '090923', '090924', '090928', '090929', '090930', '091001']
%     (clark) [   '060328', '060406', '060411', '060414', '060426', '060427', 
%                 '060428', '060502', '060503', '060509', '060511', '060531', 
%                 '060601', '060602', '060824', '060825', '060831', '060907', 
%                 '061212', '061213', '061214', '061215', '061221'          ]
% good: [ 0(artifacts), 1(no artifacts) ]
% stable: [ 0(transition), 1(stable performance), 2(both) ]
% behResp: [0(incorrect), 1(correct) ]
% rule: [ 1(identity), 2(location) ]
% epoch: [ 'base', 'sample', 'delay', 'match', 'all' ]

%next up, test the sample period to ensure it's pulling correctly
monkey='betty'
% 
i=1; %init counter
for lp=days_betty %cycle through all days
%     doSomeOperation( mystruct( lp{:} ) );
    dayy = append('d',lp{:});
    [lengths.(dayy),idx(i)] = craw(path,monkey,lp{:},1,2,0,1,'base');
    mins(i)=min(lengths.(dayy)); %find the shortest length
    i=i+1;
end

idx=idx'; mins=mins'; %easier to copy-paste
fn=fieldnames(lengths);

[lengths,idx] = craw(path,monkey,'090615',1,2,1,1,'base');


%unit test for craw counter, compare to 'lengths' vector
length(1:trial_info.CueOnset(k)-1) %epoch length
for k=8:1000 %don't account for stability
    if (trial_info.good_trials(k) == 1) && ...%artifacts/none
            (trial_info.BehResp(k) == 1) && ... %correct/incorrect
            (trial_info.rule(k) == 1) %identity/location
        break
    end
end
% 
% baseline tested fine for 090615
% 
p = inputParser;
argName = 'monkey';
monkeys = { 'toejam','earl' };
validationFcn = @(x) any(validatestring(x,monkeys));
addRequired(p,argName,validationFcn);

parse(p,'toejam')


p = inputParser;
argName = 'monkey';
monkeys = [ "toejam","earl" ];
validationFcn = @(x) validateStringParameter(x,monkeys,mfilename,argName);
addRequired(p,argName,validationFcn);

function validateStringParameter(varargin)
    validatestring(varargin{:});
end


