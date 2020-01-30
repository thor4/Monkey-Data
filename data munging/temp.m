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

monkey='betty'
% 
[lengths,idx] = craw(path,monkey,'090615',1,2,1,1,'base');
% 
% 
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


