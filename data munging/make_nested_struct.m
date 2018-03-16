load('m1GoodStableRule1PingRej-split_by_BehResp_and_Area.mat')
load('m2GoodStableRule1PingRej-split_by_BehResp_and_Area.mat')

%define nested fields
monkeys = (1:2);
behResp = ["Correct","Incorrect"];
rule = ["Rule1","Rule2"];
areas = ["9L","8B","6DR","8AD","vPFC","dPFC","LIP","MIP","PE","PG","PEC"];
%add data to structure
monkey(1).Correct.Rule1.a9L = m1CorrectRule19L;
monkey(1).Correct.Rule1.a8B = m1CorrectRule18B;
monkey(1).Correct.Rule1.adPFC = m1CorrectRule1dPFC;
monkey(1).Correct.Rule1.avPFC = m1CorrectRule1vPFC;
monkey(1).Correct.Rule1.aLIP = m1CorrectRule1LIP;
monkey(1).Correct.Rule1.aMIP = m1CorrectRule1MIP;
monkey(1).Correct.Rule1.aPEC = m1CorrectRule1PEC;
monkey(1).Correct.Rule1.aPG = m1CorrectRule1PG;
monkey(1).Incorrect.Rule1.a9L = m1IncorrectRule19L;
monkey(1).Incorrect.Rule1.a8B = m1IncorrectRule18B;
monkey(1).Incorrect.Rule1.adPFC = m1IncorrectRule1dPFC;
monkey(1).Incorrect.Rule1.avPFC = m1IncorrectRule1vPFC;
monkey(1).Incorrect.Rule1.aLIP = m1IncorrectRule1LIP;
monkey(1).Incorrect.Rule1.aMIP = m1IncorrectRule1MIP;
monkey(1).Incorrect.Rule1.aPEC = m1IncorrectRule1PEC;
monkey(1).Incorrect.Rule1.aPG = m1IncorrectRule1PG;
monkey(2).Correct.Rule1.a6DR = m2CorrectRule16DR;
monkey(2).Correct.Rule1.a8B = m2CorrectRule18B;
monkey(2).Correct.Rule1.adPFC = m2CorrectRule1dPFC;
monkey(2).Correct.Rule1.a8AD = m2CorrectRule18AD;
monkey(2).Correct.Rule1.aLIP = m2CorrectRule1LIP;
monkey(2).Correct.Rule1.aPE = m2CorrectRule1PE;
monkey(2).Correct.Rule1.aPEC = m2CorrectRule1PEC;
monkey(2).Correct.Rule1.aPG = m2CorrectRule1PG;
monkey(2).Incorrect.Rule1.a6DR = m2IncorrectRule16DR;
monkey(2).Incorrect.Rule1.a8B = m2IncorrectRule18B;
monkey(2).Incorrect.Rule1.adPFC = m2IncorrectRule1dPFC;
monkey(2).Incorrect.Rule1.a8AD = m2IncorrectRule18AD;
monkey(2).Incorrect.Rule1.aLIP = m2IncorrectRule1LIP;
monkey(2).Incorrect.Rule1.aPE = m2IncorrectRule1PE;
monkey(2).Incorrect.Rule1.aPEC = m2IncorrectRule1PEC;
monkey(2).Incorrect.Rule1.aPG = m2IncorrectRule1PG;
