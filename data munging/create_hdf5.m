%create an hdf5 file from monkey data
%first create the group structure

fid = H5F.create('mGoodStableRule1PingRej.h5');
plist = 'H5P_DEFAULT';
gid1 = H5G.create(fid,'monkey1',plist,plist,plist);
gid2 = H5G.create(fid,'monkey2',plist,plist,plist);
gid11 = H5G.create(gid1,'correct',plist,plist,plist);
gid12 = H5G.create(gid1,'incorrect',plist,plist,plist);
gid21 = H5G.create(gid2,'correct',plist,plist,plist);
gid22 = H5G.create(gid2,'incorrect',plist,plist,plist);
%monkey 1 days 1 - 23
gid111 = H5G.create(gid11,'day1',plist,plist,plist);
gid112 = H5G.create(gid11,'day2',plist,plist,plist);
>> H5G.close(gid1);
>> H5G.close(gid2);
>> H5G.close(gid);
>> H5F.close(fid);
h5disp('mGoodStableRule1PingRej.h5')
%1600 is largest days' trials, 1811 timepoints
mydata = monkey(1).day(1).correct.chan1PEC;
h5create('mGoodStableRule1PingRej.h5','/monkey1/correct/day1/chan1PEC',size(mydata))
h5write('mGoodStableRule1PingRej.h5', '/monkey1/correct/day1/chan1PEC', mydata)