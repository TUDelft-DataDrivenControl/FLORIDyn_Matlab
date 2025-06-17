% a vector.
clear all
rng(1234);
N=1000;

vec=randn(1,N);
filename = 'test/randn.mat';
save(filename, 'vec')