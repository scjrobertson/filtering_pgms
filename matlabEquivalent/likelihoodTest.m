clear all;

load('difference.mat');
load('ZInv.mat');

zDimension = 2;
targetNumber = 2;
numberOfMeasurements = 518;

tic;
ZInvStack = reshape(permute(ZInv, [1 3 2]), [zDimension*zDimension targetNumber]);
leftProduct = ZInvStack*difference;
exponentialComponent = reshape(sum(difference.*leftProduct(1:2, :), 1), [1 numberOfMeasurements targetNumber]);
toc;


tic;
diffRow = reshape(difference, [1 zDimension targetNumber*numberOfMeasurements]);
diffCol = reshape(difference, [zDimension 1 targetNumber*numberOfMeasurements]);
left = multiprod( multiprod(diffRow, ZInv), diffCol);

toc;
