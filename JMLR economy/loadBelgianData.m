% File to process Labour.csv by extracting only the columns of interest,
% normalizing and then removing some observed outliers which perturb the
% computations.

filename='Labour.csv';
dataOrig = csvread(filename,1,1);
data=dataOrig;
data(:,3)=-log(data(:,3));
data(data(:,1)>100,:)=[];
data(data(:,2)>1500,:)=[];
data(data(:,4)>100,:)=[];
dataNorm=normalize(data);
figure
scatter3(data(:,1),data(:,2),data(:,3))
figure
scatter3(dataNorm(:,1),dataNorm(:,2),dataNorm(:,3))
data=[dataNorm(:,1:3)]; %capital, labour, output, normalized
% save('dataLabourNorm', 'data');

data=data(data(:,1)<4,:);
data=data(data(:,2)<4,:);
% save('dataLabourNorm-crop', 'data');