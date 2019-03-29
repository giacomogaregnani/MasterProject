load ahom8;
for i=1:2^8+1
	seq8(i,:,:) = seq{i}(size(seq{i},1),:,:);
end
load ahom9;
for i=2:2:2^9
	seq9(i,:,:) = seq{i}(size(seq{i},1),:,:);
end
load ahom10;
for i=2:2:2^10
	seq10(i,:,:) = seq{i}(size(seq{i},1),:,:);
end

seq9(1:2:2^9+1,:,:) = seq8;
seq10(1:2:2^10+1,:,:) = seq9;
tensor = seq10;
clear seq8 seq9 seq10;

% R = (0:512)/512;

av = (tensor(:,1,2) + tensor(:,2,1)) / 2;
tensor(:,1,2) = av;
tensor(:,2,1) = av;

M = 2^10;
tensor(1:4*M+1,1,1) = ...
	[tensor(1:M,1,1);
	tensor(M+1:-1:2,2,2); 
	tensor(1:M,2,2); 
	tensor(M+1:-1:1,1,1)];
tensor(1:4*M+1,2,2) = ...
	[tensor(1:M,2,2); 
	tensor(M+1:-1:2,1,1); 
	tensor(1:M,1,1); 
	tensor(M+1:-1:1,2,2)];
tensor(1:4*M+1,1,2) = ...
	[tensor(1:M,1,2);
	tensor(M+1:-1:2,1,2);
	-tensor(1:M,1,2);
	-tensor(M+1:-1:1,1,2)];
tensor(1:4*M+1,2,1) =tensor(1:4*M+1,1,2);

R = (0:4*M)/(4*M);

 for i=1:2
 	for j=1:2
 		cs{i,j} = spline(R,tensor(:,i,j));
 	end
 end
 save('ahom_expe.mat','cs','R','tensor');