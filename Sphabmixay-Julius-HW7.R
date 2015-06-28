##### PROBLEM 1 #####

# The length of confidence intervals for the smokers dataset is 
# the same for all possible combinations, however their
# differences in means are not.

load("smokers.rda");

permHSD = function(y,g,level=.95,B=2000) {

	data = data.frame(y,g);
	Qb = numeric(0);

	obsNAMES = unique(g);
	obsMeans = computeMeans(y,g);
	obsCOUNTS = findGroupCounts(g);

	obs.S = compute.S(obsNAMES,obsCOUNTS,obsMeans);
	
	for (b in 1:B) {
		
		gPermutation = sample(g);
		permTable = table(y,gPermutation);
		means = computeMeans(permTable[,1],gPermutation);
		NAMES = unique(gPermutation);
		COUNTS = findGroupCounts(gPermutation);

		S = compute.S(NAMES,COUNTS,means);

		T = numeric(0);
		D = numeric(0);
		counter = 0;
		i = 1;
		while (i <= (length(NAMES) - 1)) {
			j = i + 1;
			while( j <= length(NAMES)) {
				diff = abs(means[i] - means[j]);
				d = sqrt((1/COUNTS[i])+(1/COUNTS[j]));
				T = append(T, (diff/(S*d)));
				D = append(D,d);
				counter = counter+1;
				j = j+1;
			}
			i = i+1;
		}
		Q = max(T);
		Qb = append(Qb,Q);
	}
	q = quantile(Qb,(1-level));
	values = numeric(0);
	for (i in 1:counter) {
		v =(q * obs.S * D[i]);
		values = append(values,v);
	}
	UPPER = numeric(0);
	LOWER = numeric(0);
	for (i in 1:(length(NAMES)-1)) {
		for (j in (i+1):length(NAMES)) {
			upper = obsMeans[i] - obsMeans[j] + values[j];
			lower = obsMeans[i] - obsMeans[j] - values[j];
			UPPER = append(UPPER,upper);
			LOWER = append(LOWER,lower);
		}
	}
	cat("The lower bounds for the CI's are:", "\n");
	print(LOWER);
	cat("\n");
	cat("The upper bounds for the CI's are:", "\n");
	print(UPPER);
}
# This is a helper method
computeMeans = function(w,g) {
	
	Names = unique(g);
	AVG = numeric(0);
	for (i in 1:length(Names)) {
		sum = 0;
		count = 0;
		for (j in 1:length(w)) {
			if (g[j] == Names[i]) {
				sum = sum + w[j];
				count = count + 1;
			}
		}
		AVG = append(AVG,(sum/count));
		AVG = unname(AVG);
	}
	return(AVG);
}
# This is a helper method
findGroupCounts = function(g) {

	Names = unique(g);
	groupCounts = numeric(0);
	for (i in 1:length(Names)) {
		count = 0;
		for (j in 1:length(g)) {
			if (g[j] == Names[i]) {
				count = count + 1;
			}
		}
		groupCounts = append(groupCounts,count);
		groupCounts = unname(groupCounts);
	}
	return(groupCounts);
}
# This is a helper method
compute.S = function(NAMES,COUNTS,means) {

	k = 0;
	Y = 0;
	for (i in 1:length(NAMES)) {
		for (j in (1+k):(COUNTS[i]+k)) {
			Y = Y + (y[j] - means[i])^2;
		}
		k = COUNTS[i];
	}
	S.Squared = Y/(length(y) - length(NAMES));
	S = sqrt(S.Squared);
	return(S);
}

y = c(smokers[,1],smokers[,2],smokers[,3],smokers[,4]);
g = c(rep("non",6),rep("light",6),rep("moderate",6),rep("heavy",6));

permHSD(y,g);

##### PROBLEM 2 #####

allpairs.t.test = function(y,g,level) {
	
	NAMES = unique(g);
	COUNTS = findGroupCounts;
	means = computeMeans(y,g);
	data = findGroupTable(y,g,NAMES);
	p.values = numeric(0);
	
	i = 2;
	count = 0;
	while (i <= length(NAMES)) {
		j = i + 1;
		while(j <= length(NAMES)+1) {
			p = t.test(data[,i],data[,j],paired=TRUE)$p.value
			p.values = c(p.values,p);
			count = count + 1;
			j = j+1;
		}
		i = i+1;
	}
	p.values = sort(p.values);
	R = 1
	for (i in 1:length(p.values)) {
		ratio = level/(count + 1 - i);
		if (p.values[i] > ratio) {
			R = i;
			break;
		}
	}
	adj.p.values = numeric(0);
	
	for (i in 1:length(p.values)) {
		values = numeric(0);
		mins = numeric(0);
		for (j in 1:i) {
			values = append(values, (count - j + 1)*p.values[j]);
		}
		for (k in 1:i) {
			mins = append(mins, min(values[k],1));
		}
		adj.p.values = append(adj.p.values, max(mins));
	}
	print(adj.p.values);
	return(adj.p.values);
}
# This is a helper method
findGroupTable = function(y,g,NAMES) {
	
	zero = rep(0,length(y));
	data = data.frame(zero);
	for (i in 1:length(NAMES)) {
		z = numeric(0);
		for (j in 1:length(g)) {
			if (g[j] == NAMES[i]) {
				z = append(z,y[j]);
			}
		}
		data = data.frame(data,z);
	}
	return(data);
}
cat("\n")
allpairs.t.test(y,g,level=.95);
result = pairwise.t.test(y,g);
print(result);

