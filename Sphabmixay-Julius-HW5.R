##### PROBLEM 1  #####
load("balance.rda");

boot2.t.test = function(x,y,B=2000) {

	T = t.test(x,y)$statistic;
	meanX = mean(x);
	meanY = mean(y);
	
	z = numeric(B);
	for (i in 1:B) {

		u = sample(x, length(x), replace=TRUE);
		v = sample(y, length(y), replace=TRUE);
		denominator = sqrt((var(u)/length(u))+(var(v)/length(v)));
		Tb = ((mean(u)-mean(v)) - (meanX-meanY))/denominator;
		z[i] = abs(Tb);
	}
	I = 1*(z >= T);
	pVal = (sum(I)+1)/(B+1);
	print(pVal);
	return(pVal);
}
boot2.t.test(balance[,2],balance[,3]);

##### PROBLEM 2 #####

boot2.t.CI = function(x,y, conf=0.95,B=2000) {

	T = t.test(x,y)$statistic;
	diff.m = mean(x)-mean(y);
	D = sqrt((var(x)/length(x))+(var(y)/length(y)));
	
	z = numeric(B);
	for (i in 1:B) {

		u = sample(x, length(x), replace=TRUE);
		v = sample(y, length(y), replace=TRUE);
		denominator = sqrt((var(u)/length(u))+(var(v)/length(v)));
		Tb = ((mean(u)-mean(v)) - (mean(x)-mean(y)))/denominator;
		z[i] = Tb;
	}

	alpha = 1 - conf;
	lower = diff.m - quantile(z,alpha/2)*D;
	upper = diff.m + quantile(z,alpha/2)*D;
	CI = c(upper,lower);
	names(CI) = c("2.5%","97.5%")
	print(CI);

	lower.p.value = mean(z <= T);
	upper.p.value = mean(z > T);
	p.value = 2*min(lower.p.value,upper.p.value);
	print(p.value);
	return(p.value);
}
boot2.t.CI(balance[,2],balance[,3]);