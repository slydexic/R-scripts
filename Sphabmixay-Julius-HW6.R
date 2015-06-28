##### PROBLEM 1 #####

median.test = function(x,y,B=2000) {

	m = median(c(x,y));
	a11 = sum(x <= m);
	a12 = sum(y <= m);
	a21 = sum(x > m);
	a22 = sum(y > m);
	u = c(rep("<= m",a11),rep("<= m",a12),rep("> m",a21),rep("> m",a22));
	v = c(rep("X",a11),rep("Y",a12),rep("X",a21),rep("Y",a22));
	obsTable = table(u,v);
	obs = chisq.test(obsTable)$statistic;
	print(obsTable);
	
	testStat = numeric(B);
	for (i in 1:B) {

		vPermutation = sample(v);
		permTable = table(u,vPermutation);
		testStat[i] = chisq.test(permTable)$statistic;
	}
	p.value = (sum(testStat >= obs) + 1)/(B + 1);
	print(p.value);
	return(p.value);
	
}

load("cloudseeding.rda");
median.test(cloudseeding[,1],cloudseeding[,2]);

##### PROBLEM 2 #####

load("smokers.rda");

# The null distribution is a chi-square distribution
# with 1 degree of freedom.

permF.test = function(y,g,B=2000) {

	y.anova = lm(y ~ g)
	obs = summary(y.anova)$fstatistic[1];
	x = numeric(B);
	
	for (i in 1:B) {
		
		permutation = sample(g);
		perm.anova = lm(y ~ permutation);
		x[i] = summary(perm.anova)$fstatistic[1];
	}
	hist(x);
	p.value = (sum(x >= obs) + 1)/(B + 1);
	print(p.value);
	return(p.value);

}

z = c(smokers[,1],smokers[,2],smokers[,3],smokers[,4]);
group1 = rep(c("non"),6);
group2 = rep(c("light"),6);
group3 = rep(c("moderate"),6);
group4 = rep(c("heavy"),6);
categories = c(group1,group2,group3,group4);

permF.test(z,categories);