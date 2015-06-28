
options(warn = -1); # suppress warnings
##### PROBLEM 1 #####
M = seq(0,1,.1);

n1 = 100;
n2 = 1000;
A = 200;
alpha = .05;
T1.vector = numeric(length(M));
T2.vector = numeric(length(M));
W1.vector = numeric(length(M));
W2.vector = numeric(length(M));

for (i in 1:length(M)) {
	
	T1 = numeric(A);
	W1 = numeric(A);
	for (j in 1:A) {
		
		X = rnorm(n1,M[i],1);
		t.pValue = t.test(X, alternative="greater",var.equal=TRUE)$p.value;
		w.pValue = wilcox.test(X, alternative="greater")$p.value;
		
		T1[j] = t.pValue;
		W1[j] = w.pValue;
	}
	T1.ratio = sum(T1 < alpha)/A;
	W1.ratio = sum(W1 < alpha)/A;
	
	T1.vector[i] = T1.ratio;
	W1.vector[i] = W1.ratio;
}

for (i in 1:length(M)) {
	
	T2 = numeric(A);
	W2 = numeric(A);
	for (j in 1:A) {
		
		X = rnorm(n2,M[i],1);
		t.pValue = t.test(X, alternative="greater",var.equal=TRUE)$p.value;
		w.pValue = wilcox.test(X, alternative="greater")$p.value;
		
		T2[j] = t.pValue;
		W2[j] = w.pValue;
	}
	T2.ratio = sum(T2 < alpha)/A;
	W2.ratio = sum(W2 < alpha)/A;
	
	T2.vector[i] = T2.ratio;
	W2.vector[i] = W2.ratio;
}

par(ask=TRUE);
plot(seq(0,1,.1),seq(0,1,.1),type="n",xlab = "M",ylab = "Power");

lines(M,T1.vector,col="red",lwd=1.5);
lines(M,W1.vector,col="blue",lwd=1.5);

lines(M,T2.vector,col="green",lwd=1.5);
lines(M,W2.vector,col="purple",lwd=1.5);

legend(.4,.5,c("t-test,100","signed-r test,100","t-test,1000","sign-test,1000"),
		lty=c(1,1),lwd=c(1.5,1.5),col=c("red","blue","green","purple"));
		
# We can see that the t-test and signed-rank test
# power curves are both similar. Thus, the signed-rank
# test is good for estimating power as well.
# If the sample size is 100, we can see that the power
# curve increases from M=0 to M=.5 .
# If the sample size is 1000, we can see that the power
# curve increases rapidly from M=0 to M=.2 .

##### PROBLEM 2 #####

# Please use the data files to make the contingency tables.
study1 = as.matrix(read.table("study1.txt",sep=",",
			row.names=c("recurrence","nonrecurrence"),
			col.names=c("treatment","control")));
study2 = as.matrix(read.table("study2.txt",sep=",",
			row.names=c("recurrence","nonrecurrence"),
			col.names=c("treatment","control")));
study3 = as.matrix(read.table("study3.txt",sep=",",
			row.names=c("recurrence","nonrecurrence"),
			col.names=c("treatment","control")));
study4 = as.matrix(read.table("study4.txt",sep=",",
			row.names=c("recurrence","nonrecurrence"),
			col.names=c("treatment","control")));
study5 = as.matrix(read.table("study5.txt",sep=",",
			row.names=c("recurrence","nonrecurrence"),
			col.names=c("treatment","control")));
study6 = as.matrix(read.table("study6.txt",sep=",",
			row.names=c("recurrence","nonrecurrence"),
			col.names=c("treatment","control")));
study7 = as.matrix(read.table("study7.txt",sep=",",
			row.names=c("recurrence","nonrecurrence"),
			col.names=c("treatment","control")));
study8 = as.matrix(read.table("study8.txt",sep=",",
			row.names=c("recurrence","nonrecurrence"),
			col.names=c("treatment","control")));

### PART A ###

# Null hypothesis: The proportions of recurrence is the
# same for both the Treatment and Control groups.
# Alternative: The proportion of recurrence in the
# Treatment group is less than the Control group.

cat("p-values obtained by PROPORTION TEST","\n");
print(study1);

study1.p = prop.test(study1,alternative="l")$p.value;
study2.p = prop.test(study2,alternative="l")$p.value;
study3.p = prop.test(study3,alternative="l")$p.value;
study4.p = prop.test(study4,alternative="l")$p.value;
study5.p = prop.test(study5,alternative="l")$p.value;
study6.p = prop.test(study6,alternative="l")$p.value;
study7.p = prop.test(study7,alternative="l")$p.value;
study8.p = prop.test(study8,alternative="l")$p.value;

p.vector1 = c(study1.p,study2.p,study3.p,study4.p,
				study5.p,study6.p,study7.p,study8.p);	
				
# Because all the studies are performed by different researchers,
# we can assume the studies are independent of each other.
# We can apply the Hochberg procedure.

alpha = .05;
m = 8
r.vector = numeric(0);
index.vector = numeric(0);
sorted.p.vector = sort(p.vector1);
R = 0;
cat("For the Hochberg procedure:","\n");
for (r in 1:m) {
	print(c(sorted.p.vector[r],alpha/(m+1-r)));
	if (sorted.p.vector[r] <= (alpha/(m+1-r))) {
		r.vector = append(r.vector,sorted.p.vector[r]);
		index.vector = append(index.vector,r);
	}
}
if (length(index.vector) > 0) {
	R = max(index.vector);
	cat("Reject all tests associated with p-values
								from index 1 to ", R);
} else {
	cat("Reject none of the tests because R = 0.","\n");
}
cat("R is", R, "\n");

# Try the Holm procedure because it is more powerful
# than the Bonferroni procedure and has
# the same assumptions as Bonferroni (independence doesn't matter).

cat("\n");
cat("For the Holm procedure:","\n");
r.vector = numeric(0);
index.vector = numeric(0);
for (k in 1:(m)) {
	print(c(sorted.p.vector[k],alpha/(m-k+1)));
	if (sorted.p.vector[k] > (alpha/(m-k+1))) {
		r.vector = append(r.vector,sorted.p.vector[k]);
		index.vector = append(index.vector,k);
	}
}
R = min(index.vector);
if (1 %in% index.vector) {
	print(index.vector);
	cat("Reject none of the tests because R = 1.","\n");
} else {
	cat("Reject all tests associated with p-values
			from index 1 to ", (R - 1),"\n");
}
cat("R is", R, "\n");
cat("\n");

# For both the Holm and Hochberg procedure, all studies
# are considered for the multiple testing.

### PART B ###

# Null hypothesis: Control and Treatment are equally
# likely to have recurrence.
# Alternative: The Treatment group has less recurrence
# than the Control group.

cat("p-values obtained by FISHER'S EXACT TEST","\n");
p.vector2 = c(fisher.test(study1,alternative="l")$p.value,
				fisher.test(study2,alternative="l")$p.value,
				fisher.test(study3,alternative="l")$p.value,
				fisher.test(study4,alternative="l")$p.value,
				fisher.test(study5,alternative="l")$p.value,
				fisher.test(study6,alternative="l")$p.value,
				fisher.test(study7,alternative="l")$p.value,
				fisher.test(study8,alternative="l")$p.value);

# Try the Hochberg procedure.

alpha = .05;
m = 8
r.vector = numeric(0);
index.vector = numeric(0);
sorted.p.vector = sort(p.vector2);
R = 0;
cat("For the Hochberg procedure:","\n");
for (r in 1:m) {
	print(c(sorted.p.vector[r],alpha/(m+1-r)));
	if (sorted.p.vector[r] <= (alpha/(m+1-r))) {
		r.vector = append(r.vector,sorted.p.vector[r]);
		index.vector = append(index.vector,r);
	}
}
if (length(index.vector) > 0) {
	R = max(index.vector);
	cat("Reject all tests associated with p-values
								from index 1 to ", R);
} else {
	cat("Reject none of the tests because R = 0.","\n");
}
cat("R is", R, "\n");

# Try the Holm procedure.

cat("\n");
cat("For the Holm procedure:","\n");
r.vector = numeric(0);
index.vector = numeric(0);
for (k in 1:(m)) {
	print(c(sorted.p.vector[k],alpha/(m-k+1)));
	if (sorted.p.vector[k] > (alpha/(m-k+1))) {
		r.vector = append(r.vector,sorted.p.vector[k]);
		index.vector = append(index.vector,k);
	}
}
R = min(index.vector);
if (1 %in% index.vector) {
	print(index.vector);
	cat("Reject none of the tests because R = 1.","\n");
} else {
	cat("Reject all tests associated with p-values
			from index 1 to ", (R - 1),"\n");
}
cat("R is", R, "\n");
cat("\n");

# For both the Holm and Hochberg procedure, all studies
# are considered for the multiple testing.

### PART C ###

fisher.multiple.test = function(p) {
	
	F = -2*sum(log(p));
	pValue = pchisq(F,df=2*length(p), lower.tail = FALSE);
	return(pValue);
}

result1 = fisher.multiple.test(p.vector1);
result2 = fisher.multiple.test(p.vector2);
cat("p-value for part A using fisher's multiple test","\n");
print(result1);
cat("p-value for part B using fisher's multiple test","\n");
print(result2);

# Null hypothesis: All studies are true.
# Alternative: At least one study is false.
# The p-values for both A and B are high. Thus, we can conclude from
# the studies that all null hypotheses for the studies are true.
# This makes sense because the rate of recurrence in the Treatment
# is mostly higher than the rate in the Control if you look
# at the studies.

##### PROBLEM 3 #####

# Please use the nottem.txt file
Data = read.table("nottem.txt");
Time = Data[,1];
Temp = Data[,2];

xc1 = cos(pi*Time/6);
xs1 = sin(pi*Time/6);

xc2 = cos(2*pi*Time/6);
xs2 = sin(2*pi*Time/6);

xc3 = cos(3*pi*Time/6);
xs3 = sin(3*pi*Time/6);

data1.lm = lm(Temp~xc1+xs1);
data1.fit = fitted(data1.lm);

data2.lm = lm(Temp~xc1+xs1+xc2+xs2);
data2.fit = fitted(data2.lm);

data3.lm = lm(Temp~xc1+xs1+xc2+xs2+xc3+xs3);
data3.fit = fitted(data3.lm);

print(anova(data3.lm));

plot(Temp ~ Time, data = Data, xlim=c(1,240));

lines(data1.fit,col="blue");

lines(data2.fit,col="red");

lines(data3.fit,col="green");

# It is clear that we should choose m = 1 because having more
# coeffcients (m > 1) will overfit the data. Thus, we should
# use data1.fit (blue curve) to fit the nottem data.

##### PROBLEM 4 #####

permLAR.test = function(x,y,B=2000) {

	data.rq = rq(y~x);
	N = length(x);
	Beta = unname(data.rq$coefficients);
	
	SAR.vector = numeric(N);
	SAF.vector1 = numeric(N);
	SAF.vector2 = numeric(N);
	for (i in 1:N) {
		
		SAR.vector[i] = abs(y[i] - Beta[1]*1);
		SAF.vector[i] = abs(y[i] - (Beta[1]*1) - (Beta[2]*x[i]));
		
	}
	SAR = sum(SAR.vector);
	SAF = sum(SAF.vector);
	
	T.obs = (SAR - SAF)/SAF;
	
	Tb.vector = numeric(B);
	for (b in 1:B) {
		
		perm.SAR.vector = numeric(N);
		perm.SAF.vector = numeric(N);
		perm.y = sample(y,N,sample=FALSE);
		perm.data.rq = rq(perm.y~x);
		perm.Beta = unname(perm.data.rq$coefficients);
		for (i in 1:N) {
		
			perm.SAR.vector[i] = abs(perm.y[i] - (perm.Beta[1]*1));
			perm.SAF.vector[i] = abs(perm.y[i] - (perm.Beta[1]*1) - (perm.Beta[2]*x[i]));
		}
		
		perm.SAR = sum(perm.SAR.vector);
		perm.SAF = sum(perm.SAF.vector);
		
		Tb.vector[b] = (perm.SAR - perm.SAF)/perm.SAF;
	}
	p = (sum(Tb.vector >= T.obs) + 1)/(B + 1);
	print(p);
	return(p);
}