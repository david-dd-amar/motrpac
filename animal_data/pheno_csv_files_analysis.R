# Set the working directory to the folder with the data
setwd("/Users/David/Desktop/MoTrPAC/march_2019/pheno_data/PASS1A.6M-RLS 0,01/3-Data Sets/")
all_csvs = list.files(".") # get all files in dir
# read all files
csv_data = list()
for(fname in all_csvs){
  csv_data[[fname]] = read.csv(fname,stringsAsFactors = F)
}

ac_test_data = csv_data[[which(grepl("Acute.Test",names(csv_data)))]]
dim(ac_test_data)
# convert the shock lengths to numbers (seconds)
parse_shocktime<-function(x){
  arr = strsplit(x,split=":")[[1]]
  if(length(arr)<2){return(NA)}
  return(as.numeric(arr[1])*60+as.numeric(arr[2]))
}
tmp_x = ac_test_data$howlongshock
tmp_x = sapply(tmp_x, parse_shocktime)
ac_test_data$howlongshock = tmp_x

# check the time differences between start and end
test_times = as.difftime(ac_test_data$t_complete) - as.difftime(ac_test_data$t_start)
# table of the values: all except for on are 0.5 hours
table(test_times)
# Get the comment of the sample that is not 0.5h
ac_test_data[test_times!=0.5,"comments"]

# histogram of distances
hist(ac_test_data$distance,col="blue",breaks=50,main = "Histogram of distances")

# Correlation between distance and number of shocks
# Get the indices of the samples with shock information -
# these the animals that did the acute test
timesshock_inds = !is.na(ac_test_data$timesshock)
# create a new dataframe with the selected animals
trained_animals_data = ac_test_data[timesshock_inds,]
sp_corr = cor(trained_animals_data$distance,
              trained_animals_data$timesshock,method="spearman")
plot(trained_animals_data$distance,trained_animals_data$timesshock,
     main=paste("Distances vs times shock given, rho=",format(sp_corr,digits = 3),sep=""),
     pch=20,ylab="Times shock given",xlab="Distance")

# A "smarter" analysis: regression of the distance using shock info
dist_lm  = lm(distance~timesshock+howlongshock+weight+days_start,data=trained_animals_data)
summary(dist_lm)
# We have some clear outliers
plot(dist_lm$residuals,main="Linear regression, raw residuals",ylab="residual")
library(MASS)
plot(studres(dist_lm),main="Linear regression, studentized residuals",ylab="residual")
outliers = abs(studres(dist_lm)) > 2
table(outliers)
trained_animals_data[outliers,"comments"]
plot(dist_lm$fitted.values,trained_animals_data$distance,lwd=2,
     main="Fitted vs real values",ylab="Distances",xlab="Fitted distances")
abline(0,1,col="red",lty=2,lwd=3)


# Load additional information about the animals
registr_data = csv_data[[which(grepl("Regist",names(csv_data)))]]
rownames(registr_data) = as.character(registr_data$participantGUID)
# make the rownames in the test data comparable
rownames(trained_animals_data) = trained_animals_data$participantGUID
# add sex to the trained animal data data frame
sex_key = c("Female","Male")
trained_animals_data$sex = sex_key[registr_data[rownames(trained_animals_data),"sex"]]

# Map site Ids to their names
site_names = c("910"="Joslin","930"="Florida")
trained_animals_data$site = site_names[as.character(trained_animals_data$siteID)]

# Sanity check: the numbers should be the same for both sites
table(ac_test_data$siteID)
table(trained_animals_data$site,trained_animals_data$sex)

run_wilcox<-function(x1,x2){
  return(wilcox.test(x1[x2==x2[1]],x1[x2!=x2[1]])$p.value)
}
# Compare the distances, shocks, and weight
par(mfrow=c(1,3),mar=c(10,4,4,4))
# Site only
p_dist = run_wilcox(trained_animals_data$distance,trained_animals_data$site)
boxplot(distance~site,data=trained_animals_data,col="cyan",ylab="Distance",
        main=paste("Site vs. distance, p<",format(p_dist,digits = 2)),
        cex.main=1,las=2)
p_timesshock = run_wilcox(trained_animals_data$timesshock,trained_animals_data$site)
boxplot(timesshock~site,data=trained_animals_data,col="red",ylab="Times shocked",
        main=paste("Site vs. times shocked, p<",format(p_timesshock,digits = 3)),
        cex.main=1,las=2)
p_w = run_wilcox(trained_animals_data$weight,trained_animals_data$site)
boxplot(weight~site,data=trained_animals_data,col="cyan",ylab="Weight",
        main=paste("Site vs. weight, p=",format(p_w,digits = 2)),
        cex.main=1,las=2)
# Site and sex
par(mfrow=c(1,3),mar=c(10,4,4,4))
boxplot(distance~site+sex,data=trained_animals_data,col="cyan",ylab="Distance",
        main="Site vs. distance",cex.main=1,las=2)
boxplot(timesshock~site+sex,data=trained_animals_data,col="red",ylab="Times shocked",
        main="Site vs. times shocked",cex.main=1,las=2)
boxplot(weight~site+sex,data=trained_animals_data,col="cyan",ylab="Weight",
        main="Site vs. weight",cex.main=1,las=2)

# Regress time shocked and distance vs. site and sex
summary(lm(timesshock~site+sex,data=trained_animals_data))
summary(lm(distance~site+sex,data=trained_animals_data))

