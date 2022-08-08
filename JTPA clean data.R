library(haven)

# the original data can be found at: https://www.upjohn.org/data-tools/employment-research-data-center/national-jtpa-study
# the 1 or 2 indicates from which expfup folder the data come from

jtpapartfu1<-read_dta("mainfup1.dta")
jtpapartfu2<-read_dta("mainfup2.dta")
covariates<-read_dta("expbif.dta")
jobf1<-read_dta("bjobfup1.dta")
jobf2<-read_dta("bjobfup2.dta")
headf1<-read_dta("head_fp1.dta")
headf2<-read_dta("head_fp2.dta")

# JTPA participation

jtpa1 = jtpapartfu1
jtpa2 = jtpapartfu2
jtpa1 = jtpa1[,c("recid","f1d7")]
jtpa2 = jtpa2[,c("recid","f2d7")]
jtpa2 = jtpa2[(jtpa2[,2]== 1 | jtpa2[,2]== 2),]
jtpa1 = jtpa1[jtpa1$f1d7 %in% c(1,2),] # 1 = yes, they participated in JTPA training / 2 = no participation


# ID + covariates to save
cov_names = c("recid",# ID
              "white", # 1 = white, 0 = non-white
              "sex", #1 = male, 2 = female
              "age",
              "ra_stat", # 1 = treatment, 2 = control
              "hsged", # 1 = has high school degree or GED, 9 = unknown
              "bfmarsta" # 2 = married, 8 = missing, 9 = unknown
) 
children = apply(covariates[,c("chunder4","ch4or5", "ch6to12","chover18","ch13to18")],1,function(x){
  count = sum(as.numeric(x)) # sum of the children indicators 
  return(ifelse(count>=9,9, # if there is a 9 in the row ==> count >= 9 ==> return 9 (unknown),
                ifelse(count==0,0,1))) # otherwise check if the count is 0 (no children)
})
covariates$children = children
cov_names = c(cov_names, "children")

# filter the covariates

covariates = as.data.frame(covariates[,names(covariates) %in% cov_names])
head(covariates)

# add covariate to data

data = covariates

# add JTPA participation for each follow up interview

data = merge(data,jtpa1, all = T)
data = merge(data,jtpa2, all = T)

# transform to numeric all the data

data[] <- lapply(data, function(x) as.numeric(as.character(x)))
head(data)

# merge the JTPA participation information

data$jtpa = 2

for (i in 1:nrow(data)) {
  if (is.na(data$f1d7[i]) & is.na(data$f2d7[i])){
      data$jtpa[i] = NA
  }
  if ((!is.na(data$f1d7[i]) & data$f1d7[i]==1) | (!is.na(data$f2d7[i]) & data$f2d7[i]==1)){
    data$jtpa[i]=1
  }
}

data= subset(data, select = -c(9,10))

# only complete rows

data = data[complete.cases(data),]

# remove rows with undefined values

data = data[ (data$children != 9  &
               data$bfmarsta != 9 &
               data$bfmarsta != 8 &
               data$hsged != 9)
             ,]

head(data)

# transform variables to binary

data$married = as.numeric(data$bfmarsta==2)
data$hsged = as.numeric(data$hsged == 1)
data$male = as.numeric(data$sex == 1)
data$jtpa = as.numeric(data$jtpa == 1)
data$treatment = as.numeric(data$ra_stat == 1)
data = data[,-which(names(data) %in% c("bfmarsta","ra_stat","f1d7", "sex"))]
head(data)

# calculate Y=min(T,C)

beginendfu1<-cbind.data.frame(as.numeric(headf1$recid),headf1$f1beg,headf1$f1ref) # beginning of study and first follow up date
colnames(beginendfu1)<-c("recid","beg","end1")
datejobfu1<-cbind.data.frame(jobf1$recid,jobf1$f1b2bdd,jobf1$f1b2bmm,jobf1$f1b2byy) # date of first job (first follow up interview)

endfu2<-cbind.data.frame(as.numeric(headf2$recid), headf2$fdat_fin) # date of second follow up interview
colnames(endfu2)<-c("recid","end2")
datejobfu2<-cbind.data.frame(jobf2$recid,jobf2$f2b2bdd,jobf2$f2b2bmm,jobf2$f2b2byy) # date of first job (second follow up interview)

# look at data from first follow up interview and impute values / clean unknowns

datejobfu1[datejobfu1$`jobf1$f1b2bdd` == 92,c("jobf1$f1b2bdd")] = "01" # 92 == start of the month
datejobfu1[datejobfu1$`jobf1$f1b2bdd` == 93,c("jobf1$f1b2bdd")] = "14" # 93 == middle of the month
datejobfu1[datejobfu1$`jobf1$f1b2bdd` == 94,c("jobf1$f1b2bdd")] = "28" # 94 == end of the month

datejobfu1 = datejobfu1[!(datejobfu1$`jobf1$f1b2bdd` %in% c(98,99)),] # remove 98 = don't known, 99 = missing
datejobfu1 = datejobfu1[(datejobfu1$`jobf1$f1b2byy` != 99),] #99 = unknown

datejobfu1 = datejobfu1[datejobfu1$`jobf1$recid` %in% data$recid, ] # only the subjects in our data 
head(datejobfu1)

# transform the dates to easier format

datejobfu1$date = sapply(1:nrow(datejobfu1),function(i){
  year = datejobfu1$`jobf1$f1b2byy`[i]
  month = datejobfu1$`jobf1$f1b2bmm`[i]
  day = datejobfu1$`jobf1$f1b2bdd`[i]
  return(paste("19",year,"-",month,"-",day,sep = ""))
})

datejobfu1_clean = c()
done_recid = c()
for(i in 1:nrow(datejobfu1)){
  id = datejobfu1$`jobf1$recid`[i]
  if(!(id%in%done_recid)){
    done_recid = c(done_recid, id)
    current = datejobfu1[datejobfu1$`jobf1$recid` == id, c("date")]
    datejobfu1_clean = rbind(datejobfu1_clean, c(id,min(current)))
  }
  
}
datejobfu1_clean = as.data.frame(datejobfu1_clean)
names(datejobfu1_clean) = c("recid", "date")
head(datejobfu1_clean)

data = merge(x = data, y = datejobfu1_clean, by = "recid", all.x = TRUE)

# look at data from second follow up interview and impute values / clean unknowns

datejobfu2[datejobfu2$`jobf2$f2b2bdd` == 92,c("jobf2$f2b2bdd")] = "01"# 92 == start of the month
datejobfu2[datejobfu2$`jobf2$f2b2bdd` == 93,c("jobf2$f2b2bdd")] = "14"# 93 == middle of the month
datejobfu2[datejobfu2$`jobf2$f2b2bdd` == 94,c("jobf2$f2b2bdd")] = "28"# 93 == end of the month

datejobfu2 = datejobfu2[!(datejobfu2$`jobf2$f2b2bdd` %in% c(98,99)),] # remove 98 = don't known, 99 = missing
datejobfu2 = datejobfu2[(datejobfu2$`jobf2$f2b2byy` != 99),] # 99 = unknown

# only the subjects in our data that have NA as date, meaning they did not have a job when the first follow up interview was taken

done_recid = c()
datejobfu2_cleaned = c()
for(i in 1:nrow(datejobfu2)){
  id = datejobfu2$`jobf2$recid`[i]
  if( !(id %in% done_recid) ){
    if( (id %in% data$recid) && is.na(data[data$recid == id, c("date")])){
      datejobfu2_cleaned = rbind(datejobfu2_cleaned, datejobfu2[i,])
    }
  }
}
datejobfu2_cleaned = as.data.frame(datejobfu2_cleaned)
names(datejobfu2_cleaned) = names(datejobfu2)

# transform the dates to different format

datejobfu2_cleaned$date = sapply(1:nrow(datejobfu2_cleaned),function(i){
  year = datejobfu2_cleaned$`jobf2$f2b2byy`[i]
  month = datejobfu2_cleaned$`jobf2$f2b2bmm`[i]
  day = datejobfu2_cleaned$`jobf2$f2b2bdd`[i]
  return(paste("19",year,"-",month,"-",day,sep = ""))
})

datejobfu2_cleaned2 = c()
done_recid = c()
for(i in 1:nrow(datejobfu2_cleaned)){
  id = datejobfu2_cleaned$`jobf2$recid`[i]
  if(!(id%in%done_recid)){
    done_recid = c(done_recid, id)
    current = datejobfu2_cleaned[datejobfu2_cleaned$`jobf2$recid` == id, c("date")]
    datejobfu2_cleaned2 = rbind(datejobfu2_cleaned2, c(id,min(current)))
  }
  
}
datejobfu2_cleaned2 = as.data.frame(datejobfu2_cleaned2)
names(datejobfu2_cleaned2) = c("recid", "date")
head(datejobfu2_cleaned2)

for(i in 1:nrow(data)){
  id = data$recid[i]
  if(is.na(data$date[i])){
    if(id %in% datejobfu2_cleaned2$recid){
      data[i,c("date")] = datejobfu2_cleaned2[datejobfu2_cleaned2$recid == id, c("date")]
    }
  }
}

sum(is.na(data$date))/nrow(data) # gives the censoring percentage

data$delta = as.numeric(!is.na(data$date)) # assigns the censoring indicator

# assigning the censoring date for the participants of follow up 2

months = as.data.frame(rbind(c("JAN","01"),
                             c("FEB","02"),
                             c("MAR","03"),
                             c("APR","04"),
                             c("MAY","05"),
                             c("JUN","06"),
                             c("JUL","07"),
                             c("AUG","08"),
                             c("SEP","09"),
                             c("OCT","10"),
                             c("NOV","11"),
                             c("DEC","12")))
names(months) = c("name", "number")
endfu2$date = sapply(1:nrow(endfu2),function(i){
  str =  endfu2$end2[i]
  if(str == ""){
    return(NA)
  }
  day = substr(str,1,2)
  month = months[months$name == substr(str, 3,5), c("number") ]
  year = substr(str, 6,7)
  return(paste("19",year,"-",month,"-",day,sep = ""))})

for(i in 1:nrow(data)){
  if(data$delta[i] == 0){
    id = data$recid[i]
    if(id %in% endfu2$recid){
      data[i,c("date")] = endfu2$date[endfu2$recid == id]
    }
  }
}

# assigning the censoring date for the participants of follow up 1

beginendfu1$date_end = sapply(1:nrow(beginendfu1),function(i){
  str = beginendfu1$end1[i]
  return(paste("19",substr(str,1,2),"-",substr(str,3,4),"-",substr(str,5,6),sep = ""))
})

beginendfu1$date_beg = sapply(1:nrow(beginendfu1),function(i){
  str = beginendfu1$beg[i]
  return(paste("19",substr(str,1,2),"-",substr(str,3,4),"-",substr(str,5,6),sep = ""))
})

for(i in 1:nrow(data)){
  if(is.na(data$date[i])){
    id = data$recid[i]
    if(id %in% beginendfu1$recid){
      data[i,c("date")] = beginendfu1$date_end[beginendfu1$recid == id]
    }
  }
}

# assigning the date of treatment assignment

data$date_start = NA
for(i in 1:nrow(data)){
  id = data$recid[i]
  data[i,c("date_start")] = beginendfu1[beginendfu1$recid == id, c("date_beg")]
}

# calculate the number of days between treatment assignment and finding a job / being censored

data$date=as.Date.character(data$date)
data$date_start=as.Date.character(data$date_start)

data$days =  sapply(1:nrow(data), function(i){
  return(difftime(data$date[i],data$date_start[i], units = "days"))
})

# only the people who start unemployed

data = data[data$days>0,]
data=na.omit(data)
data=subset(data, select = -c(10,12))

# first look at the clean data

head(data)
par(mfrow=c(1,1))
hist(data$days, breaks = sqrt(nrow(data)), freq = F, ylim = c(0,0.009))

par(mfrow=c(2,1))
hist(data[data$delta == 1, c("days")], freq = F, breaks = sqrt(sum(as.numeric(data$delta == 1))))
hist(data[data$delta == 0, c("days")], freq = F, breaks = sqrt(sum(as.numeric(data$delta == 0))))

write.csv(data,paste("clean_dataset_JTPA.csv",sep = ""), row.names = F)


