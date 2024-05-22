library(dplyr)
# Used for pulling all the Kraken data together and estimating correlation with Alzheimer's

# gEt the annotation file
metad <- read.table('onecohort.csv', sep = ',', header=TRUE)
virus_reads <- read.table('virus_report.virus_report', sep = '\t', header=TRUE)

# mutate columns
virus_reads <- virus_reads %>% as_tibble() %>%
    mutate(file= tools::file_path_sans_ext(basename(file))) %>%
    rename(percent= X.)

# remove empty rna_id
metad <- metad[metad$rna_id!="",]
joined <- merge(metad, virus_reads, by.x="rna_id", by.y="file", all.x=TRUE)

joined$v_present <- "yes"
joined[is.na(joined$reads),'v_present'] <- 'no' 

# Divide numbers in this file by 4
size_of_files <- read.table('/data/CARD_AA/users/willeyml/ViralRNA_Blast/Files/FastaCount.csv', sep=',', header=FALSE)
colnames(size_of_files) <- c('reads', 'file')
size_of_files['reads'] <- size_of_files['reads'] / 4
size_of_files['reads'] <- size_of_files['reads'] + 1

# Test 
log( 26 / (size_of_files[size_of_files$file == '01_120405',]['reads']+1))

size_of_library <- merge(size_of_files, joined, by.x='file', by.y='rna_id')
size_of_library['realperclog'] <- log(size_of_library['reads.y'] / size_of_library['reads.x'])

# By percent of reads
t.test(formula=percent~clinAD, data=joined)
# By read count
t.test(formula=reads~clinAD, data=joined)
# By whether any virus or none
chisq.test(joined$v_present, joined$clinAD, correct=FALSE)

joined <- within(joined, {coglev.ct <- C(coglev.f, sum)
       print(attributes(coglev.ct))
})

# By percent of reads
t.test(formula=realperclog~clinAD, data=size_of_library)
# By whether any virus or none
chisq.test(joined$v_present, joined$clinAD, correct=FALSE)

joined <- within(joined, {coglev.ct <- C(coglev.f, sum)
       print(attributes(coglev.ct))
})

# NOW USING LIBRARY SIZE ADJUSTMENT
# By percent of reads adjusted
t.test(formula=realperclog~clinAD, data=size_of_library)
# By whether any virus or none
chisq.test(joined$v_present, joined$clinAD, correct=FALSE)

joined <- within(joined, {coglev.ct <- C(coglev.f, sum)
       print(attributes(coglev.ct))
})

#model <- lm(factor(coglev.ct) ~ percent + age_death + female, data = joined) 
#model <- lm(factor(coglev.ct) ~ realperclog + age_death + female, data = size_of_library) 
model <- lm(braaksc ~ realperclog + age_death + female, data = size_of_library) 
real_perc_model <- lm(braaksc ~ realperclog, data = size_of_library) 

joined_dementia <- joined[joined$coglev != "MCI",]
# By percent of reads
t.test(formula=percent~coglev, data=joined_dementia)
# By read count
t.test(formula=reads~coglev, data=joined_dementia)
# By whether any virus or none
chisq.test(joined_dementia$v_present, joined_dementia$coglev, correct=FALSE)


# LEt's generate a heatmap
# The goal here is to get the relevant virus rows from each file and make them into a matrix
kreport <- read.table('reports/606_120523.kreport', sep='\t', header=TRUE)
colnames(kreport)[1] <- '%'
# Get the row where the viruses start
v_point <- which(kreport[['taxID']] == 10239)
# The only way to slice in r is with tail so
# get the subtraction point
len <- dim(kreport)[1]
tail_len <- len - v_point
tail(kreport, tail_len+1)
kreport['sample'] <- sample


