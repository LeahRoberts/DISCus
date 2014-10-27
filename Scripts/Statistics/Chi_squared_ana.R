################################################################################
#                                                                              #
## AUTHOR: Melinda M Ashcroft                                                  #
## AFFILIATION: Beatson Lab | SCMB - University of Queensland                  #
## DATE: April 2014                                                            #
## PURPOSE: Statistical evaluation of association betweetn clades and reads    # 
## orientation               												   #
## TEST: Chi-Square Test                                                       #
#                                                                              #
## THEORY: The chi-square test is used to detect whether there is a            #
## statistically significant difference between the expected frequencies and   #
## observed frequencies in one or multiple categories. For example, do the     #
## counts of your variables that fall in each category differ significantly    #
## from the numbers you expect? Is this difference due to sampling variation   #
## or change, or is it a real difference?                                      #
#                                                                              #
################################################################################

# Set working directory
setwd('~/Dropbox')

# Read in the data
# Use sep='\t' for tab delimited or sep=',' for comma delimited data
reads = read.csv(file.choose(), header=T, sep='\t')

# Attach variable (column) names
attach(reads)
# See the first few observations
head(reads)
# Show the names of each column
names(reads)

# View a side by side boxplot for reads in on vs off orientations
# for each clade or sub clade
op = par(mfrow=c(1,4)) 
boxplot(reads$OFF_1 ~ reads$Clade, las = 2, main = "Off orientation reads 1")
boxplot(reads$OFF_2 ~ reads$Clade, las = 2, main = "Off orientation reads 2")
boxplot(reads$ON_1 ~ reads$Clade, las = 2, main = "On orientation reads 1")
boxplot(reads$ON_2 ~ reads$Clade, las = 2, main = "On orientation reads 2")
par(op)

# Chi square tables 3x2 contingency table using 
# Off_1 vs On_1 for all major clades
C <- matrix(c(660,1008,5515,28,248,170), nrow=3)
C

# Chi square summaries
# Perform Chi-square test on data contained in variable C
chisq.test(C)          
# The ratio of the difference between the observed count and the expected
# count to the standard deviation of the expected count is called the 
# standardized residual
chisq.test(C)$stdres   # standardized residuals
chisq.test(C)$expected # expected numbers

# Subsetting data - multiple Chi-square tests - subgroups of each clade
# Currently only using off_1 vs on_1
# Does not take into account multiple testing
# A to B1
A_B1 <- matrix(c(1,398, 7, 345), nrow=2) # insert the data
A_B1 # view the data
chisq.test(A_B1) # perform Chi-square test on the data

# A to B2
A_B2 <- matrix(c(617,155, 17, 97), nrow=2) # insert the data
A_B2 # view the data
chisq.test(A_B2) # perform Chi-square test on the data

# A to C1
A_C1 <- matrix(c(617,1513, 17, 63), nrow=2) # insert the data
A_C1 # view the data
chisq.test(A_C1) # perform Chi-square test on the data

# A to C2
A_C2 <- matrix(c(617,4132, 17, 126), nrow=2) # insert the data
A_C2 # view the data
chisq.test(A_C2) # perform Chi-square test on the data

# B1 to B2
B1_B2 <- matrix(c(745,155, 62, 97), nrow=2) # insert the data
B1_B2 # view the data
chisq.test(B1_B2) # perform Chi-square test on the data

# B1 to C1
B1_C1 <- matrix(c(745,1513, 62, 63), nrow=2) # insert the data
B1_C1 # view the data
chisq.test(B1_C1) # perform Chi-square test on the data

# B1 to C2
B1_C2 <- matrix(c(745,4132, 62, 126), nrow=2) # insert the data
B1_C2 # view the data
chisq.test(B1_C2) # perform Chi-square test on the data

# B2 to C1
B2_C1 <- matrix(c(155,1513, 97, 63), nrow=2) # insert the data
B2_C1 # view the data
chisq.test(B2_C1) # perform Chi-square test on the data

# B2 to C2
B2_C2 <- matrix(c(155,4132, 97, 126), nrow=2) # insert the data
B2_C2 # view the data
chisq.test(B2_C2) # perform Chi-square test on the data

# C1 to C2
C1_C2 <- matrix(c(1513,4132, 63, 126), nrow=2) # insert the data
C1_C2 # view the data
chisq.test(C1_C2) # perform Chi-square test on the data
