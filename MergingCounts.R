
### Use R to merge all the counts data into a single file. 

## Have all of your *Counts.txt files in a single folder and set that folder as the working directory
setwd("~/Desktop/CountsData")

rm(final.data)
rm(data)

filenames <- list.files(path = ".", pattern = "*.txt", all.files = FALSE, full.names = FALSE, recursive = FALSE, ignore.case = FALSE)

for (i in filenames)
{
	data <- read.table(i)
	if(!exists("final.data"))
	{
		final.data <- data[,c(1,2,3)]
		colnames(final.data) <- c("gene", "length",i)
		next
	}
	if(exists("final.data"))
	{
	colnames(data)[c(1,3)] <- c("gene",i)
	final.data <- merge(final.data,data[,c(1,3)],by="gene",all.x = TRUE)
}
}

#write.csv(final.data, "combined.data.csv")
write.table (final.data, file="combined.data.txt", row.names=FALSE, quote=FALSE,  sep="\t")

