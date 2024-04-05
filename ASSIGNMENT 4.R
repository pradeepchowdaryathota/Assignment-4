install.packages("seqinr")
install.packages("stringr")
library(stringr)
library(seqinr)
# Load the file into a variable"
sequence <- read.fasta(file = "chr1_GL383518v1_alt.fa", "r", seqtype="DNA")
sequence <- as.character(sequence[[1]])

# Join lines to form the sequence string
sequence <- paste(sequence, collapse = "")

# Now, 'sequence' contains the content of the file as a string
#Q1.
# Print the 10th letter of the sequence
print(paste("10th letter of the sequence:", substr(sequence,10,10)))

# Print the 758th letter of the sequence
print(paste("758th letter of the sequence:",substr(sequence,758,758)))

#Q2
# reverse the sequence
reverse_sequence <- rev(strsplit(sequence, "")[[1]])

#complement the data
complement_dict <- list("A"="T", "T"="A", "G"="C", "C"="G", "a"="t", "t"="a", "g"="c", "c"="g")

#join the reversed sequence with complemented sequence
complemented_sequence <- paste(sapply(reverse_sequence, function(base) complement_dict[(base)]), collapse = "")

# Convert the complemented sequence to a single string
complemented_sequence <- paste(complemented_sequence, collapse = "")

#Print the 79th letter of this sequence.
print(paste("79th letter of the sequence:", substr(complemented_sequence,79,79)))

#Print the 500th through the 800th letters of this sequence.
print(paste("500 through 800letters of the sequence:", substr(complemented_sequence,500,800)))

#Q.3 Create a list that contains the number of times each letter appears in the downloaded sequence, as a function of which kilobase of the sequence you are looking at.

# Define the window size (in bases)
window_size <- 1000

# Initialize a list to store the counts for each kilobase
kilobase_counts <- list()

# Iterate over the sequence in kilobase windows
for (i in seq(1, nchar(sequence), by = window_size)) {
  start <- i
  end <- min(i + window_size - 1, nchar(sequence))
  
  window <- substr(sequence, start, end)
  
  # Count the occurrences of each letter in the window
  window_counts <- c(A = 0, T = 0, G = 0, C = 0, a = 0, t = 0, g = 0, c = 0)
  for (base in strsplit(window, "")[[1]]) {
    window_counts[base] <- window_counts[base] + 1
  }
  
  # Add the counts to the list
  kilobase_counts[[paste0("Kilobase", i/window_size)]] <- window_counts
}
print(kilobase_counts)

#Q4.
#A. Create a data frame with 4 columns, containing the number of times each nucleotide (A,C,G,T) is contained in the first 1000 base pairs.
# Extract counts for the first kilobase and convert to numeric
first_kilobase_counts <- as.numeric(kilobase_counts[[1]])

# Extract the last four counts
last_four_counts <- tail(first_kilobase_counts, 4)

# Create a data frame with counts for each nucleotide
NUCLEOTIDE_counts <- data.frame(
  Nucleotide = c("A", "T", "C", "G"),
  Count = last_four_counts
)

print(NUCLEOTIDE_counts)

#4B.
# Initialize an empty list to store the nucleotide counts for each kilobase
list1 <- list()

# Iterate over each kilobase and extract the nucleotide counts
for (kilobase in names(kilobase_counts)) {
  nucleotide_counts <- kilobase_counts[[kilobase]]
  last_four_counts <- tail(nucleotide_counts, 4)
  list1[[kilobase]] <- unlist(last_four_counts)
}
# Print the nucleotides in each kilobase
print("nucleotides in each kilobase:")
print(list1)

#4c.
# Combine the individual lists into a single data frame
combined_df <- do.call(rbind, list1)

# Print the combined data frame
print(combined_df)

#4d.
# Calculate the sum of each row
row_sums <- rowSums(combined_df)

# Print the row sums
print(row_sums)

#4e.
# 1. Expected sum for the each list is 1000.
# 2.yes, sum of the list of every kilobase is 1000 except the last kilobase is with only the sum of 439.
# 3. The difference is because the we assigned a value of 1000 for kilobase and the in the last kilobase we only had the 439 seuences where the sequence of teh chromosome ended.so,it has different sum than expected.


















