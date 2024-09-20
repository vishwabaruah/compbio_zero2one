To calculate TPM (Transcripts Per Million) from gene count or transcript count using `edgeR` in R, follow the script below. It involves several steps: loading the data, normalizing it to library size, and then calculating TPM.

### TPM Calculation Process:
1. **Load count data**: This is your raw gene count or transcript count matrix.
2. **Convert counts to CPM (Counts Per Million)**: CPM is used as an intermediate step.
3. **Calculate gene/transcript lengths**: TPM also accounts for the length of genes/transcripts, so length data is necessary.
4. **Calculate TPM**: Use the CPM values adjusted by gene/transcript lengths.

Here is the complete script:

```r
# Load necessary libraries
if (!requireNamespace("edgeR", quietly = TRUE)) {
    install.packages("edgeR")
}
library(edgeR)

# Load your count matrix (gene counts or transcript counts)
# Ensure that your data is in the correct format: a matrix with genes/transcripts as rows and samples as columns
# Assuming `counts_matrix` is your gene/transcript count matrix
counts_matrix <- read.table("gene_counts.txt", header = TRUE, row.names = 1)

# Assuming you have a vector of gene/transcript lengths in base pairs
# Replace 'gene_lengths.txt' with the appropriate file containing gene lengths
gene_lengths <- read.table("gene_lengths.txt", header = TRUE, row.names = 1)

# Step 1: Create DGEList object for normalization
dge <- DGEList(counts = counts_matrix)

# Step 2: Calculate Counts Per Million (CPM)
cpm_values <- cpm(dge)

# Step 3: Normalize CPM by gene/transcript length
# TPM formula: TPM = CPM / (gene_length / 1e6)
# We need to normalize the CPM by dividing it by the gene length (in kilobases)

# Convert gene lengths from base pairs to kilobases (kb)
gene_lengths_kb <- gene_lengths / 1000

# Step 4: Calculate TPM
# This step divides CPM by gene length (in kb)
tpm_matrix <- sweep(cpm_values, 1, gene_lengths_kb[,1], FUN = "/")

# Step 5: Normalize TPM so that for each sample the sum of TPMs is 1 million
# Calculate the sum of TPMs for each sample
scaling_factors <- colSums(tpm_matrix)

# Divide each sample's TPM by its scaling factor and multiply by 1 million
tpm_matrix <- sweep(tpm_matrix, 2, scaling_factors, FUN = "/") * 1e6

# Step 6: Output the TPM values
write.table(tpm_matrix, file = "tpm_values.txt", sep = "\t", quote = FALSE, col.names = NA)
```

### Explanation:
- **Counts Per Million (CPM)**: The `cpm()` function from `edgeR` converts raw counts to CPM, which adjusts for differences in sequencing depth.
- **Gene lengths**: TPM requires gene lengths (in base pairs). Ensure that your `gene_lengths.txt` file matches the genes in your count matrix.
- **TPM Calculation**: First, divide CPM by gene lengths (in kilobases), then normalize across all genes so that the sum of TPMs in each sample is 1 million.
- **Output**: The result is written to a file called `tpm_values.txt`.

### Assumptions:
- `counts_matrix`: Your gene count matrix (genes as rows, samples as columns).
- `gene_lengths`: A file with gene lengths corresponding to the genes in the count matrix.
