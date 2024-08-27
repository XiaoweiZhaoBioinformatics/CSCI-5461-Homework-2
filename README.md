# CSCI 5461 Homework 2

The cell lines were derived from a set of breast tumors or normal tissue. Each of the cell
lines was exposed independently to five drugs, and the GI50 was measured. GI50 represents the
drug concentration at which growth is inhibited by 50%. The GI50 values have been – log10
transformed so higher values reflect greater sensitivity. For simplicity, in this assignment we
have binarized the drug sensitivity values such that “1” corresponds to a drug-sensitive cell line
and “0” corresponds to a drug-resistant cell line.

**The goal is to build a kNN classifier that predicts which cell lines are sensitive to each drug.**

## Data files:

The drug sensitivity and gene expression data have been merged into a single text file,
called DREAM_data.txt. This is a tab-delimited file with the cell line IDs across the first
row, and the drug sensitivity and gene expression data in the rows that follow. Please read the
following important notes:
- Drug sensitivity data (2nd to 6th rows of the file): the first column is the drug name and the rest
of the columns are the drug sensitivity values (1=drug sensitive, 0=not sensitive). “NA” values
correspond to missing values, which are instances where the drug sensitivity was not measured
well—these should not be included in either the sensitive or the resistant groups.
- Gene expression data (7th row to the end): the first column is the gene name and the rest of the
columns are gene-level summaries across breast cancer cell lines. Gene-level summaries of
expression are already quantile normalized and log2-transformed.

## Questions:
**1. k-NN implementation:**

Implement a kNN classifier, using the Pearson
correlation coefficient as the similarity metric between cell lines’ expression profiles. You
should create a classifier for each drug that takes a cell line expression profile as input and
produces a score that predicts whether it is sensitive or resistant to the given drug. Your
classification score for each cell line should be the fraction of the k-nearest neighbors that are
sensitive to the corresponding drug.

**2. k-NN performance evaluation:**

Set k=5 for all of the evaluations in this
problem. Apply leave-one-out cross-validation (LOOCV) to measure the performance of
your classifier. For each of the 5 drugs, plot an ROC curve with sensitivity on the y-axis and
(1-specificity) on the x-axis. Use the number of drug-sensitive neighbors (among the k total)
for each cell line to rank the predictions in order to draw the ROC curve. For each ROC
curve, be sure to plot the performance expected by a random classifier. Answer the following
questions:
- (a) Does your classifier work better than a random classifier? For which drugs? Refer to
specific evidence from your analysis to justify your answer.
- (b) For the drug on which your approach appears to work the best, how good is the
classification performance? Pick a point on the ROC curve, and report the relevant
metrics (number of true positives, false positives).

**3. Exploration of parameters affecting kNN performance:** 

For this problem, we
will explore how the performance is influenced by the key parameter in the kNN classifier, k,
and the scoring function.
- (a) Rerun your classification results with k=3,5,7. For each drug, create a single figure, but
plot the ROC curves for all values of k on the same curve. Again, sort the predicted cell
lines by the number of their nearest neighbors that are drug-sensitive. Discuss the results
of your analysis. Does the choice of k affect the performance of the classifier?
- (b) Set k=5. Instead of scoring each cell line with the number of its nearest neighbors that
are drug-sensitive, use the following weighted score:

<img width="685" alt="Screenshot 2024-08-27 at 9 36 41 AM" src="https://github.com/user-attachments/assets/05922fa5-0b2f-45ba-9661-549c4270e748">

- Where yi are the nearest neighbors and 𝑠𝑖𝑔𝑛(𝑦𝑖) = +1 for drug-sensitive cell lines and
𝑠𝑖𝑔𝑛(𝑦𝑖) = −1 for drug-resistant cell lines, and 𝑃𝐶𝐶(𝑥, 𝑦𝑖) is the Pearson correlation
coefficient. Sort the predictions of drug-sensitivity by this score. Plot new ROC curves
for each drug comparing the approach from Problem (2) with this approach.
