# cginess_and_aging

These Perl codes enable to calculate 'distribution shift' and 'net expression change' from RNA-seq data.
(reference will be added soon.)

# How to use
Usage: perl cal_distshift_allpairstesting_rmBoth0_multithreads.pl [read count table]

Example: perl cal_distshift_allpairstesting_rmBoth0_multithreads.pl normcnt_deseq2_GSE121539_mouse_infoAdded.txt

[read count table]: A table of normalized read counts is recommended. DESeq2-normalized counts were used in the study.

To calculate distribution shift, use 'cal_distshift_allpairstesting_rmBoth0_multithreads.pl'.

To calculate net expression change, use 'cal_netexpchange_allpairstesting_rmBoth0_multithreads.pl'.

Their usages are identical.

# Note
1. These codes need 'transpose.pl', which is involved in this project.
2. Input file format : Rows = genes, Columns = samples.
   Each sample should have "y" or "o" as the first word to distinguish groups, which mean 'young' or 'old', respectively.
   (e.g. y.6.M.GSM12345; y=group; 6=age; M=sex) 
3. You will need some Perl packages: 'Statistics::Test::WilcoxonRankSum', 'List::Util' and 'Parallel::ForkManager'.
4. They will use multithreads and you can change the number of threads you want to use (default = 20).
