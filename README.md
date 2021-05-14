# cginess_and_aging

These Perl codes enable to calculate 'distribution shift' and 'net expression change' from RNA-seq data.
(reference will be added soon.)


# How to use
Usage: perl cal_distshift_allpairstesting_rmBoth0_multithreads.pl [#threads] [file name]

Example: perl cal_distshift_allpairstesting_rmBoth0_multithreads.pl 16 normcnt_deseq2

[#threads]: These codes run with multithreads. you can choose the number of threads you want to use.

[file name]: Table(s) of read counts. Normalized read counts are recommended. DESeq2-normalized counts were used in the study.

You can specify a file name for a single run, or can use common file name that put together your target files (multithreading).  

To calculate distribution shift, use 'cal_distshift_allpairstesting_rmBoth0_multithreads.pl'.

To calculate net expression change, use 'cal_netexpchange_allpairstesting_rmBoth0_multithreads.pl'.


# Note
1. These codes need 'transpose.pl', which is involved in this project.
2. Input file format : Rows = genes, Columns = samples.
   Each sample should have "y" or "o" as the first word to distinguish groups, which mean 'young' or 'old', respectively.
   (e.g. y.6.M.GSM12345; y=group; 6=age; M=sex) 
3. You will need some Perl packages: 'Statistics::Test::WilcoxonRankSum', 'List::Util' and 'Parallel::ForkManager'.
