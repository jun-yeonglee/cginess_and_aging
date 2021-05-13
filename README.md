# cginess_and_aging

These Perl codes enable to calculate 'distribution shift' and 'net expression change' from RNA-seq data.
(reference will be added soon.)

# cal_distshift_allpairstesting_rmBoth0_multithreads.pl
For distribution shift.

# cal_netexpchange_allpairstesting_rmBoth0_multithreads.pl
For net expression change.

# Note
1. These codes need 'transpose.pl', which is involved in this project.
2. Input file format : Rows = genes, Columns = samples.
   Each sample should have "y" or "o" as the first word to distinguish groups, which mean 'young' or 'old', respectively, to distinguish their group.
   (e.g. y.6.M.GSM12345; y=group; 6=age; M=sex) 
3. You will need to load ‘Statistics::Test::WilcoxonRankSum’, 'List::Util' and 'Parallel::ForkManager' packages.
   It will use multithreads and you can change the number of threads you want to use (default = 20).
