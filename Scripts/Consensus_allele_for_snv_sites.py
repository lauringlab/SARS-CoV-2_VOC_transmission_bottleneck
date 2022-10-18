

# ====== Import Things ======
import os
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


# ================================ Main =========================================

consensus = "Results/consensus_for_sample_pairs_aligned.fasta"
sites =  (29708, 29440, 6352, 6470, 12094, 5759, 11020, 25792, 27416, 29467, 7113, 1960, 17392, 23123, 14484, 1684, 26256, 10852, 13665, 24999, 28201, 28079, 14448, 18907, 21114, 4708, 17550, 29769, 7279, 8825, 22104, 5770, 14599, 1623, 18788, 23902, 24793, 24794, 27688, 9757, 8964, 15172, 3287, 17746, 27671, 1347, 18877, 500, 22997, 29751, 291, 12076, 13743, 27213, 11633, 27230, 3611, 13038, 15042, 4905, 9070, 10889, 27522, 196, 10798, 15346, 14094, 18417, 2062, 9737, 28809, 28657, 6372, 29224, 11750, 29185, 805, 24724, 420, 25452, 1878, 884, 10802, 4230, 27244, 15654, 16466, 10135, 12379, 2147, 21595, 25708, 29235, 29301, 29774, 3283, 25354, 27619, 29022, 19476)
f = open("Results/consensus_allele_forsnv_sites.txt", "a")
header = "sample Allele POS"
f.write (header + '\n')

for record in SeqIO.parse(consensus, "fasta"):
    for site in sites:
         bp=(record.id, record.seq[site-1], site)
         f.write (' '.join(str(s) for s in bp) + '\n')
f.close ()    
