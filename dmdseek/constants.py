"""Fixed GRCh38 coordinates and read filters (0-based half-open intervals)."""
CHROM = 'chrX'
CHRX_LEN_GRCH38 = 156_040_895
DMD_START = 31_119_221          # DMD gene, MANE NM_004006.3 (exon 79 start .. exon 1 end)
DMD_END = 33_339_388
TILE = 200                      # depth bin size (bp)
WGS_FLANK = 1_000_000           # WGS: bins over DMD +/- 1 Mb, flanks used for normalisation
SV_WINDOW = (29_000_000, 35_500_000)   # regional BAM given to Sawfish
MIN_MAPQ = 20
EXCLUDE_FLAGS = 0x4 | 0x100 | 0x200 | 0x400   # unmapped, secondary, QC-fail, duplicate (= mosdepth 1796)
TARGETED_MIN_FRACTION = 0.01    # mapped reads inside DMD / all mapped reads: WGS ~0.0007, capture >0.1
FIRST_LAST_EXONS = {1, 79}
REFERENCE_FORMAT = 2
