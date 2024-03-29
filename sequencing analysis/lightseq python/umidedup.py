import glob
import os
import pprint
import pysam
import statistics
from collections import defaultdict

GENE_FILES = sorted(glob.glob('outFiles_trimmed_merged/*Aligned.sortedByCoord.out.bam'))
GFF_FILE = 'Gallus_gallus.GRCg6a.105.gtf'

OUT_DIR = 'outFiles_trimmed_merged/'

#Two mitochondrial reads have to be ignored because they soak up over 1 million reads and slow down the UMI dedup
GENES_TO_FILTER = ['ENSGALG00000036956', 'ENSGALG00000043598']

def createSmallerGeneFile(filePrefix, geneFile):
  outGeneFile = '%s.bam' % filePrefix

  gF = pysam.AlignmentFile(geneFile)
  outGF = pysam.AlignmentFile(outGeneFile, 'w', template=gF)
  
  filteredCounts = defaultdict(int)

  for read in gF.fetch(until_eof=True):
    if read.has_tag('XT'):
      xtTag = read.get_tag('XT')

      if xtTag in GENES_TO_FILTER:
        filteredCounts[xtTag] += 1
      else:
        outGF.write(read)

  gF.close()
  outGF.close()

  print('Filtered out %s' % ', '.join(['%d of %s' % \
                                       (filteredCounts[gene], gene) \
                                       for gene in GENES_TO_FILTER]))

  return outGeneFile

def deduplicate(filePrefix, geneFile):
  assignedGeneFile = '%s_GeneAssigned' % filePrefix
  countedGeneFile = '%s.featureCounts.bam' % geneFile  
  sortedGeneFile = '%s_Sorted.bam' % filePrefix
  dedupGeneFile = '%s_Dedup.bam' % filePrefix
  dedupLogFile = '%s_Dedup.log' % filePrefix

  # Use feature-counts to count genes
  print('  Counting genes with featureCounts...')
  os.system('featureCounts -a %s -o %s -R BAM %s' % \
            (GFF_FILE, assignedGeneFile, geneFile))

  # sort and Index the featureCounts Files

  os.system('samtools sort %s -o %s' % (countedGeneFile,
                                        countedGeneFile))
  os.system('samtools index %s' % countedGeneFile)

  # Now, filter out reads from the two problem loci
  filteredCountsFile = createSmallerGeneFile(filePrefix, countedGeneFile)

  # De-duplicate
  print('  Deduping...')

  # Filter, sort, index the bam files then deduplicate
  # NOTE: the filtering is not random
  os.system('samtools sort %s -o %s' % (filteredCountsFile,
                                        sortedGeneFile))

  os.system('samtools index %s' % sortedGeneFile)

  os.system(('umi_tools dedup --per-gene --gene-tag=XT ' \
             '--assigned-status-tag=XS -I %s -S %s -L %s') \
             % (sortedGeneFile, dedupGeneFile,
                dedupLogFile))

def main():  
  # Deduplicate each file individually
  for geneFile in GENE_FILES:
    filePrefix = '%s%s' % (OUT_DIR, geneFile.split('/')[1].split('_R1_')[0])
    deduplicate(filePrefix, geneFile)

if __name__ == '__main__':
  main()
