import glob
import os

IN_FILES = sorted(glob.glob('outFiles_trimmed_merged/*Merged_R1.fastq.gz'))
OUT_DIR = 'outFiles_trimmed_merged/'

CHICKEN_INDEX_PATH = '../Gallus_Indexed'
STAR_PATH = '../STAR-2.7.3a/source/STAR'

def mapChickenSeq(trimmedR1File, chickenR1Prefix):
    # Map to Chicken
    print('  Mapping to Chicken genome...')
    os.system(('%s --runThreadN 12 --quantMode TranscriptomeSAM GeneCounts ' \
               'GeneCounts --genomeDir %s --readFilesIn %s ' \
               '--outFileNamePrefix %s --readFilesCommand zcat ' \
               '--outFilterMultimapNmax 1 ' \
               '--outSAMtype BAM SortedByCoordinate') % \
                 (STAR_PATH, CHICKEN_INDEX_PATH, trimmedR1File, chickenR1Prefix))

def main():
  parseFiles = []
  for inFile in IN_FILES:
    filePrefix = '%s%s' % (OUT_DIR, inFile.split('/')[-1].split('_')[0])
    #trimmedR1File = '%s_Merged_R.fastq.gz' % filePrefix
    chickenR1Prefix = '%s_R1_Chicken' % filePrefix
    parseFiles.append((inFile, chickenR1Prefix))

  for pFiles in parseFiles:
    mapChickenSeq(*pFiles)

if __name__ == '__main__':
  main()
