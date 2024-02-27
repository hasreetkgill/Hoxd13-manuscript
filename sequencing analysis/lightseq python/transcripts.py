from collections import defaultdict 
import glob
import pandas as pd
import pprint
import pysam
import statistics

DEDUP_FILES = sorted(glob.glob('outFiles_trimmed_merged/*_Dedup.bam'))
GFF_FILE = 'Gallus_gallus.GRCg6a.105.gtf'
OUT_FILE_FREQ = 'Hasreet_TranscriptFrequencies_MidHindHoxD13.csv'
OUT_FILE = 'Hasreet_ReorderedLightSeq_MidHindHoxD13.csv'

def readGFF():
  genes = {}
  
  print('Reading GFF file...')
  with open(GFF_FILE) as gffHandle:
    for line in gffHandle:
      if line[0] != '#' and line.split('\t')[2] == 'gene':
        genes[line.split('gene_id "')[1].split('";')[0]] = line
  
  return genes

def main():
  filesToAnalyze = DEDUP_FILES

  genes = readGFF()
  allGenes = []
  allConditionBarcodes = []
  frequencies = defaultdict(lambda: defaultdict(int)) # default to 0 count

  # Create dictionary of gene frequencies
  for geneFile in filesToAnalyze:  
    print('Reading %s...' % geneFile)

    condition = geneFile.split('/')[-1].split('_Dedup')[0]
    samFile = pysam.AlignmentFile(geneFile, "rb")
    samIter = samFile.fetch(until_eof=True)

    for read in samIter:
      barcodeSeq = read.query_name.split('_')[1]
      geneName = read.get_tag('XT')

      condBar = '%s/%s' % (condition, barcodeSeq)

      frequencies[geneName][condBar] += 1

      if not condBar in allConditionBarcodes:
        allConditionBarcodes.append(condBar)

  # Now, create CSV
  allData = defaultdict(list)

  allData['Gene'] = frequencies.keys()

  print('Compiling data...')
  for gene in allData['Gene']:
    try:
      allData['Gene name'].append(genes[gene].split('gene_name "')[1] \
                                  .split('";')[0])
    except:
      allData['Gene name'].append(genes[gene].split('gene_id "')[1] \
                                  .split('";')[0])
    for condBar in allConditionBarcodes:
      allData[condBar].append(frequencies[gene][condBar])

  lightFreqs = pd.DataFrame.from_dict(allData)
  lightFreqs.to_csv(OUT_FILE_FREQ, index=False)

  # Now, re-order data according to layer, replicate information
  del lightFreqs['Gene']


  lightFreqs.rename(columns={'Gene name': 'Gene',
                             'HindA/AGGGTA': 'InnerMes1',
                             'HindB/AGGGTA': 'InnerMes2',
                             'HindC/AGGGTA': 'InnerMes3',
                             'HindD/AGGGTA': 'InnerMes4',
                             'HindA/GTTAGG': 'MuscMuco1',
                             'HindB/GTTAGG': 'MuscMuco2',
                             'HindC/GTTAGG': 'MuscMuco3',
                             'HindD/GTTAGG': 'MuscMuco4',
   			     'HindA/TATGGA': 'InnerCirc1',
                             'HindB/TATGGA': 'InnerCirc2',
                             'HindC/TATGGA': 'InnerCirc3',
                             'HindD/TATGGA': 'InnerCirc4',
                             'HoxA/AGGGTA': 'InnerMes1',
                             'HoxB/AGGGTA': 'InnerMes2',
                             'HoxC/AGGGTA': 'InnerMes3',
                             'HoxD/AGGGTA': 'InnerMes4',
 			     'HoxA/GTTAGG': 'MuscMuco1',
                             'HoxB/GTTAGG': 'MuscMuco2',
                             'HoxC/GTTAGG': 'MuscMuco3',
                             'HoxD/GTTAGG': 'MuscMuco4',
			     'HoxA/TATGGA': 'InnerCirc1',
                             'HoxB/TATGGA': 'InnerCirc2',
                             'HoxC/TATGGA': 'InnerCirc3',
                             'HoxD/TATGGA': 'InnerCirc4',
                             'MidA/AGGGTA': 'InnerMes1',
                             'MidB/AGGGTA': 'InnerMes2',
                             'MidC/AGGGTA': 'InnerMes3',
                             'MidD/AGGGTA': 'InnerMes4',
                             'MidA/GTTAGG': 'MuscMuco1',
                             'MidB/GTTAGG': 'MuscMuco2',
                             'MidC/GTTAGG': 'MuscMuco3',
                             'MidD/GTTAGG': 'MuscMuco4',
                             'MidA/TATGGA': 'InnerCirc1',
                             'MidB/TATGGA': 'InnerCirc2',
                             'MidC/TATGGA': 'InnerCirc3',
    		             'MidD/TATGGA': 'InnerCirc4',
                             },
                    inplace=True)

  lightFreqs = lightFreqs[['Gene', 'InnerMes1', 'InnerMes2', 'InnerMes3', 'InnerMes4', 'MuscMuco1', 'MuscMuco2', 
                           'MuscMuco3', 'MuscMuco4', 'InnerCirc1', 'InnerCirc2', 'InnerCirc3', 'InnerCirc4','InnerMes1',
			   'InnerMes2','InnerMes3','InnerMes4', 'MuscMuco1', 'MuscMuco2','MuscMuco3', 'MuscMuco4', 
                           'InnerCirc1', 'InnerCirc2', 'InnerCirc3', 'InnerCirc4', 'InnerMes1', 'InnerMes2', 'InnerMes3', 
                           'InnerMes4', 'MuscMuco1', 'MuscMuco2','MuscMuco3', 'MuscMuco4', 'InnerCirc1', 'InnerCirc2', 
                           'InnerCirc3', 'InnerCirc4']]

  lightFreqs.set_index('Gene', inplace=True)
  # Add together counts with duplicate indices
  lightFreqs = lightFreqs.groupby(lightFreqs.index).sum()

  lightFreqs.to_csv(OUT_FILE)

if __name__ == '__main__':
  main()
