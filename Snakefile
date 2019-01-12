import os, shutil
from glob import glob

#configfile: "config.yaml"
CONDAENV='envs'


## check input
if not 'input' in config:
    raise IOError("need a 'input' argument in the config file")

if ',' in config['input']:
    input_files= config['input'].strip().split(',')
elif os.path.isdir(config['input']):
    input_files = glob(config['input'] +"/*.*")
else:
    input_files =[config['input']]


extensions = set(os.path.splitext(f)[-1] for f in input_files)
assert len(extensions)==1, f"You different fasta extensions {extensions}"

extension = extensions.pop()
del extensions


if not extension in ['.fna','.fasta','.fa_nt']:
    raise IOError(f"Files have {extension} extension which is not a correct fasta nucleotide extension")

for f in input_files:
    assert os.path.exists(f), "{f} doesn't exist"


genome_files = dict( ( os.path.basename(f).replace(extension,''),f)  for f in input_files )
GENOMES = list(genome_files.keys())
print("{} Genomes".format(len(GENOMES)))

def get_genome(wildchards):
    return genome_files[wildchards.genome]


rule all:
    input:
        #expand("annotations/genes/{genome}.faa",genome=GENOMES),
        expand("annotations/genomeproperties/{outfiles}_{genome}",genome=GENOMES, outfiles=['SUMMARY_FILE','TABLE'])

include: "rules/genomeproperties.smk"
