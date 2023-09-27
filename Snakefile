import os, shutil
from glob import glob
import yaml
#configfile: "config.yaml"
CONDAENV='envs'

default_config = yaml.load(open(os.path.join(os.path.dirname(workflow.snakefile),'config.yaml')))
update_config(default_config, config)
config = default_config

#HACK:

genprop_flat_file='/Users/silas/Documents/GitHub/genome-properties/flatfiles/genomeProperties.txt'



## check input
if not 'input' in config:
    raise IOError("need a 'input' argument in the config file")

if ',' in config['input']:
    input_files= config['input'].strip().split(',')
elif os.path.isdir(config['input']):
    input_files = glob(config['input'] +"/*.*")
else:
    input_files =[config['input']]


for f in input_files:
    assert os.path.exists(f), f"{f} doesn't exist"


extensions = set(os.path.splitext(f)[-1] for f in input_files)
assert len(extensions)==1, f"You different fasta extensions {extensions}"

extension = extensions.pop()
del extensions


if not extension in ['.fna','.fasta','.fa_nt']:
    raise IOError(f"Files have {extension} extension which is not a correct fasta nucleotide extension")



genome_files = dict( ( os.path.basename(f).replace(extension,''),f)  for f in input_files )
GENOMES = list(genome_files.keys())
print("{} Genomes".format(len(GENOMES)))

def get_genome(wildchards):
    return genome_files[wildchards.genome]


rule all:
    input:
        "annotations/GenomeProperties.tsv"

include: "rules/genomeproperties.smk"
