


rule predict_genes:
    input:
        get_genome
    output:
        fna = "annotations/genes/{genome}.fna",
        faa = "annotations/genes/{genome}.faa",
        gff = temp("annotations/genes/{genome}.gff")
    conda:
        "envs/required_packages.yaml"
    log:
        "logs/predict_genes/{genome}.txt"
    threads:
        1
    shell:
        """
        prodigal -i {input} -o {output.gff} -d {output.fna} \
            -a {output.faa} -p meta -f gff 2> >(tee {log})
        """

localrules: remove_asterix
rule remove_asterix:
    input:
        "annotations/genes/{genome}.faa",
    output:
        temp("annotations/genes/{genome}_without_asterix.faa")
    run:
        with open(input[0]) as f, open(output[0],'w') as fo:
            for line in f:
                fo.write(line.replace('*',''))

faa_input=rules.remove_asterix.output

if config['function_predicton']=='interproscan':

    rule interproscan:
        input:
            faa_input
        output:
            "annotations/interproscan/{genome}.tsv"
        params:
            appl='Pfam,TIGRFAM',#,PRINTS,Gene3D,ProSiteProfiles
            extra = config['interproscan_extra']
        log:
            "logs/interproscan/{genome}.txt"
        singularity:
            "docker://continuumio/miniconda3:4.4.10"
        threads:
            config['threads_interproscan']
        shell:
            " interproscan.sh "
            " --input {input} "
            " -appl {params.appl}"
            " --disable-precalc "
            " --cpu {threads} "
            " --formats tsv "
            " --outfile {output} "
            " &> >(tee {log}) "

    inut_genome_properties=rules.interproscan.output




rule genomeproperties:
    input:
        inut_genome_properties
    output:
        expand("annotations/genomeproperties/{outfiles}_{{genome}}",
               outfiles=config['genomeproperties']['outputformats'].values() ),
    params:
        outfiles= ' '.join([f"-outfiles {f}" for f in
                            config['genomeproperties']['outputformats'].keys() ]),
        gpdir='~/Documents/GitHub/genome-properties/flatfiles',
        out_dir= lambda wc, output: os.path.dirname(output[0])
    conda:
        "envs/genomeproperties.yaml"
    threads:
        1
    shell:
        "assign_genome_properties.pl -matches {input} -all  -name {wildcards.genome} -outdir {params.out_dir} "
        "-gpdir {params.gpdir} -gpff genomeProperties.txt {params.outfiles}"

localrules: combine_genome_properties
rule combine_genome_properties:
    input:
        expand("annotations/genomeproperties/SUMMARY_FILE_{genome}",genome=GENOMES),
    output:
        "annotations/Combined_genomeproperties.tsv"
    run:
        import pandas as pd
        G={}
        for i,genome in enumerate(GENOMES):
            G[genome] = pd.read_table(input[i],index_col=[0,1],squeeze=True,header=None)

        G= pd.DataFrame(G).replace({'NO':0,'PARTIAL':0.5,'YES':1})

        G.index.names= ['ID','Description']
        G.to_csv(output[0],sep='\t')





#
#
# usage: java -XX:+UseParallelGC -XX:ParallelGCThreads=2 -XX:+AggressiveOpts
#             -XX:+UseFastAccessorMethods -Xms128M -Xmx2048M -jar
#             interproscan-5.jar
#
#
# Please give us your feedback by sending an email to
#
# interhelp@ebi.ac.uk
#
#  -appl,--applications <ANALYSES>            Optional, comma separated list
#                                             of analyses.  If this option
#                                             is not set, ALL analyses will
#                                             be run.
#  -b,--output-file-base <OUTPUT-FILE-BASE>   Optional, base output filename
#                                             (relative or absolute path).
#                                             Note that this option, the
#                                             --output-dir (-d) option and
#                                             the --outfile (-o) option are
#                                             mutually exclusive.  The
#                                             appropriate file extension for
#                                             the output format(s) will be
#                                             appended automatically. By
#                                             default the input file
#                                             path/name will be used.
#  -cpu,--cpu <CPU>                           Optional, number of cores for
#                                             inteproscan.
#  -d,--output-dir <OUTPUT-DIR>               Optional, output directory.
#                                             Note that this option, the
#                                             --outfile (-o) option and the
#                                             --output-file-base (-b) option
#                                             are mutually exclusive. The
#                                             output filename(s) are the
#                                             same as the input filename,
#                                             with the appropriate file
#                                             extension(s) for the output
#                                             format(s) appended
#                                             automatically .
#  -dp,--disable-precalc                      Optional.  Disables use of the
#                                             precalculated match lookup
#                                             service.  All match
#                                             calculations will be run
#                                             locally.
#  -dra,--disable-residue-annot               Optional, excludes sites from
#                                             the XML, JSON output
#  -f,--formats <OUTPUT-FORMATS>              Optional, case-insensitive,
#                                             comma separated list of output
#                                             formats. Supported formats are
#                                             TSV, XML, JSON, GFF3, HTML and
#                                             SVG. Default for protein
#                                             sequences are TSV, XML and
#                                             GFF3, or for nucleotide
#                                             sequences GFF3 and XML.
#  -goterms,--goterms                         Optional, switch on lookup of
#                                             corresponding Gene Ontology
#                                             annotation (IMPLIES -iprlookup
#                                             option)
#  -help,--help                               Optional, display help
#                                             information
#  -i,--input <INPUT-FILE-PATH>               Optional, path to fasta file
#                                             that should be loaded on
#                                             Master startup. Alternatively,
#                                             in CONVERT mode, the
#                                             InterProScan 5 XML file to
#                                             convert.
#  -iprlookup,--iprlookup                     Also include lookup of
#                                             corresponding InterPro
#                                             annotation in the TSV and GFF3
#                                             output formats.
#  -ms,--minsize <MINIMUM-SIZE>               Optional, minimum nucleotide
#                                             size of ORF to report. Will
#                                             only be considered if n is
#                                             specified as a sequence type.
#                                             Please be aware of the fact
#                                             that if you specify a too
#                                             short value it might be that
#                                             the analysis takes a very long
#                                             time!
#  -o,--outfile <EXPLICIT_OUTPUT_FILENAME>    Optional explicit output file
#                                             name (relative or absolute
#                                             path).  Note that this option,
#                                             the --output-dir (-d) option
#                                             and the --output-file-base
#                                             (-b) option are mutually
#                                             exclusive. If this option is
#                                             given, you MUST specify a
#                                             single output format using the
#                                             -f option.  The output file
#                                             name will not be modified.
#                                             Note that specifying an output
#                                             file name using this option
#                                             OVERWRITES ANY EXISTING FILE.
#  -pa,--pathways                             Optional, switch on lookup of
#                                             corresponding Pathway
#                                             annotation (IMPLIES -iprlookup
#                                             option)
#  -t,--seqtype <SEQUENCE-TYPE>               Optional, the type of the
#                                             input sequences (dna/rna (n)
#                                             or protein (p)).  The default
#                                             sequence type is protein.
#  -T,--tempdir <TEMP-DIR>                    Optional, specify temporary
#                                             file directory (relative or
#                                             absolute path). The default
#                                             location is temp/.
#  -version,--version                         Optional, display version
#                                             number
#  -vtsv,--output-tsv-version                 Optional, includes a TSV
#                                             version file along with any
#                                             TSV output (when TSV output
#                                             requested)
# Copyright © EMBL European Bioinformatics Institute, Hinxton, Cambridge,
# UK. (http://www.ebi.ac.uk) The InterProScan software itself is provided
# under the Apache License, Version 2.0
# (http://www.apache.org/licenses/LICENSE-2.0.html). Third party components
# (e.g. member database binaries and models) are subject to separate
# licensing - please see the individual member database websites for
# details.
#
# Available analyses:
#                       TIGRFAM (15.0) : TIGRFAMs are protein families based on Hidden Markov Models or HMMs
#                          SFLD (4) : SFLDs are protein families based on Hidden Markov Models or HMMs
#                   SUPERFAMILY (1.75) : SUPERFAMILY is a database of structural and functional annotation for all proteins and genomes.
#                        Gene3D (4.2.0) : Structural assignment for whole genes and genomes using the CATH domain structure database
#                         Hamap (2018_03) : High-quality Automated and Manual Annotation of Microbial Proteomes
#                         Coils (2.2.1) : Prediction of Coiled Coil Regions in Proteins
#               ProSiteProfiles (2018_02) : PROSITE consists of documentation entries describing protein domains, families and functional sites as well as associated patterns and profiles to identify them
#                         SMART (7.1) : SMART allows the identification and analysis of domain architectures based on Hidden Markov Models or HMMs
#                           CDD (3.16) : Prediction of CDD domains in Proteins
#                        PRINTS (42.0) : A fingerprint is a group of conserved motifs used to characterise a protein family
#               ProSitePatterns (2018_02) : PROSITE consists of documentation entries describing protein domains, families and functional sites as well as associated patterns and profiles to identify them
#                          Pfam (31.0) : A large collection of protein families, each represented by multiple sequence alignments and hidden Markov models (HMMs)
#                        ProDom (2006.1) : ProDom is a comprehensive set of protein domain families automatically generated from the UniProt Knowledge Database.
#                    MobiDBLite (2.0) : Prediction of disordered domains Regions in Proteins
#                         PIRSF (3.02) : The PIRSF concept is being used as a guiding principle to provide comprehensive and non-overlapping clustering of UniProtKB sequences into a hierarchical order to reflect their evolutionary relationships.
#
# Deactivated analyses:
#                   SignalP_EUK (4.1) : Analysis SignalP_EUK is deactivated, because the resources expected at the following paths do not exist: bin/signalp/4.1/signalp
#                       PANTHER (12.0) : Analysis Panther is deactivated, because the resources expected at the following paths do not exist: data/panther/12.0/panther.hmm, data/panther/12.0/names.tab
#                         TMHMM (2.0c) : Analysis TMHMM is deactivated, because the resources expected at the following paths do not exist: bin/tmhmm/2.0c/decodeanhmm, data/tmhmm/2.0c/TMHMM2.0c.model
#                       Phobius (1.01) : Analysis Phobius is deactivated, because the resources expected at the following paths do not exist: bin/phobius/1.01/phobius.pl
#         SignalP_GRAM_NEGATIVE (4.1) : Analysis SignalP_GRAM_NEGATIVE is deactivated, because the resources expected at the following paths do not exist: bin/signalp/4.1/signalp
#         SignalP_GRAM_POSITIVE (4.1) : Analysis SignalP_GRAM_POSITIVE is deactivated, because the resources expected at the following paths do not exist: bin/signalp/4.1/signalp
