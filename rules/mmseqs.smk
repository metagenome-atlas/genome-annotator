


rule stockholm2msa:
    input:
        "{database}.stk"
    output:
        'mmseqs/input/msa_trimmed.db'
    conda:
        '../envs/mmseqs.yaml'
    threads: 1
    shell:
        'mmseqs convertmsa {input} {output}'


rule msa2profile:
    input:
        rules.stockholm_to_msa.output
    output:
        "{config['database_folder']}/profiles/{database}"
    conda:
        '../envs/mmseqs.yaml'
    shell:
        'mmseqs msa2profile {input} {output} '
        '--match-mode 1 --msa-type 2 --threads {threads}'

rule make_db:
    input:
        "proteins/{file}.fasta"
    output:
        'queries/{file}.db'
    conda:
        "../envs/mmseqs.yaml"
    threads: 1
    shell:
        'mmseqs createdb {input} {output}'

rule search_mmseqs:
    input:
        profile = "{config['database_folder']}/profiles/{database}",
        fasta = 'queries/{query}.db'
    output:
        db = temp('matches/{database}/{query}.db'),
        index = temp('matches/{database}/{query}.index'),
        tsv = 'matches/{database}/{query}.m8',
    params:
        extra=config.get("mmseqs_search_commands","")
    threads: config['threads']
    conda:
        "../envs/mmseqs.yaml"
    shell:
        """
            mmseqs search {params.extra} --threads {threads} {input.fasta} {input.profile} {output.db} {config[tmpdir]}

            mmseqs convertalis --threads {threads} {input.fasta} {input.profile} {output.db} {output.tsv}
        """
