


rule download_tigrfam:
    output:
        temp("TIGRFAMs_14.0_SEED.tar.gz")
    params:
        tigfam_url= "ftp://ftp.jcvi.org/pub/data/TIGRFAMs/14.0_Release/TIGRFAMs_14.0_SEED.tar.gz"
    shell:
        """
        wget {params.tigfam_url}
        """
