from sys import argv
import GEOparse

if len(argv)<3:
    print("""Usage: python 00_download_GEO.py GEO_Accession E-mail
Please add Red-C GEO accession and your e-mail address for SRA download as arguments.
Requirements: 
  GEOParse https://GEOparse.readthedocs.org
  gzip 
  sra-tools
""")
    exit()

geo_id = argv[1]
email = argv[2]

fastq_dump_options={'split-files': None,
    'read-filter': 'pass',
    'dumpbase': None,
    'gzip': None}

gse = GEOparse.get_GEO(geo=geo_id, destdir="./")
gsms = gse.gsms
downloaded_paths = gse.download_SRA(email,
                                    filetype='fastq')
                                    #fastq_dump_options=fastq_dump_options,
                                    #nproc=nproc,
                                    #silent=True)