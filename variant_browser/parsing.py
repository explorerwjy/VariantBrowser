
"""
Utils for reading flat files that are loaded into database
"""
import re
import traceback
from utils import *
import copy

POPS = {
    'AFR': 'African',
    'AMR': 'Latino',
    'EAS': 'East Asian',
    'FIN': 'European (Finnish)',
    'NFE': 'European (Non-Finnish)',
    'SAS': 'South Asian',
    'OTH': 'Other'
}

#======================================================================================
# This Func is Alpha version, read variants from tsv file
# Next version should be read variants from vcf file
#======================================================================================
def get_variants_from_sites_vcf(sites_vcf):
    sites_vcf = sites_vcf[0]
    print sites_vcf
    sites_vcf = open(sites_vcf, 'rb')

    for line in sites_vcf:
        try:

            line = line.strip('\n')
            line = line.strip('\r')
            # Header line, get index of fields we interested
            if line.startswith('Chr'):
                print line
                headers = line.split('\t')
                Chr = headers.index('Chr')
                Start = headers.index('Start')
                End = headers.index('End')
                Ref = headers.index('Ref')
                Alt = headers.index('Alt')
                Func_refgene = headers.index('Func.refgene')
                Gene_refgene = headers.index('Gene.refgene')
                ExonicFunc_refgene = headers.index('ExonicFunc.refgene')
                AAChange_refgene = headers.index('AAChange.refgene')
                ExAC_Freq = headers.index('ExAC_Freq')
                dbSNP = headers.index('dbSNP')
                RadialSVM_pred = headers.index('RadialSVM_pred')
                CADD_phred = headers.index('CADD_phred')
                QUAL = headers.index('QUAL')
                Filter = headers.index('Filter')
                Disease = headers.index('Disease')

            # Get variants
            else:
                print line
                info = line.split('\t')

                alt_alleles = info[Alt].split(',')
                for i, alt_allele in enumerate(alt_alleles):
                    variant = {}

                    pos, ref, alt = get_minimal_representation(info[Start], info[Ref], alt_allele)

                    variant['chrom'] = info[Chr]
                    variant['ref'] = ref
                    variant['alt'] = alt
                    variant['pos'] = pos
                    variant['xpos'] = get_xpos(variant['chrom'], variant['pos'])
                    variant['xstart'] = variant['xpos']
                    variant['xstop'] = variant['xpos'] + len(variant['alt']) - len(variant['ref'])
                    variant['variant_id'] = '{}-{}-{}-{}'.format(variant['chrom'], variant['pos'], variant['ref'], variant['alt'])
                    variant['End'] = info[End]
                    variant['Func_refgene'] = info[Func_refgene]
                    variant['genes'] = info[Gene_refgene]
                    variant['ExonicFunc_refgene'] = info[ExonicFunc_refgene]
                    variant['AAChange_refgene'] = info[AAChange_refgene]
                    variant['ExAC_Freq'] = info[ExAC_Freq]
                    variant['rsid'] = info[dbSNP]
                    variant['RadialSVM_pred'] = info[RadialSVM_pred]
                    variant['CADD_phred'] = info[CADD_phred]
                    variant['site_quality'] = info[QUAL]
                    variant['filter'] = info[Filter]
                    variant['Disease'] = info[Disease].split('=')[1]

                    yield variant
        except Exception:
            print("Error parsing vcf line: " + line)
            traceback.print_exc()
            break


def get_canonical_transcripts(canonical_transcript_file):
    for line in canonical_transcript_file:
        gene, transcript = line.strip().split()
        yield gene, transcript


def get_omim_associations(omim_file):
    for line in omim_file:
        fields = line.strip().split('\t')
        if len(fields) == 4:
            yield fields
        else:
            yield None

def get_genes_from_gencode_gtf(gtf_file):
    """
    Parse gencode GTF file;
    Returns iter of gene dicts
    """
    for line in gtf_file:
        if line.startswith('#'):
            continue
        fields = line.strip('\n').split('\t')

        if fields[2] != 'gene':
            continue

        chrom = fields[0]
        chrs = map(str,[i for i in xrange(1,23)])
        chrs.extend(['X','Y'])
        if chrom not in chrs:
            continue

        start = int(fields[3]) + 1  # bed files are 0-indexed
        stop = int(fields[4]) + 1
        try:
            info = dict((x.strip().split()[0], ' '.join(x.strip().split()[1:])) for x in fields[8].split(';') if x != '')
            info = {k: v.strip('"') for k, v in info.items()}
        except ValueError:
            print fields[8]
            print fields[8].split(';')
            exit()
        gene_id = info['gene_id'].split('.')[0]

        gene = {
            'gene_id': gene_id,
            'gene_name': info['gene_name'],
            'gene_name_upper': info['gene_name'].upper(),
            'chrom': chrom,
            'start': start,
            'stop': stop,
            'strand': fields[6],
            'xstart': get_xpos(chrom, start),
            'xstop': get_xpos(chrom, stop),
        }
        yield gene


def get_transcripts_from_gencode_gtf(gtf_file):
    """
    Parse gencode GTF file;
    Returns iter of transcript dicts
    """
    for line in gtf_file:
        if line.startswith('#'):
            continue
        fields = line.strip('\n').split('\t')

        if fields[2] != 'transcript':
            continue

        chrom = fields[0]
        chrs = map(str,[i for i in xrange(1,23)])
        chrs.extend(['X','Y'])
        if chrom not in chrs:
            continue

        start = int(fields[3]) + 1  # bed files are 0-indexed
        stop = int(fields[4]) + 1
        try:
            info = dict((x.strip().split()[0], ' '.join(x.strip().split()[1:])) for x in fields[8].split(';') if x != '')
            info = {k: v.strip('"') for k, v in info.items()}
        except ValueError:
            print fields[8]
            print fields[8].split(';')
            exit()
        transcript_id = info['transcript_id'].split('.')[0]
        gene_id = info['gene_id'].split('.')[0]

        gene = {
            'transcript_id': transcript_id,
            'gene_id': gene_id,
            'chrom': chrom,
            'start': start,
            'stop': stop,
            'strand': fields[6],
            'xstart': get_xpos(chrom, start),
            'xstop': get_xpos(chrom, stop),
        }
        yield gene


def get_exons_from_gencode_gtf(gtf_file):
    """
    Parse gencode GTF file;
    Returns iter of transcript dicts
    """
    for line in gtf_file:
        if line.startswith('#'):
            continue
        fields = line.strip('\n').split('\t')

        if fields[2] not in ['exon', 'CDS', 'UTR']:
            continue

        chrom = fields[0]
        chrs = map(str,[i for i in xrange(1,23)])
        chrs.extend(['X','Y'])
        if chrom not in chrs:
            continue


        feature_type = fields[2]
        start = int(fields[3]) + 1  # bed files are 0-indexed
        stop = int(fields[4]) + 1
        try:
            info = dict((x.strip().split()[0], ' '.join(x.strip().split()[1:])) for x in fields[8].split(';') if x != '')
            info = {k: v.strip('"') for k, v in info.items()}
        except ValueError:
            print fields[8]
            print fields[8].split(';')
            exit()
        transcript_id = info['transcript_id'].split('.')[0]
        gene_id = info['gene_id'].split('.')[0]

        exon = {
            'feature_type': feature_type,
            'transcript_id': transcript_id,
            'gene_id': gene_id,
            'chrom': chrom,
            'start': start,
            'stop': stop,
            'strand': fields[6],
            'xstart': get_xpos(chrom, start),
            'xstop': get_xpos(chrom, stop),
        }
        yield exon



def get_dbnsfp_info(dbnsfp_file):
    """
    Parse dbNSFP_gene file;
    Returns iter of transcript dicts
    """
    header = dbnsfp_file.next().split('\t')
    fields = dict(zip(header, range(len(header))))
    for line in dbnsfp_file:
        line = line.split('\t')
        other_names = line[fields["Gene_old_names"]].split(';') if line[fields["Gene_old_names"]] != '.' else []
        if line[fields["Gene_other_names"]] != '.':
            other_names.extend(line[fields["Gene_other_names"]].split(';'))
        gene_info = {
            'gene_name': line[fields["Gene_name"]],
            'ensembl_gene': line[fields["Ensembl_gene"]],
            'gene_full_name': line[fields["Gene_full_name"]],
            'gene_other_names': other_names
        }
        yield gene_info

def get_snp_from_dbsnp_file(dbsnp_file):
    for line in dbsnp_file:
        fields = line.split('\t')
        if len(fields) < 3: continue
        rsid = int(fields[0])
        chrom = fields[1].rstrip('T')
        if chrom == 'PAR': continue
        start = int(fields[2]) + 1
        snp = {
            'xpos': get_xpos(chrom, start),
            'rsid': rsid
        }
        yield snp
