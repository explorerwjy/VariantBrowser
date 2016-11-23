#===============================================================================================
# Variant Browser v1
#===============================================================================================
import pymongo
import os
from flask import Flask, request, session, g, redirect, url_for, abort, render_template, flash, jsonify, send_from_directory
import time
import gzip
from flask import Response
from collections import defaultdict, OrderedDict
from werkzeug.contrib.cache import SimpleCache
import json
import itertools
import random
import sys
from multiprocessing import Process
import glob	

from parsing import *
import lookups


app = Flask(__name__)

DATA_FILES_DIRECTORY = '../data/'

app.config['COMPRESS_DEBUG'] = True
# Load default config and override config from an environment variable
app.config.update(dict(
    DB_HOST='localhost',
    DB_PORT=27017,
    DB_NAME='VariantBrowser', 
    DEBUG=True,
    SECRET_KEY='development key',
    LOAD_DB_PARALLEL_PROCESSES = 4,  # contigs assigned to threads, so good to make this a factor of 24 (eg. 2,3,4,6,8)
    SITES_VCFS=glob.glob(os.path.join(os.path.dirname(__file__), DATA_FILES_DIRECTORY, 'Variants.vcf.txt')),
    GENCODE_GTF=os.path.join(os.path.dirname(__file__), DATA_FILES_DIRECTORY, 'Homo_sapiens.GRCh38.86.gtf.gz'),
    CANONICAL_TRANSCRIPT_FILE=os.path.join(os.path.dirname(__file__), DATA_FILES_DIRECTORY, 'canonical_transcripts.txt.gz'),
    OMIM_FILE=os.path.join(os.path.dirname(__file__), DATA_FILES_DIRECTORY, 'omim_info.txt.gz'),
    DBNSFP_FILE=os.path.join(os.path.dirname(__file__), DATA_FILES_DIRECTORY, 'dbNSFP2.6_gene.gz')
))


#===============================================================================================
# Connect to mongoDB
#===============================================================================================
def connect_db():
    client = pymongo.MongoClient(host=app.config['DB_HOST'], port=app.config['DB_PORT'])
    return client[app.config['DB_NAME']]


#===============================================================================================
# Load variants from vcf file
#===============================================================================================
def load_variants_file():
    def load_variants(sites_file, db):
        variants_generator = get_variants_from_sites_vcf(sites_file)
        try:
            db.variants.insert(variants_generator, w=0)
        except pymongo.errors.InvalidOperation:
            pass  # handle error when variant_generator is empty

    db = get_db()
    db.variants.drop()
    print("Dropped db.variants")

    # grab variants from sites VCF
    db.variants.ensure_index('xpos')
    db.variants.ensure_index('xstart')
    db.variants.ensure_index('xstop')
    db.variants.ensure_index('rsid')
    db.variants.ensure_index('genes')
    #db.variants.ensure_index('transcripts')

    sites_vcfs = app.config['SITES_VCFS']
    if len(sites_vcfs) == 0:
        raise IOError("No vcf file found")
    elif len(sites_vcfs) > 1:
        raise Exception("More than one sites vcf file found: %s" % sites_vcfs)

    load_variants(sites_vcfs, db)
    print "Done"
    return []

#===============================================================================================
# Load Genes from GTF file
#===============================================================================================
def load_gene_models():
    db = get_db()

    db.genes.drop()
    db.transcripts.drop()
    db.exons.drop()
    print 'Dropped db.genes, db.transcripts, and db.exons.'


    print "Loading Meta data"
    start_time = time.time()
    canonical_transcripts = {}
    with gzip.open(app.config['CANONICAL_TRANSCRIPT_FILE']) as canonical_transcript_file:
        for gene, transcript in get_canonical_transcripts(canonical_transcript_file):
            canonical_transcripts[gene] = transcript

    omim_annotations = {}
    with gzip.open(app.config['OMIM_FILE']) as omim_file:
        for fields in get_omim_associations(omim_file):
            if fields is None:
                continue
            gene, transcript, accession, description = fields
            omim_annotations[gene] = (accession, description)

    dbnsfp_info = {}
    with gzip.open(app.config['DBNSFP_FILE']) as dbnsfp_file:
        for dbnsfp_gene in get_dbnsfp_info(dbnsfp_file):
            other_names = [other_name.upper() for other_name in dbnsfp_gene['gene_other_names']]
            dbnsfp_info[dbnsfp_gene['ensembl_gene']] = (dbnsfp_gene['gene_full_name'], other_names)

    print 'Done loading metadata. Took %s seconds' % int(time.time() - start_time)


    # grab genes from GTF
    print "Loading Genes"
    start_time = time.time()
    with gzip.open(app.config['GENCODE_GTF']) as gtf_file:
        for gene in get_genes_from_gencode_gtf(gtf_file):
            gene_id = gene['gene_id']
            if gene_id in canonical_transcripts:
                gene['canonical_transcript'] = canonical_transcripts[gene_id]
            if gene_id in omim_annotations:
                gene['omim_accession'] = omim_annotations[gene_id][0]
                gene['omim_description'] = omim_annotations[gene_id][1]
            if gene_id in dbnsfp_info:
                gene['full_gene_name'] = dbnsfp_info[gene_id][0]
                gene['other_names'] = dbnsfp_info[gene_id][1]
            db.genes.insert(gene, w=0)

    print 'Done loading genes. Took %s seconds' % int(time.time() - start_time)

    start_time = time.time()
    db.genes.ensure_index('gene_id')
    db.genes.ensure_index('gene_name_upper')
    db.genes.ensure_index('gene_name')
    db.genes.ensure_index('other_names')
    db.genes.ensure_index('xstart')
    db.genes.ensure_index('xstop')
    print 'Done indexing gene table. Took %s seconds' % int(time.time() - start_time)


    # and now transcripts
    print "Loading transcripts"
    start_time = time.time()
    with gzip.open(app.config['GENCODE_GTF']) as gtf_file:
        db.transcripts.insert((transcript for transcript in get_transcripts_from_gencode_gtf(gtf_file)), w=0)
    print 'Done loading transcripts. Took %s seconds' % int(time.time() - start_time)

    start_time = time.time()
    db.transcripts.ensure_index('transcript_id')
    db.transcripts.ensure_index('gene_id')
    print 'Done indexing transcript table. Took %s seconds' % int(time.time() - start_time)

    # Building up gene definitions
    print "Loading exons"
    start_time = time.time()
    with gzip.open(app.config['GENCODE_GTF']) as gtf_file:
        db.exons.insert((exon for exon in get_exons_from_gencode_gtf(gtf_file)), w=0)
    print 'Done loading exons. Took %s seconds' % int(time.time() - start_time)

    start_time = time.time()
    db.exons.ensure_index('exon_id')
    db.exons.ensure_index('transcript_id')
    db.exons.ensure_index('gene_id')
    print 'Done indexing exon table. Took %s seconds' % int(time.time() - start_time)

    return []


#=====================================================================================
# Load database
# Currently only Variants and Genes
#=====================================================================================
def load_db():
    # Initialize database
    # Don't need to explicitly create tables with mongo, just indices
    confirm = raw_input('This will drop the database and reload. Are you sure you want to continue? [no] ')
    if not confirm.startswith('y'):
        print('Exiting...')
        sys.exit(1)
    all_procs = []
    for load_function in [load_variants_file, load_gene_models]:
        procs = load_function()
        all_procs.extend(procs)
        print("Started %s processes to run %s" % (len(procs), load_function.__name__))

    [p.join() for p in all_procs]
    print('Done!')


#=====================================================================================
# Opens a new database connection if there is none yet for the
# current application context.
#=====================================================================================
def get_db():
    if not hasattr(g, 'db_conn'):
        g.db_conn = connect_db()
    return g.db_conn


@app.route('/')
def homepage():
    return render_template('homepage.html')

@app.route('/autocomplete/<query>')
def awesome_autocomplete(query):
    print query
    if not hasattr(g, 'autocomplete_strings'):
        g.autocomplete_strings = [s.strip() for s in open(os.path.join(os.path.dirname(__file__), 'autocomplete_strings.txt'))]
    suggestions = lookups.get_awesomebar_suggestions(g, query)
    return Response(json.dumps([{'value': s} for s in suggestions]),  mimetype='application/json')


@app.route('/awesome')
def awesome():
    db = get_db()
    query = request.args.get('query')
    datatype, identifier = lookups.get_awesomebar_result(db, query)

    print "Searched for %s: %s" % (datatype, identifier)
    if datatype == 'gene':
        return redirect('/gene/{}'.format(identifier))
    elif datatype == 'variant':
        return redirect('/variant/{}'.format(identifier))
    else:
        raise Exception

@app.route('/variant/<variant_str>')
def variant_page(variant_str):
    db = get_db()
    try:
        chrom, pos, ref, alt = variant_str.split('-')
        pos = int(pos)
        # pos, ref, alt = get_minimal_representation(pos, ref, alt)
        xpos = get_xpos(chrom, pos)
        variant = lookups.get_variant(db, xpos, ref, alt)

        if variant is None:
            variant = {
                'chrom': chrom,
                'pos': pos,
                'xpos': xpos,
                'ref': ref,
                'alt': alt
            }


        print 'Rendering variant: %s' % variant_str
        return render_template(
            'variant.html',
            variant=variant
        )
    except Exception:
        print 'Failed on variant:', variant_str, ';Error=', traceback.format_exc()
        abort(404)


@app.route('/gene/<gene_id>')
def gene_page(gene_id):
    return get_gene_page_content(gene_id)

def get_gene_page_content(gene_id):
    db = get_db()
    try:
        gene = lookups.get_gene(db, gene_id)
        if gene is None:
            abort(404)
        gene_name = gene['gene_name']
        print gene_name
        variants_in_gene = lookups.get_variants_in_gene(db, gene_name)
        transcripts_in_gene = lookups.get_transcripts_in_gene(db, gene_id)
        transcript_id = gene['canonical_transcript']
        transcript = lookups.get_transcript(db, transcript_id)
        variants_in_transcript = lookups.get_variants_in_transcript(db, transcript_id)
        add_transcript_coordinate_to_variants(db, variants_in_transcript, transcript_id)

        print gene
        print variants_in_gene
        print transcript
        print transcripts_in_gene
        print "variants_in_transcript\n",variants_in_transcript

        t = render_template(
            'gene.html',
            gene=gene,
            transcript=transcript,
            variants_in_gene=variants_in_gene,
            variants_in_transcript=variants_in_transcript,
            transcripts_in_gene=transcripts_in_gene
        )
    	print 'Rendering gene: %s' % gene_id
    	return t
    except Exception, e:
        print 'Failed on gene:', gene_id, ';Error=', traceback.format_exc()
        abort(404)


@app.route('/not_found/<query>')
@app.errorhandler(404)
def not_found_page(query):
    return render_template(
        'not_found.html',
        query=query
    ), 404


@app.route('/error/<query>')
@app.errorhandler(404)
def error_page(query):
    return render_template(
        'error.html',
        query=query
    ), 404



if __name__ == "__main__":
    app.run()









