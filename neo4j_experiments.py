#!/usr/bin/env python

import datetime
import pytz
import sys
from neomodel import config, db
from os.path import expanduser
import pprint
from BCBio.GFF import GFFExaminer
from BCBio import GFF
from combat_tb_neomodel.model.core import *

config.DATABASE_URL = 'bolt://neo4j:neo4j@localhost:7687'

def has_qualifier(feature, key, value):
    if key in feature.qualifiers and feature.qualifiers[key][0] == value:
        return True
    else:
        return False

def set_if_has(thing, feature, key, transform=None):
    if key in feature.qualifiers:
        value = feature.qualifiers[key][0]
        if transform is not None:
            value = transform(value)
        setattr(thing, key.lower(), feature.qualifiers[key][0])

def strip_id_colon(id):
    if ':' in id:
        parts = id.split(':')
        assert len(parts) == 2,\
                "Expected ID with 2 parts, like thing:ID, got " + id
        return parts[1]
    else:
        return id

def get_parent(feature, thing, lookup_dict):
    """Get parent object of a thing

    :param feature: Bio.SeqFeature.SeqFeature
    :param thing: combat_tb_neomodel.model.core.Feature
    :param lookup_dict: dict
    :return: combat_tb_neomodel.model.core.Feature
    """
    assert 'Parent' in feature.qualifiers,\
        "No Parent attribute found in feature for {} {}".format(
            thing.__class__.__name__,
            thing.uniquename
        )
    parent_id = feature.qualifiers['Parent'][0]
    assert parent_id in lookup_dict,\
        "Parent with ID {} not found for {} {} with biotype {}".format(
            parent_id,
            thing.__class__.__name__,
            thing.uniquename,
            thing.biotype
        )
    parent_thing = lookup_dict[parent_id]
    return parent_thing

mtb_gff_file = expanduser('~/Data/gff/Mycobacterium_tuberculosis_h37rv.GCA_000195955.2.30.gff3')
in_handle = open(mtb_gff_file)
# examiner = GFFExaminer()
# pprint.pprint(examiner.available_limits(in_handle))

gff_limits = dict(gff_type=('gene', 'transcript', 'pseudogene',
                            'CDS', 'tRNA_gene', 'ncRNA_gene',
                            'repeat_region'))
with db.transaction:
    tb = Organism(genus='Mycobacterium',
                  species='tuberculosis',
                  strain='H37Rv',
                  common_name='M. tuberculosis').save()

    record = next(GFF.parse(in_handle, target_lines=100))
    assert 'sequence-region' in record.annotations,\
        'Could not find sequence-region in GFF3 annotations'
    (name, start, end) = record.annotations['sequence-region'][0]
    chrom_length = end
    chrom = Chromosome(seqlen=chrom_length, uniquename='1')
    chrom.start = 0
    chrom.end = 0
    chrom.strand = '1'
    chrom.save()
    chrom.belongs_to.connect(tb)
    in_handle.seek(0)  # rewind file

    count = 0
    old_count = 1000
    quit = False
    gene_dict = dict()
    pseudogene_dict = dict()
    transcript_dict = dict()
    trna_dict = dict()
    ncrna_dict = dict()
    cds_dict = dict()
    repeat_region_count = 0
    while not quit:
        if count == old_count:
            break
        old_count = count
        for record in GFF.parse(in_handle, limit_info=gff_limits,
                                target_lines=1):
            count += 1
            feature = record.features[0]
            f_type = feature.type
            if 'ID' in feature.qualifiers:
                id = feature.qualifiers['ID'][0]
            elif f_type == 'repeat_region':
                repeat_region_count += 1
                id = 'repeat:' + str(repeat_region_count)
            else:
                exit("Feature found with no ID and not a repeat_region: " + repr(feature))
            start = feature.location.start
            end = feature.location.end
            strand = str(feature.location.strand)
            feature_length = end - start
            assert feature_length >= 1,\
                "Feature found with length < 1 after processing {} records".format(count)
            if f_type == 'gene':
                thing = Gene(uniquename=id)
                set_if_has(thing, feature, 'biotype')
                set_if_has(thing, feature, 'description')
                gene_dict[id] = thing
            elif f_type == 'pseudogene':
                if 'transcript_id' in feature.qualifiers:
                    thing = Transcript(uniquename=id)
                    set_if_has(thing, feature, 'biotype')
                else:
                    thing = PseudoGene(uniquename=id)
                    set_if_has(thing, feature, 'description')
                    pseudogene_dict[id] = thing
            elif (f_type == 'transcript'):
                thing = Transcript(uniquename=id)
                set_if_has(thing, feature, 'biotype')
                transcript_dict[id] = thing
            elif f_type == 'CDS':
                thing = CDS(uniquename=id)
                cds_dict[id] = thing
            elif f_type == 'tRNA_gene':
                thing = Trna(uniquename=id)
                trna_dict[id] = thing
            elif f_type == 'ncRNA_gene':
                thing = NCrna(uniquename=id)
                ncrna_dict[id] = thing
            elif f_type == 'repeat_region':
                thing = RepeatRegion(uniquename=id)
                set_if_has(thing, feature, 'description')
            else:
                print("Unknown feature type:", f_type, file=sys.stderr)
                continue

            set_if_has(thing, feature, 'Name', strip_id_colon)

            thing.seqlen = feature_length
            thing.timelastmodified = datetime.datetime.now(tz=pytz.utc)
            thing.start = start
            thing.end = end
            thing.strand = strand
            thing.save()
            if isinstance(thing, Transcript):
                if thing.biotype == 'protein_coding':
                    lookup_dict = gene_dict
                elif thing.biotype == 'tRNA':
                    lookup_dict = trna_dict
                elif thing.biotype == 'ncRNA':
                    lookup_dict = ncrna_dict
                elif thing.biotype == 'pseudogene':
                    lookup_dict = pseudogene_dict
                else:
                    exit("unknown biotype:" + thing.biotype)
                parent_thing = get_parent(feature, thing, lookup_dict)
                thing.part_of.connect(parent_thing)
            elif f_type == 'CDS':
                parent_transcript = get_parent(feature, thing, transcript_dict)
                thing.part_of.connect(parent_transcript)

            thing.location.connect(chrom)
            print('.', end='', file=sys.stderr)
            # print(thing.uniquename)
            # location = thing.location.connect(chrom)
            # location.start = start
            # location.end = end
            # location.strand = strand
            # location.save()

            thing.belongs_to.connect(tb)
            break
print(file=sys.stderr)
