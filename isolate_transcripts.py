#!/usr/bin/env python3
import argparse
import errno
import os
from operator import itemgetter

from utils import reverse_complement, parse_rest, parse_fasta, merge_intervals, collapse_N

def parse_gtf(path):
    trs = {}
    rest = []
    with open(path, 'r') as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            data = line.split('\t')

            fields = {
                'scaffold': data[0],
                'feature': data[2],
                'start': data[3],
                'end': data[4],
                'strand': data[6],
            }
            if fields['feature'] in ['exon', 'UTR', 'stop_codon']:
                fields['rest'] = parse_rest(data[8])
                tr = fields['rest']['transcript_id']
                if tr not in trs:
                    trs[tr] = {
                            'name': tr,
                            'exons': [fields],
                            'scaffold': fields['scaffold'],
                            'strand': fields['strand']
                    }
                else:
                    trs[tr]['exons'].append(fields)
    return trs

def split_assembled_genome(gtf_path, fasta_path, od='.', L=98):
    trs = parse_gtf(gtf_path)
    print('GTF parsed')
    scaffolds = parse_fasta(fasta_path)
    print('Scaffolds parsed')

    wrong_scaffolds = 0
    for tr, data in trs.items():
        processed = []
        unprocessed = []
        if data['scaffold'] in scaffolds:
            sequence = scaffolds[data['scaffold']]
        else:
            print(f'{data["scaffold"]} not in FASTA file {fasta_path}')
            wrong_scaffolds += 1
            continue
        for exon in data['exons']:
            processed.append((int(exon['start']) - 1, int(exon['end']) -1))
            unprocessed.append((int(exon['start']) - L, int(exon['end']) + L - 2))

        processed = list(merge_intervals(processed))
        unprocessed = merge_intervals(unprocessed)

        processed = ''.join(map(lambda iv: sequence[iv[0]:iv[1]], processed))
        splice_junctions = []
        for iv in unprocessed:
            if iv[1] - iv[0] < 3 * L - 3:
                # If the length of the exon is < L-1
                splice_junctions.append(iv)
            else:
                splice_junctions.append((iv[0], iv[0] + 2 * L - 2))
                splice_junctions.append((iv[1] - (2 * L) + 2, iv[1]))

        splice_junctions = [
            collapse_N(sequence[iv[0]:iv[1]].upper()) for iv in splice_junctions
        ]

        processed = collapse_N(processed).upper()
        if data['strand'] == '-':
            processed = reverse_complement(processed)
            splice_junctions = map(reverse_complement, splice_junctions)

        with open(f'{od}/processed_transcripts.fasta', 'a') as fh:
            fh.write(f'>{tr}\n')
            fh.write('\n'.join([processed[i:i+80] for i in range(0, len(processed), 80)]))
            fh.write('\n')

        with open(f'{od}/splice_junctions.fasta', 'a') as fh:
            for i, sj in enumerate(splice_junctions):
                fh.write(f'>{tr}:{i}\n')
                fh.write(f'{sj}\n')
    print('DONE!')
    print(f'{wrong_scaffolds} scaffolds were not found, and the corresponding annotations were ignored.')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="""Given a GTF annotation
	and a corresponding fasta file, prints out a file containing the
	processed transcripts, and a window of size 2L-2 around each splice
	junction of each transcript.""")
    parser.add_argument('gtf', metavar='annotation.gtf', type=str,
                        help='Path to the GTF annotation file.')
    parser.add_argument('fasta', metavar='scaffolds.fasta', type=str,
                        help='Path to the GTF annotation file.')
    parser.add_argument('od', metavar='output_directory', type=str,
                        help='Directory in which results will be written.')
    parser.add_argument('-L', '--length', type=int, default=98, help='Read length.')
    args = parser.parse_args()

    try:
        os.mkdir(args.od)
    except OSError as exc:
        if exc.errno != errno.EEXIST:
            raise
        pass
    
    split_assembled_genome(args.gtf, args.fasta, args.od, args.length)
