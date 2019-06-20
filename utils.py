import re
from operator import itemgetter

def reverse_complement(string):
    comp = {
        'A': 'T',
        'T': 'A',
        'G': 'C',
        'C': 'G',
        'M': 'K',
        'K': 'M',
        'Y': 'R',
        'R': 'Y',
        'S': 'S',
        'W': 'W',
        'N': 'N'
    }
    return ''.join(list(map(lambda c: comp[c], string))[::-1])

def parse_rest(rest, file_format='.gtf'):
    # Parses the 'Attributes' field of the gff
    if file_format == '.gff3':
        return dict(map(lambda x: x.rstrip('\n').split('='), rest.split(';')))
    elif file_format == '.gtf':
        # Sorry, this is a disgusting file format
        # Who the heck puts semicolons in a semicolon-delimited file?!
        PATTERN = re.compile(r'''((?:[^;"']|"[^"]*"|'[^']*')+)''')
        rest = PATTERN.split(rest.rstrip('\n'))[1::2]
        try:
            return dict(map(lambda x: x.rstrip('\n')\
                                       .lstrip(' ')\
                                       .replace('"', '')\
                                       .split(' ', 1), rest))
                            # rest.rstrip(';').split(';')[:-1]))
        except Exception as e:
            print([r.rstrip('\n').lstrip(' ').replace('"', '').split(' ', 1) for r in rest])
            print(rest)


def parse_fasta(path):
    scaffolds = {}
    sequences = []
    current = ''
    with open(path, 'r') as fh:
        for line in fh:
            line = line.rstrip('\n')
            if line.startswith('>'):
                if len(current) > 0:
                    scaffolds[current] = ''.join(sequences)
                sequences = []
                current = line.split()[0][1:]
            else:
                sequences.append(line)

        scaffolds[current] = ''.join(sequences)
    return scaffolds

def merge_intervals(intervals):
    sorted_intervals = sorted(intervals, key=itemgetter(0))
    if not sorted_intervals:  # no intervals to merge
        return
    low, high = sorted_intervals[0]
    for iv in sorted_intervals[1:]:
        if iv[0] <= high:
            high = max(high, iv[1])
        else:
            yield low, high
            low, high = iv
    yield low, high

def len_overlap(a, b):
    return max(0, min(a[1], b[1]) - max(a[0], b[0]))

def overlap(iv1, iv2):
    start = end = 0
    if iv2[0] <= iv1[0] <= iv2[1]:
        start = iv1[0]
    elif iv1[0] <= iv2[0] <= iv1[1]:
        start = iv2[0]
    if iv2[0] <= iv1[1] <= iv2[1]:
        end = iv1[1]
    elif iv1[0] <= iv2[1] <= iv1[1]:
        end = iv2[1]
    return (start, end)


def collapse_N(string, k=31):
    # Collapses arbitrary length regions of N-bases into regions of length k
    curr = False
    mod = []
    counter = 0
    for c in string:
        if c == 'N' and curr and counter >= k:
            continue
        elif c == 'N':
            curr = True
            counter += 1
        else:
            curr = False
            counter = 0
        mod.append(c)
    return ''.join(mod)
