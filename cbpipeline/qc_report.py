#!/usr/bin/env python
"""Collecting quality reports from new pipeline."""

from __future__ import print_function

__author__ = 'Manuel Holtgrewe <manuel.holtgrewe@charite.de>'

import argparse
import os
import pkg_resources
import sys
import zipfile

import jinja2

class AlignmentStatsFile(object):
    """Information stored in a .sta file.

    @ivar file_name         path to .sta file
    @ivar bed               path to BED file for statistics
    @ivar reads             number of reads
    @ivar readsAligned      number of aligned reads
    @ivar readsUnaligned    number of unaligned reads
    @ivar duplicates        number of duplicates
    @ivar duplicationRate   percentage of duplicated reads
    @ivar mappedNucleotides total number of aligned nucleotides in target intervals
    @ivar fractionWithCov   dict mapping integer coverage lower bound to
                            float percentage
    @ivar totalBedLength    sum of interval lengths in BED file
    """

    def __init__(self):
        self.file_name = None
        self.bed = None
        self.reads = None
        self.readsAligned = None
        self.readsUnaligned = None
        self.duplicates = None
        self.duplicationRate = None
        self.mappedNucleotides = None
        self.fractionWithCov = {}
        self.totalBedLength = None

    @property
    def alignedRate(self):
        return 1.0 * self.readsAligned / self.reads

    @property
    def meanCoverage(self):
        if not self.totalBedLength or not self.mappedNucleotides:
            return 0.0
        else:
            return 1.0 * self.mappedNucleotides / self.totalBedLength

    def parse(self, f):
        """Parse statistics from given file."""
        # first line is path to BED file
        self.file_name = f.name
        line = f.readline().strip()
        if not line.startswith('#'):
            raise Exception('First line of .sta file must start '
                            'with hash (#) and contain path to BED file.')
        self.bed = line[1:]
        # parse alignment and duplication lines
        self.reads = self._parseLine(f.readline().strip(), 'reads:', int)
        self.readsAligned = self._parseLine(f.readline().strip(), 'reads aligned:', int)
        self.readsUnaligned = self._parseLine(f.readline().strip(), 'reads unaligned:', int)
        self.duplicates = self._parseLine(f.readline().strip(), 'duplicates:', int)
        self.duplicationRate = self._parseLine(f.readline().strip(), 'duplication rate:', float)
        self.mappedNucleotides = self._parseLine(
            f.readline().strip(), 'Mapped nucleotides on target:', int)
        # parse fraction lines
        for perc in range(10, 110, 10):
            key = 'Fraction with cov. >= %d :' % perc
            self.fractionWithCov[perc] = self._parseLine(f.readline().strip(), key, float)

    def _parseLine(self, line, starts_with, type_):
        """Parse one line, check that line starts with starts_with."""
        if not line.startswith(starts_with):
            raise Exception('Expected line for "%s" but got "%s"' %
                            (starts_with, line))
        return type_(line.split('\t')[1])

    @classmethod
    def loadFromFile(self, f):
        """Parse stats file and get interval lengths from BED file."""
        result = AlignmentStatsFile()
        result.parse(f)
        result.totalBedLength = totalBedIntervalLength(result.bed)
        return result


class BedRecord(object):
    """One entry from a BED file.

    @ivar ref       name of contig/chromosome
    @ivar begin_pos 0-based begin position
    @ivar end_pos   0-based end position
    """

    def __init__(self, ref=None, beginPos=None, endPos=None):
        self.ref = ref
        self.beginPos = beginPos
        self.endPos = endPos

    def null(self):
        return self.ref is None

    def overlapsWith(self, other):
        return (self.ref == other.ref and self.beginPos < other.endPos and
                other.beginPos < self.endPos)

    def merge(self, other):
        assert self.ref == other.ref
        self.beginPos = min(self.beginPos, other.beginPos)
        self.endPos = max(self.endPos, other.endPos)

    def length(self):
        return (self.endPos - self.beginPos)


def readBedFile(path, mergeOverlapping=False):
    """Read BED file and yield BedRecord by BedRecord.

    @param path             str with path to BED file to read
    @param mergeOverlapping bool, if True then overlapping adjacent records
                            are merged
    """
    current = BedRecord()
    with open(path, 'rb') as f:
        for line in f:
            line = line.strip()
            arr = line.split('\t', 3)
            ref, beginPos, endPos = arr[:3]
            record = BedRecord(ref, int(beginPos), int(endPos))
            if not mergeOverlapping:
                yield record
            else:
                if current.null():
                    current = record
                elif not current.null() and current.overlapsWith(record):
                    current.merge(record)
                else:
                    yield current
                    current = record
    if mergeOverlapping and not current.null():
        yield current


def totalBedIntervalLength(path):
    """Returns total length of all intervals in path."""
    result = 0
    for record in readBedFile(path):
        result += record.length()
    return result


class FastQCReport(object):
    """QC report key metrics.

    Valid keys are:

    basic_statistics, per_base_sequence_quality, per_tile_sequence_quality,
    per_sequence_quality_scores, per_base_sequence_content, per_sequence_gc_content,
    per_base_n_content, sequence_length_distribution, sequence_duplication_levels,
    overrepresented_sequences, adapter_content, kmer_content

    @ivar summary dict mapping the statistics names to "pass", "warn",
                  "fail", see above for valid keys
    """

    def __init__(self, file_name):
        self.file_name = file_name
        self.summary = {}

    def parse(self, report_str):
        """Parse from string with report."""
        for line in report_str.splitlines():
            if not line.startswith('>>') or line.startswith('>>END_MODULE'):
                continue
            line = line[2:]
            arr = line.split('\t')
            key = arr[0].lower().replace(' ', '_')
            value = arr[1].lower()
            self.summary[key] = value

    @classmethod
    def loadFromFile(self, path):
        """Load FastQC report from the ZIP file at path."""
        # load contents of fastqc_data.txt
        f = zipfile.ZipFile(path, 'r')
        for name in f.namelist():
            if name.endswith('fastqc_data.txt'):
                fastqc_data = f.read(name)
                break
        report = FastQCReport(os.path.basename(path))
        report.parse(fastqc_data)
        return report


class VcfStatsGenerator(object):
    """Helper class for building VcfStats objects."""

    def __init__(self, f):
        self.f = f
        self.vcf_stats = VcfStats()
        self.vcf_stats.file_name = f.name

    def run(self):
        """Perform parsing."""
        for line in self.f:
            if line.startswith('#'):
                continue  # skip header
            (chrom, pos, id_, ref, alt, qual, filter_,
             info_, format_, individual) = line.strip().split('\t')[:10]

            # split FORMAT field and get index of genotype
            fmt = dict([(key, i) for i, key in enumerate(format_.split(':'))])
            if fmt.get('GT') is None:
                print('WARNING: no genotype!', file=sys.stderr)
                continue  # skip records without calls
            # split INFO field
            info = set(info_.split(';'))
            # split call individual
            call = individual.split(':')
            # split alternatives
            alts = alt.split(',')

            # check that is variant
            gt = self._getGenotypes(ref, fmt, alts, call)
            if ref == gt[0] and ref == gt[1]:
                print('Warning: not a variant\n%s\n' % line,
                      file=sys.stderr)
                continue

            #print("%s\t%s" % (chrom, pos))

            self._countDbSnp(ref, fmt, info, alts, call)
            self._countHomHet(ref, fmt, alts, call)
            self._countVariantKind(ref, fmt, alts, call)
            self._countTrans(ref, fmt, alts, call)
            self._countJannovar(info)
        return self.vcf_stats

    def _getGenotypes(self, ref, fmt, alts, calls):
        """Return actual genotypes"""
        genotypes = map(int, calls[fmt['GT']].split('/'))
        #print(genotypes)
        return [([ref] + alts)[i] for i in genotypes]

    def _countDbSnp(self, ref, fmt, info, alts, calls):
        """Update db_snp_count, db_snp_rare_counts."""
        # check that record describes a SNV in our individual
        genotypes = self._getGenotypes(ref, fmt, alts, calls)
        min_len = min([len(x) for x in genotypes])
        max_len = max([len(x) for x in genotypes])
        if min_len != 1 or max_len != 1:
            return  # skip, is not a SNV
        
        self.vcf_stats.dbsnp_count += int('DB' in info)
        self.vcf_stats.no_dbsnp_count += int('DB' not in info)
        # TODO(holtgrew): annotation field for db_snp_rare_counts so far

    def _countHomHet(self, ref, fmt, alts, calls):
        """Update heterozygous, homozygous counts."""
        genotypes = self._getGenotypes(ref, fmt, alts, calls)
        is_hom = (genotypes[0] == genotypes[1])
        self.vcf_stats.heterozygous += int(not is_hom)
        self.vcf_stats.homozygous += int(is_hom)

    def _countVariantKind(self, ref, fmt, alts, calls):
        """Update snv_count, ins_count, del_count."""
        gts = self._getGenotypes(ref, fmt, alts, calls)
        for i, gt in enumerate(gts):
            if i == 1 and gts[0] == gts[1]:
                continue  # skip calling twice
            if gt == ref:
                continue  # ALT not differing from REF
            if len(gt) == 1 and len(ref) == 1:  # SNP
                self.vcf_stats.snv_count += 1
            elif len(gt) > len(ref):  # INSertion
                self.vcf_stats.ins_count += 1
            elif len(gt) < len(ref):  # DELetion
                self.vcf_stats.del_count += 1
            else:
                pass  # substitute more than one char, WARN?

    def _countJannovar(self, info):
        """Update counters based on Jannovar annotations."""
        for field in info:
            if field.startswith('EFFECT'):
                key, values = field.split('=', 1)
                for v in values.split(','):
                    self.vcf_stats.num_annos.setdefault(v, 0)
                    self.vcf_stats.num_annos[v] += 1

    def _countTrans(self, ref, fmt, alts, calls):
        """Update transitions, transversions."""
        TRANS = set(['GA', 'AG', 'TC', 'CT'])
        gts = self._getGenotypes(ref, fmt, alts, calls)
        for i, gt in enumerate(gts):
            if len(ref) != 1 or len(gt) != 1:
                continue  # ignore, not SNV
            if gt == ref:
                continue  # skip, not variant
            var = ('%s%s' % (gt, ref)).upper()
            self.vcf_stats.transitions += int(var in TRANS)
            self.vcf_stats.transversions += int(var not in TRANS)

class VcfStats(object):
    """Statistics on a VCF file.

    @ivar file_name        str, name of analyzed file
    @ivar no_dbsnp_count   number of SNVs not in dbSNP
    @ivar dbsnp_count      number of SNVs in dbSNP
    @ivar dbsnp_rare_count number of rare variants from dbSNP
    @ivar del_count        number of called deletions
    @ivar heterozygous     number of heterozygous variants
    @ivar homozygous       number of homozygous variants
    @ivar ins_count        number of called insertions
    @ivar num_annos        number of annotated 
    @ivar snv_count        number of SNVs
    @ivar transitions      number of transitions
    @ivar transversions     number of transversions
    """

    def __init__(self):
        self.file_name = None
        self.snv_count = 0
        self.no_dbsnp_count = 0
        self.dbsnp_count = 0
        self.dbsnp_rare_count = 0
        self.del_count = 0
        self.heterozygous = 0
        self.homozygous = 0
        self.ins_count = 0
        self.num_annos = {}
        self.snv_count = 0
        self.transitions = 0
        self.transversions = 0

    @classmethod
    def generateForFile(self, vcf_file):
        return VcfStatsGenerator(vcf_file).run()
    

class TrafficLights(object):
    """Handler for traffic lights thresholds."""

    def __init__(self, args):
        self.args = args
        self.thresholds = {
            'aligned_reads': {'green': args.aligned_reads_green,
                              'yellow': args.aligned_reads_yellow},
            'mean_coverage': {'green': args.mean_coverage_green,
                              'yellow': args.mean_coverage_yellow},
            'target_coverage_10x': {'green': args.target_coverage_10x_green,
                                    'yellow': args.target_coverage_10x_yellow},
            'target_coverage_20x': {'green': args.target_coverage_20x_green,
                                    'yellow': args.target_coverage_20x_yellow},
            }

    def colorPassFor(self, value):
        """Return color for pass, warn, and fail."""
        mapping = {'pass': 'green',
                   'warn': 'yellow',
                   'fail': 'red'}
        return mapping.get(value, 'black')

    def colorPassFunc(self):
        """Return a function that returns the color for pass, warn, fail."""
        def result(value):
            return self.colorPassFor(value)
        return result

    def colorThreshFor(self, category, value):
        """Return color for value and category, determined by threshold."""
        if not self.thresholds.get(category):
            return 'unknown category %s' % category
        if not self.thresholds[category].get('green'):
            return 'missing entry green for %s' % category
        if not self.thresholds[category].get('yellow'):
            return 'missing entry yellow for %s' % category
        if self.thresholds[category]['green'] <= value:
            return 'green'
        elif self.thresholds[category]['yellow'] <= value:
            return 'yellow'
        else:
            return 'red'

    def colorThreshFunc(self):
        """Return a function that returns the color for a value in a category, by threshold."""
        def result(category, value):
            return self.colorThreshFor(category, value)
        return result

    @classmethod
    def addArguments(self, parser):
        """Add arguments to parser."""
        group = parser.add_argument_group('quality thresholds')
        group.add_argument('--mean-coverage-green', type=float, default=0.85,
                           help='lower bound for green color on mean coverage')
        group.add_argument('--mean-coverage-yellow', type=float, default=0.70,
                           help='lower bound for green color on mean coverage')

        group.add_argument('--target-coverage-10x-green', type=float, default=0.85,
                           help='lower bound for green color on 10x target coverage')
        group.add_argument('--target-coverage-10x-yellow', type=float, default=0.80,
                           help='lower bound for green color on 10x target coverage')

        group.add_argument('--target-coverage-20x-green', type=float, default=0.85,
                           help='lower bound for green color on 20x target coverage')
        group.add_argument('--target-coverage-20x-yellow', type=float, default=0.80,
                           help='lower bound for green color on 20x target coverage')

        group.add_argument('--aligned-reads-green', type=float, default=0.97,
                           help='lower bound for green color on aligned reads percentage')
        group.add_argument('--aligned-reads-yellow', type=float, default=0.90,
                           help='lower bound for yellow color on aligned reads percentage')


def renderReport(args, traffic_lights, sta, fastqcs, vcf_stats):
    """Render reports as configured in args from the collected statistics."""
    def formatNumber(num):
        import locale
        locale.setlocale(locale.LC_ALL, 'en_US.utf-8')
        return locale.format('%d', num, 1)
    # load HTML template file
    html = pkg_resources.resource_string('cbpipeline', 'qc_report.html')
    # render template, logic is there
    env = jinja2.Environment(loader=jinja2.PackageLoader('cbpipeline', '.'))
    env.filters['format_number'] = formatNumber
    env.filters['basename'] = os.path.basename
    tpl = env.get_template('qc_report.html')
    args.out_html.write(tpl.render(sta=sta,
                                   fastqcs=fastqcs,
                                   vcf_stats=vcf_stats,
                                   color_pass=traffic_lights.colorPassFunc(),
                                   color_thresh=traffic_lights.colorThreshFunc(),
                                   black_and_white=args.black_and_white,
                                   version=pkg_resources.require('cbpipeline')[0].version
                               ))


def main():
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('--black-and-white', action='store_true', default=False,
                        help='do not use colors in output')

    file_group = parser.add_argument_group('input/output paths')
    file_group.add_argument('--sta-file', type=argparse.FileType('rb'),
                            help='Picard statistics file', required=True)
    file_group.add_argument('--fastqc-zip', type=str, required=True,
                            action='append', default=[], dest='fastqc_zips',
                            help='path to FASTQC zip file')
    file_group.add_argument('--vcf-file', type=argparse.FileType('rb'), required=True,
                            help='path to VCF file to analyze')
    file_group.add_argument('--out-html', type=argparse.FileType('wb'),
                            help='output HTML file', required=True)

    TrafficLights.addArguments(parser)

    args = parser.parse_args()
    sta = AlignmentStatsFile.loadFromFile(args.sta_file)
    fastqcs = [FastQCReport.loadFromFile(f) for f in args.fastqc_zips]
    vcf_stats = VcfStats.generateForFile(args.vcf_file)

    renderReport(args, TrafficLights(args), sta, fastqcs, vcf_stats)


if __name__ == '__main__':
    sys.exit(main())
