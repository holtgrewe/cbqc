#!/usr/bin/env python
"""Collecting quality reports from new pipeline."""

from __future__ import with_statement

import argparse
import pkg_resources
import sys

import jinja2


class AlignmentStatsFile(object):
    """Information stored in a .sta file.

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
    def meanCoverage(self):
        if not self.totalBedLength or not self.mappedNucleotides:
            return 0.0
        else:
            return 1.0 * self.mappedNucleotides / self.totalBedLength

    def parse(self, f):
        """Parse statistics from given file."""
        # first line is path to BED file
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


class TrafficLights(object):
    """Handler for traffic lights thresholds."""

    def __init__(self, args):
        self.args = args
        self.thresholds = {
            'mean_coverage': {'green': args.mean_coverage_green,
                              'yellow': args.mean_coverage_yellow},
            'target_coverage_10x': {'green': args.target_coverage_10x_green,
                                    'yellow': args.target_coverage_10x_yellow},
            'target_coverage_20x': {'green': args.target_coverage_20x_green,
                                    'yellow': args.target_coverage_20x_yellow},
            }

    def colorFor(self, category, value):
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

    def colorFunc(self):
        """Return a function that returns the color for a value in a category."""
        def result(category, value):
            return self.colorFor(category, value)
        return result

    @classmethod
    def addArguments(self, parser):
        """Add arguments to parser."""
        parser.add_argument('--mean-coverage-green', type=float, default=0.85,
                            help='lower bound for green color on mean coverage')
        parser.add_argument('--mean-coverage-yellow', type=float, default=0.70,
                            help='lower bound for green color on mean coverage')

        parser.add_argument('--target-coverage-10x-green', type=float, default=0.85,
                            help='lower bound for green color on 10x target coverage')
        parser.add_argument('--target-coverage-10x-yellow', type=float, default=0.80,
                            help='lower bound for green color on 10x target coverage')

        parser.add_argument('--target-coverage-20x-green', type=float, default=0.85,
                            help='lower bound for green color on 20x target coverage')
        parser.add_argument('--target-coverage-20x-yellow', type=float, default=0.80,
                            help='lower bound for green color on 20x target coverage')


def renderReport(args, traffic_lights, sta):
    """Render reports as configured in args from the collected statistics."""
    # load HTML template file
    html = pkg_resources.resource_string('cbpipeline', 'qc_report.html')
    # render template, logic is there
    tpl = jinja2.Template(html)
    args.out_html.write(tpl.render(sta=sta, color=traffic_lights.colorFunc()))


def main():
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('--sta-file', type=argparse.FileType('rb'),
                        help='Picard statistics file', required=True)
    parser.add_argument('--out-html', type=argparse.FileType('wb'),
                        help='output HTML file', required=True)

    TrafficLights.addArguments(parser)

    args = parser.parse_args()
    sta = AlignmentStatsFile.loadFromFile(args.sta_file)

    renderReport(args, TrafficLights(args), sta)


if __name__ == '__main__':
    sys.exit(main())
