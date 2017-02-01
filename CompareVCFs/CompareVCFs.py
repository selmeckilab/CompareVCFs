__author__ = "cwf08523"
__date__ = "$Jan 4, 2017 4:23:14 PM$"

import vcf
from vcf import utils
import argparse
import csv
import re
from intervaltree import IntervalTree
from fuzzywuzzy import process
from abc import ABCMeta, abstractmethod


class Features(object):
    def __init__(self, gff):
        self.features_dict = gff

    @staticmethod
    def get_name(info):
        match = re.search("Gene.+?;", info)
        if not match:
            match = re.search("Name.+?;", info)
        return match.group()[5:-1]

    @property
    def features_dict(self):
        return self._features_dict

    @features_dict.setter
    def features_dict(self, gff):
        # TODO make sure this contains all the genes
        fd = {}
        with open(gff, "r") as infile:
            reader = csv.reader(infile, delimiter="\t")
            for r in reader:
                if ("#" not in r[0]) and (len(r) > 1) and (r[2] == "gene"):
                    fd.setdefault(r[0], IntervalTree())
                    fd[r[0]][int(r[3]):int(r[4])] = self.get_name(r[8])
        self._features_dict = fd


class Child(object):
    def __init__(self, call, name=""):
        """

        :type call: object
        """
        if call.is_filtered:
            self.filtered = True
        self.name = call.sample
        if name:
            self.name = name
        self.mutation = None
        self.al = call.gt_bases[-1]
        self.gt = call.gt_type

    @property
    def effect(self):
        #TODO
        return ""

    def __str__(self):
        return "Name:" + self.name + " Allele:" + str(self.al) + " Genotype:" + str(self.gt) + " MutationType:" + str(
            self.mutation)


class Row(object):
    """A row represents all records at a given loci.
    The most attributes of row are data, transformed_data and indel_ref.
    Data is a list of tuples, each tuple represents the alternate allele and the genotype of a child.
    "." is written at avery data point that is identicle to the parent.
    transformed_data is the same as data, except that every data point
    belonging to the wrong report is filled in with an "!".
    the indel ref is a list, the same length as the data, it keeps tracks of which tuples represent indels
    and which do not. This is important for transforming the data.
    """
    parent_name = None
    features = None
    file_names = []

    def __init__(self):
        self.parent = None
        self.chrom = None
        self.pos = None
        self.children = []

    @staticmethod
    def extract_name(file_name):
        """Given a long file name, it extracts the essential name
        :type file_name: String
        """
        start_index = file_name.rfind("/") + 1
        end_index = file_name.find(".")
        return file_name[start_index:end_index]

    @property
    def gene(self):
        gene = ""
        feature = self.features[self.chrom][self.pos]
        if feature:
            gene = feature.pop()[2]
        return gene

    def compare_no_null(self, child):
        """

        :param child:
        :return:
        """
        mut = ""
        if self.parent.gt > 1:
            if child.gt > 1:
                mut = "hom"
            else:
                mut = "het"
        if self.parent.gt == 1:
            if child.gt > 1:
                mut = "het"
            else:
                mut = "loh"
        return mut

    @property
    def compare_null_child(self):
        mut = ""
        if self.parent.gt > 1:
            mut = "hom"
        else:
            mut = "loh"
        return mut

    @property
    def compare_null_parent(self, child):
        mut = ""
        if child.gt > 1:
            mut = "hom"
        else:
            mut = "het"
        return mut

    def compare_children(self):

        for child in self.children:
            if self.parent and child:
                mut = self.compare_no_null(child)
            if self.parent and not child:
                mut = self.compare_null_child
            if not self.parent and child:
                mut = self.compare_null_parent(child)
            child.mutation = mut

    def __str__(self):
        return "Reference:" + self.ref + " parent_al:" + self.parent_al + " parent_gt:" + str(self.parent_gt)


class MultiSampleRow(Row):
    def __init__(self, recs):
        Row.__init__(self)

    def make_row(self, raw_data, parent_call_name):
        if not raw_data[0].FILTER:
            self.make_children(parent=parent_call_name)

    def melt_vcf_row(self, rec, parent=""):
        if rec:
            for sample in rec.samples:
                if sample.name == self.parent:
                    self.parent = Child(sample)
                else:
                    self.children.append(Child(sample))

    def make_children(self, parent):
        pass


class MultiFileRow(Row):
    parent_index = -1

    def __init__(self, recs):
        Row.__init__(self)
        self.make_row(recs)

    def make_children(self, recs):
        for idx, rec in enumerate(recs):
            if rec and idx != self.parent_index:
                self.children.append(Child(rec.samples[0], name=self.extract_name(self.file_names[idx])))
                if not (self.pos and self.chrom):
                    self.pos = int(rec.POS)
                    self.chrom = rec.CHROM

    def make_row(self, recs):
        if recs[self.parent_index] and not recs[self.parent_index].FILTER:
            parent_rec = recs[self.parent_index]
            self.parent = Child(parent_rec.samples[0])
            self.pos = parent_rec.POS
            self.chrom = parent_rec.CHROM
            self.make_children(recs)

        elif not recs[self.parent_index]:
            self.parent = None
            self.make_children(recs)

    def melt_vcf_row(self, rec):
        # type: (vcf.Model.Record) -> None
        for sample in rec.samples:
            if sample.sample == self.parent:
                self.parent = Child(sample)
            else:
                self.children.append(Child(sample))


class CustomOpen:
    def __init__(self, infiles):
        self.files = [open(f, "r") for f in infiles]

    def __enter__(self):
        return self.files

    def __exit__(self, exc_type, exc_value, traceback):
        for f in self.files:
            f.close()


class Writer(object):
    __metaclass__ = ABCMeta

    def __init__(self, row_generator):
        self.row_generator = row_generator

    @abstractmethod
    def write_file(self):
        pass

    @abstractmethod
    def convert_row(self):
        pass


class VCFWriter(Writer):

    file_names = []
    readers = []

    def __init__(self):
        Writer.__init__(self)

        self.writers = [vcf.Writer(filename[:-4] + "_compared.vcf", self.readers[idx]) for idx, filename in
                        enumerate(self.filenames)]

    def write_file(self):
        for row in self.row_generator():
            recs = self.convert_row(row)
            for idx, rec in enumerate(recs):
                self.writers[idx].write_record(rec)

    def convert_row(self, rec):
        new_recs = []
        for idx, child in enumerate(rec.children):
            alt = [vcf.model._Substitution([child.al])]
            protein_effect = child.effect
            if child.effect:
                protein_effect = "ProteinEffect"+child.effect
            new_rec=vcf.model._Record(rec.chrom, rec.pos, "", rec.parent.al, alt, child.quality, "PASS",
                                      child.effect, None)
            # TODO format data
            # TODO check call name
            call = vcf.model._call(new_rec, child.name, child.gt)

            new_rec.samples = [call]
            new_recs.append(new_rec)
        return new_recs


class CSVWriter(Writer):

    def __init__(self):
        Writer.__init__(self)


    def make_csv(self, outfile):
        with open(outfile, "w")as out:
            writer = csv.writer(out)
            for row in self.row_generator():
                out_row = row.csv_row
                writer.writerows(out_row)

    def write_row(self):
        pass

    def convert_row(self, row):
        base_row = [row.chrom, row.pos, row.extract_name(row.file_names[row.parent_index])]
        rows = []
        for child in row.children:
            if not child.filtered:
                row = base_row + [child.name] + [child.mutation]
                if row.features:
                    row.append(row.gene)
                rows.append(row)
        return rows

class Comparator:
    __metaclass__ = ABCMeta

    def __init__(self, args, csv=False):
        self.infiles = args.infiles
        self.location = args.location
        if csv:
            self.writer = CSVWriter()
        else:
            self.writer = VCFWriter()

    def get_location(self, reader):
        adjusted_reader = reader
        if self.location:
            loc_chr = self.location[0]
            start = int(self.location[1])
            end = int(self.location[2])
            adjusted_reader = reader.fetch(loc_chr, start, end)
        return adjusted_reader

    def write_file(self):
        for row in self.row_generator:
            self.writer.write_row(row)

    @abstractmethod
    def row_generator(self):
        pass


class MultiSampleComparator(Comparator):
    def __init__(self, args):
        Comparator.__init__(self, args)

    def row_generator(self):
        with CustomOpen(self.infiles) as infiles:
            readers = [self.get_location(vcf.Reader(infile)) for infile in infiles]
            reader = readers[0]
            for rec in reader:
                yield MultiSampleRow(rec).compare_children()


class MultiFileComparator(Comparator):
    def __init__(self, args):
        # type: (list) -> None
        Comparator.__init__(self, args)

    def row_generator(self):
        with CustomOpen(self.infiles) as infiles:
            readers = [self.get_location(vcf.Reader(infile)) for infile in infiles]
            VcfWriter.readers=readers
            walker = utils.walk_together(*readers)
            for recs in walker:
                row = MultiFileRow(recs)
                row.compare_children()
                yield row


def main():
    parser = argparse.ArgumentParser(
        description="Compare one parent vcf to one or more child vcfs. All files should be bgzipped and tabix indexed.")
    parser.add_argument("infiles", nargs="+", help='''The locations of all files to be compared.
                        If only one file is given, the file is assumed to be a multisample vcf.
                        In this case, specify the call name of the parent. If there are multiple files,
                        the parent is assumed to be the first, unless a different file name is given by the parent option.''')
    parser.add_argument("--parent", help="The call name or filename of the parent")
    parser.add_argument("--pedigree", help=(
        "The relationship between parents and children can also be specified with a\n"
        "     json file of the form:\n"
        "    {parent:/file/path/to/parent.vcf,children:[file/path/to/first/child.vcf,file/path/to/second/child.vcf]}"))
    parser.add_argument("--location", nargs=3,
                        help='''The chromesome, start position and end position that you would like to pull out. For
                        example: 5 44567 44789. If you specify this option the file must be tabix indexed.''')
    parser.add_argument("--gff", default="", help='''Path of GFF for features''')
    parser.add_argument("--csv", action='store_true', help='''If you want a vcf output, otherwise will be a csv''')
    parser.add_argument("--outfile",
                        help='''Name of the outfile. By default the outfile is the name of the parent followed by the
                        extension''')
    path = "/Users/cwf08523/Desktop/test_data/"
    args = parser.parse_args([path + "SC5314_parental_trimmed_sorted_rmdup_realigned.bam.gatk.vcf.gz",
                              path + "aleeza1_trimmed_sorted_rmdup_realigned.bam.gatk.vcf.gz",
                              path + "aleeza3_trimmed_sorted_rmdup_realigned.bam.gatk.vcf.gz",
                              "--parent", "SC5314",
                              "--gff",
                              path + "C_albicans_SC5314_version_A21-s02-m09-r08_features_with_chromosome_sequences.gff",
                              ])

    if args.gff:
        feat = Features(args.gff)
        Row.features = feat.features_dict

    if len(args.infiles) > 1:
        MultiFileRow.parent_index = args.infiles.index(process.extractOne(args.parent, args.infiles)[0])
        MultiFileRow.file_names = args.infiles
        comp = MultiFileComparator(args)
    else:
        comp = MultiSampleComparator(args)

    if args.csv:
        writer = CSVWriter()
    else:
        writer = VCFWriter()

    writer.write_file(comp.row_generator())


if __name__ == "__main__":
    main()
