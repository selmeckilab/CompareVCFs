import vcf
from vcf import utils
import argparse
import csv
import re
from intervaltree import IntervalTree
from fuzzywuzzy import process
from abc import ABCMeta, abstractmethod


class Genes(object):

    """
    The features class stores genes taken from a gff file
    """
    def __init__(self, gff):
        """

        :param gff: the file name of the gff file you want to import
        """
        self.gene_dict = gff

    @staticmethod
    def get_name(info):
        """
        import
        :param info:
        :return: Either the name of the Gene or the Orf number
        """
        match = re.search("Gene.+?;", info)
        if not match:
            match = re.search("Name.+?;", info)
        return match.group()[5:-1]

    @property
    def gene_dict(self):
        """

        :return: gene_dict
        """
        return self._gene_dict

    @gene_dict.setter
    def gene_dict(self, gff):
        """
        Sets the the features to a range tree of genes
        :param gff:

        """
        gd = {}
        with open(gff, "r") as infile:
            reader = csv.reader(infile, delimiter="\t")
            for r in reader:
                if ("#" not in r[0]) and (len(r) > 1) and (r[2] == "gene"):
                    gd.setdefault(r[0], IntervalTree())
                    gd[r[0]][int(r[3]):int(r[4])] = self.get_name(r[8])
        self._gene_dict = gd


class Child(object):
    def __init__(self, name, al, gt, gt_data, filter):
        """

        :param name:
        :param al:
        :param gt:
        :param gt_data:
        """
        self.name = name
        self.al = al
        self.gt = gt
        self.gt_data = gt_data
        self.filter = filter
        self.mutation = None
    @property
    def effect(self):
        #TODO
        """

        :rtype: str
        """
        return ""

    def __str__(self):
        return "Name:" + self.name + " Allele:" + str(self.al) + " Genotype:" + str(self.gt) + " MutationType:" + str(
            self.mutation)


class Row(object):

    features = None
    file_names = []
    parent_name = None

    def __init__(self):

        self.parent = self.extract_name(self.file_names[0])
        self.chrom = None
        self.pos = None
        self.children = []
        self.compare_children()

    @staticmethod
    def extract_name(file_name):
        """

        :param file_name:
        :return:
        """
        start_index = file_name.rfind("/") + 1
        end_index = file_name.find(".")
        return file_name[start_index:end_index]

    @property
    def gene(self):
        """
        The name of the gene that contains this position.
        :return: the gene that that contains this position. Returns a null string if the position is intergenic.
        """
        gene = ""
        feature = self.features[self.chrom][self.pos]
        if feature:
            gene = feature.pop()[2]
        return gene

    def compare_no_null(self, child):
        """
        compares the child if both the child and parent vary from the reference
        :param child:
        :return: str: hom, het or loh
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

    def compare_null_child(self):
        """
        compares the child if the child is homozygeous to the reference but the parent is not
        :return: str: hom or loh
        """
        mut = ""
        if self.parent.gt > 1:
            mut = "hom"
        else:
            mut = "loh"
        return mut

    @staticmethod
    def compare_null_parent(child):
        """
        compares the child if the parent is homozygeous to the reference but the child is not
        :param child:
        :return: str: hom or het
        """
        mut = ""
        if child.gt > 1:
            mut = "hom"
        else:
            mut = "het"
        return mut

    def compare_children(self):
        """
        determines the each child's mutation type and adds it to the child instance
        :return:
        """
        for child in self.children:
            if not child.filter:
                if self.parent and child:
                    child.mutation = self.compare_no_null(child)
                if self.parent and not child:
                    child.mutation = self.compare_null_child()
                if not self.parent and child:
                    child.mutation = self.compare_null_parent(child)


class MultiSampleRow(Row):
    def __init__(self, recs):
        """
        The multisample row represents a row taken from a vcf with multiple sample per file.
        :param recs:
        """
        Row.__init__(self)
        #if not raw_data[0].FILTER:
        #    self.make_children(parent=parent_call_name)

    def make_children(self, parent):
        """

        :param parent:
        :return:
        """
        # TODO
        pass


class MultiFileRow(Row):
    parent_index = 0

    def __init__(self, recs):
        """
        Represents a row created from multiple vcfs
        :param recs:
        """
        Row.__init__(self)

        if recs[self.parent_index]:
            parent_rec = recs[self.parent_index]
            name = self.extract_name(self.file_names[self.parent_index])
            al = parent_rec.samples[0].gt_bases[-1]
            gt = parent_rec.samples[0].gt_type
            gt_data = parent_rec.samples[0].data

            self.parent = Child(name, al, gt, gt_data, None)
            self.ref = parent_rec.REF
            self.pos = int(parent_rec.POS)
            self.chrom = parent_rec.CHROM
            self.make_children(recs)

        elif not recs[self.parent_index]:
            self.parent = None
            self.make_children(recs)

    def make_children(self, recs):
        """
        Adds a Child() for each child in the row.
        :param recs:
        :return:
        """

        for idx, rec in enumerate(recs):
            if idx != self.parent_index:
                name = self.extract_name(self.file_names[idx])
                if rec:
                    al = rec.samples[0].gt_bases[-1]
                    gt = rec.samples[0].gt_type
                    gt_data = rec.samples[0].data
                    filt = rec.FILTER
                    if not (self.pos and self.chrom):
                        self.pos = int(rec.POS)
                        self.chrom = rec.CHROM
                        self.ref = rec.REF
                    self.children.append(Child(name, al, gt, gt_data,filt))


class VCFWriter:

    file_names = []
    readers = []

    def __init__(self, file_names, readers):
        """
        Writes the data to a vcf
        :param file_names:
        :param readers:
        """
        self.outfiles = [open(f[:-3]+".comparevcf.vcf", "w") for f in file_names[1:]]
        self.writers = [vcf.Writer(filename, readers[idx+1]) for idx, filename in
                        enumerate(self.outfiles)]

    def __enter__(self):
        """

        :return:
        """
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """
        Closes all the files
        :param exc_type:
        :param exc_val:
        :param exc_tb:
        :return:
        """
        for f in self.outfiles:
            f.close()

    def write_row(self, row):
        """
        writes a Row() to the vcf
        :param row:
        :return:
        """
        formatted_row = self.format_row(row)
        for idx, rec in enumerate(formatted_row):
            self.writers[idx].write_record(rec)


    @staticmethod
    def format_row(row):
        """

        :param row:
        :return: a list of PyVCF Record objects. There is a different record for every child.
        """
        new_recs = []
        for idx, child in enumerate(row.children):

            if child.mutation:
                alt = [vcf.model._Substitution([child.al])]

                if not row.parent:
                    ref = row.ref
                else:
                    ref = row.parent.al

                new_rec=vcf.model._Record(row.chrom, row.pos, ref, "", alt, "", "PASS",
                                         "", "", None)

                call = vcf.model._Call(new_rec, child.name, child.gt_data)

                new_rec.samples = [call]
                new_recs.append(new_rec)
        return new_recs


class CSVWriter():

    def __init__(self, outfile):
        """
        Writes the data to a csv file
        :param outfile:
        """
        if outfile:
            self.outfile = open(outfile, "w")
        else:
            self.outfile = open(Row.file_names[0][:-4]+".comparevcf.csv", "w")

        self.csv_writer = csv.writer(self.outfile)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """
        closes the outfile
        :param exc_type:
        :param exc_val:
        :param exc_tb:
        :return:
        """
        self.outfile.close()

    def write_row(self,row):
        """

        :param row:
        :return:
        """

        formatted_rows = self.format_row(row)
        self.csv_writer.writerows(formatted_rows)

    @staticmethod
    def format_row(row):
        """
        Returns a list of lists. Where each list represents the csv row for a child.
        :param row:
        :return:
        """
        formatted_rows = []
        for child in row.children:

            if child.mutation:
                csv_row = [row.parent_name, child.name, row.chrom, row.pos, row.gene, child.mutation]
                formatted_rows.append(csv_row)

        return formatted_rows


class Comparator():
    __metaclass__ = ABCMeta

    def __init__(self, args):
        """
        Abstract class that iterates through all the VCF and generates Rows
        :param args:
        """
        self.infiles = args.infiles
        self.location = args.location
        self.files = [open(f, "r") for f in args.infiles]
        self.readers = [self.get_location(vcf.Reader(f)) for f in self.files]

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        """
        Closes all the files
        :param exc_type:
        :param exc_value:
        :param traceback:
        :return:
        """
        for f in self.files:
            f.close()

    def __iter__(self):
        return self

    @abstractmethod
    def next(self):

        pass

    def __next__(self):
        return self.next()

    def get_location(self, reader):
        """
        Uses an index to select only a portion of the genome
        :param reader:
        :return: a new reader that one contains a specific portion of the genome
        """
        adjusted_reader = reader
        if self.location:
            loc_chr = self.location[0]
            start = int(self.location[1])
            end = int(self.location[2])
            adjusted_reader = reader.fetch(loc_chr, start, end)
        return adjusted_reader


class MultiSampleComparator(Comparator):

    def __init__(self, args):
        """

        :param args:
        """
        Comparator.__init__(self, args)

    def next(self):
        """

        :return: a MultiSampleRow()
        """
        rec=next(self.readers[0])
        return MultiSampleRow(rec)


class MultiFileComparator(Comparator):

    def __init__(self, args):
        """

        :param args:
        """
        Comparator.__init__(self, args)
        self.walker = vcf.utils.walk_together(*self.readers)

    def next(self):
        """

        :return: a Row() creates from the walker's row
        """
        while True:
            recs = self.walker.next()
            if self.valid_recs(recs):
                row = MultiFileRow(recs)
                return row

    @staticmethod
    def valid_recs(recs):
        valid = True
        if not any(isinstance(x, vcf.model._Record) for x in recs):
            valid = False
        if recs[0] and recs[0].FILTER:
            valid = False
        return valid


def main():
    """

    :return:
    """
    parser = argparse.ArgumentParser(
        description="Compare one parent vcf to one or more child vcfs. All files should be bgzipped and tabix indexed.")
    parser.add_argument("infiles", nargs="+", help='''The locations of all files to be compared.
                        If only one file is given, the file is assumed to be a multisample vcf.
                        In this case, specify the call name of the parent. If there are multiple files,
                        the parent is assumed to be the first, unless a different file name is given by the parent
                        option.''')
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
    path = "/Users/cwf08523/PycharmProjects/CompareVCFs/test_data/"
    args = parser.parse_args([path + "SC5314_parental_trimmed_sorted_rmdup_realigned.bam.gatk.vcf.gz",
                              path + "aleeza1_trimmed_sorted_rmdup_realigned.bam.gatk.vcf.gz",
                              path + "aleeza3_trimmed_sorted_rmdup_realigned.bam.gatk.vcf.gz",
                              "--parent",
                              "SC5314",
                              "--gff",
                              path + "C_albicans_SC5314_version_A21-s02-m09-r08_features_with_chromosome_sequences.gff",
                              "--csv",])

    Row.parent_name = args.parent
    Row.file_names = args.infiles

    if args.gff:
        feat = Genes(args.gff)
        Row.features = feat.gene_dict

    if len(args.infiles) > 1:
        MultiFileRow.parent_index = args.infiles.index(process.extractOne(args.parent, args.infiles)[0])
        MultiFileRow.file_names = args.infiles
        comp = MultiFileComparator
    else:
        comp = MultiSampleComparator

    with comp(args) as cp:

        if args.csv:
            writer = CSVWriter
            writer_args = [args.outfile]
        else:
            writer = VCFWriter
            writer_args = [args.infiles, cp.readers]
        with writer(*writer_args) as wr:
            for row in cp:
                row.compare_children()
                wr.write_row(row)

if __name__ == "__main__":
    main()