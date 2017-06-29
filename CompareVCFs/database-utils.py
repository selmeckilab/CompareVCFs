#!/anaconda/bin/python
import vcf
import psycopg2
import csv
import glob
try:
    import config
except ImportError:
    print 'Please add the config file'
import argparse
from psycopg2.extensions import AsIs
from psycopg2.extras import execute_values


def sequenza_add_strain(conn, reads_file):
    # Connect to the database
    c = conn.cursor()

    # extract the strain name
    strain_name = reads_file[reads_file.rfind('/') + 1: reads_file.find('_trimmed')]
    print "Adding strain: " + strain_name

    # get the positions of all mutations in the strain
    c.execute(
        "SELECT chromosome, position FROM mutations WHERE child = (%s) AND chromosome = (%s);",
        [strain_name, 'Ca21chrR_C_albicans_SC5314']
    )
    positions = c.fetchall()
    print 'Updating ' + str(len(positions)) + ' positions.'
    next_position = positions.pop(0)
    # print "The first position is", next_position
    updated_positions = 0
    skipped_positions = 0
    with open(reads_file, 'r') as infile:
        reader = csv.DictReader(infile, delimiter='\t')
        for row in reader:
            if next_position[0] == row["chr"] and str(next_position[1]) == row["n_base"]:
                # print next_position
                vals = (str(row['A']), str(row['C']), str(row['T']), str(row['G']), strain_name, row['chr'],
                        str(row['n_base']))
                c.execute(
                    "UPDATE mutations SET areads = (%s), creads = (%S), treads = (%S), greads = (%S)"
                    "WHERE child = (%s) AND chromosome = (%s) AND position = (%s);",
                    vals)
                updated_positions += 1
                # print "The next position is:", next_position
                # print "We have " + str(len(positions)) + " to go."
                try:
                    next_position = positions.pop(0)
                except IndexError:
                    break
            elif next_position[0] == row["chr"] and next_position[1] < int(row["n_base"]):
                # print "skipping:", next_position
                # print "We have " + str(len(positions)) + " to go."
                skipped_positions += 1
                try:
                    next_position = positions.pop(0)
                except IndexError:
                    break
                    # print "The next position is: ", next_position
        print "Updated " + str(updated_positions) + " positions."
        print "Skipped " + str(skipped_positions) + " positions."
        allele_dict = {'A': 'areads', 'C': 'creads', 'G': 'greads', 'T': 'treads'}
        for allele, column_name in allele_dict.items():
            query = "UPDATE mutations SET percentalt = CASE WHEN (areads+treads+creads+greads) != 0 " \
                    "THEN (%s/(areads+treads+creads+greads)) ELSE NULL END WHERE childallele = (%S)"
            c.execute(query, (AsIs(column_name), allele))

    c.close()
    conn.commit()


def bamreadcount_add_strain(conn, reads_file):
    c = conn.cursor()
    strain_name = reads_file[reads_file.rfind('/') + 1: reads_file.find('_trimmed')]
    print "Adding strain: " + strain_name
    c.execute(
        "SELECT chromosome, position, childallele FROM mutations WHERE child = (%s) ORDER BY chromosome, position;",
        [strain_name]
    )
    positions = c.fetchall()

    c.execute(
        "SELECT count(*) FROM read_info WHERE child = (%s);",
        [strain_name]
    )
    already_added = c.fetchone()
    if already_added[0] == 0 and len(positions) > 0:
        print 'Updating ' + str(len(positions)) + ' positions.'
        next_position = positions.pop(0)
        mutations_values = []
        read_info_values = []
        with open(reads_file, 'r') as infile:
            reader = csv.reader(infile, delimiter='\t')
            for row in reader:
                if next_position[0] == row[0] and str(next_position[1]) == row[1]:
                    base_values = [strain_name, row[0], row[1]]
                    for info in row[5:10]:
                        split_info = info.split(':')
                        read_info_values.append(tuple(base_values + split_info))

                        if split_info[0] == next_position[2]:
                            mutations_values.append([row[3], float(split_info[1]) / float(row[3])] + base_values)
                    try:
                        next_position = positions.pop(0)
                    except IndexError:
                        break
        mutations_sql = "UPDATE mutations SET readdepth = vals.readdepth::INT, bamreads_percentalt = vals.percentalt " \
                        "FROM (VALUES %s) AS vals (readdepth,percentalt, child, chromosome, position)" \
                        "WHERE mutations.child = vals.child AND mutations.chromosome = vals.chromosome AND mutations.position = vals.position::INT"
        temp = "(%s, %s, %s, %s, %s)"
        execute_values(c, mutations_sql, mutations_values, template=temp)
        read_info_sql = "INSERT INTO read_info (child, chromosome, position, allele, count, avg_mapping_quality, avg_base_quality," \
                        " avg_se_mapping_quality, num_plus_strand, num_minus_strand, avg_pos_as_fraction," \
                        " avg_num_mismatches_as_fraction, avg_sum_mismatch_qualities, num_q2_containing_reads, " \
                        "avg_distance_to_q2_start_in_q2_reads,avg_clipped_length, avg_distance_to_effective_3p_end) " \
                        "VALUES %s"
        temp = "(%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s)"
        execute_values(c, read_info_sql, read_info_values, template=temp)

    c.close()
    conn.commit()



def vcf_add_strain(conn, filename):
    c = conn.cursor()
    child_name = filename[filename.rfind("/") + 1:filename.find("_trimmed")]
    print "Adding:", child_name
    reader = vcf.Reader(open(filename, "r"))
    for row in reader:
        vals = (child_name, row.CHROM, row.POS, row.REF, str(row.ALT[0]))
        c.execute(
            "INSERT INTO public.mutations(child,chromosome,position,parentallele,childallele) VALUES (%s, %s, %s, %s, %s);",
            vals)
        if "ANN" in row.INFO:
            for annotation in row.INFO["ANN"]:
                annotation_list = [None if a == '' else a for a in annotation.split("|")]
                vals = tuple(annotation_list + [child_name, row.CHROM, row.POS])
                c.execute(
                    "INSERT INTO public.annotations(allele, annotation, annotationimpact, genename, geneid, featuretype,"
                    "featureid, transcriptbiotype, rank, hgvsc, hgvsp, dnapos, cdspos, aapos, distance, notes,"
                    "child, chromosome, position) VALUES "
                    "(%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s);",
                    vals)
    c.close()
    conn.commit()

def process_gff(args):
    # TODO Add code to add gff data into db
    pass


def varscan_generator(filename):
    child_name = filename[filename.rfind("/") + 1:filename.find("_trimmed")]

    with open(filename, 'r') as infile:
        reader = csv.reader(infile)
        next(reader,None)
        for row in reader:
            newrow = [child_name] + row
            newrow[6] = float(newrow[6].strip('%'))/100
            newrow[10] = float(newrow[10].strip('%'))/100
            yield row


def process_varscan(conn,filename):
    c = conn.cursor()


    read_info_sql = "INSERT INTO varscan (child, chromosome, position, refernce, variant, normal_reads1, normal_reads2,"\
                    " normal_var_freq, normal_gt, tumor_reads1, tumor_reads2," \
                    " tumor_reads2, tumor_var_freq, tumor_gt, somatic_status, variant_p_value, somatic_p_value, " \
                    "tumor_reads1_plus, tumor_reads1_minus, tumor_reads2_plus, tumor_reads2_minus, normal_reads1_plus,"\
                    " normal_reads1_minus, normal_reads2_plus, normal_reads2_minus) "\
                    "VALUES %s"
    temp = "("+"%s,"* 23 + "%s)"
    execute_values(c, read_info_sql, varscan_generator(filename), template=temp)

def add_read_data(args):
    try:
        conn = psycopg2.connect(
            "dbname='postgres' host='"+config.host+"' user='" + config.username + "' password='" + config.password + "'")
    except psycopg2.OperationalError:
        print "Unable to connect to the database"
    if args.directory:
        infiles = glob.glob(args.name + '/*')
        for infile in infiles:
            args.func(conn, infile)
    else:
        infile = args.name
        args.func(conn, infile)
    conn.close()


def main():
    # TODO Add parser for database parameters
    # db_info_parser = argparse.ArgumentParser(description='Parses database info')

    parser = argparse.ArgumentParser(prog="data")
    subparsers = parser.add_subparsers(help='sub-command help')
    sequenza_parser = subparsers.add_parser('sequenzareadcount', help=' Adds read counts and alternate allele frequency'
                                                                      ' using sequenza utils output')
    sequenza_parser.add_argument('name', help='The name of the file or directory you want to process')
    sequenza_parser.add_argument('--directory', action='store_true', help='If this flag is present the input will be '
                                                                          'treated as a directory and all .allelecount'
                                                                          'files are added. Otherwise the input is '
                                                                          'assumed to be one file.')
    sequenza_parser.set_defaults(func=sequenza_add_strain)
    bam_readcount_parser = subparsers.add_parser('bamreadcount',
                                                 help='Adds mutect')
    bam_readcount_parser.add_argument('name', help='The name of the file or directory you want to process')
    bam_readcount_parser.add_argument('--directory', action='store_true',
                                      help='If this flag is present the input will be '
                                           'treated as a directory and all .vcf or .vcf.gz'
                                           'files will be added.')
    bam_readcount_parser.set_defaults(func=bamreadcount_add_strain)

    vcf_parser = subparsers.add_parser('vcf',
                                       help='Adds read counts and alternate allele frequency and other variables'
                                            'using bam-readcounts outpu')
    vcf_parser.add_argument('name', help='The name of the file or directory you want to process')
    vcf_parser.add_argument('--directory', action='store_true', help='If this flag is present the input will be '
                                                                     'treated as a directory and all .allelecount'
                                                                     'files are added. Otherwise the input is assumed '
                                                                     'to be one file.')

    varscan_parser = subparsers.add_parser('vcf',
                                       help='Adds both snps and indels from varscan. The files show share the same prefix'
                                            'the snp file should end with .snp and the indel file with .indel')
    varscan_parser.add_argument('name', help='The name of the file or directory you want to process')
    varscan_parser.add_argument('--directory', action='store_true', help='If this flag is present the input will be '
                                                                     'treated as a directory and all .allelecount'
                                                                     'files are added. Otherwise the input is assumed '
                                                                     'to be one file.')
    varscan_parser.set_defaults(func=vcf_add_strain)

    # TODO Make gff processer
    # gff_parser = subparsers.add_parser('gff', help="Adds gff data to DB")
    # gff_parser.add_argument('--directory', action='store_true')
    # gff_parser.set_defaults(func=process_gff())

    args = parser.parse_args()
    if args.func == bamreadcount_add_strain or args.func == sequenza_add_strain or args.func == vcf_add_strain:
        add_read_data(args)
    else:
        args.func(args)


if __name__ == "__main__":
    main()
