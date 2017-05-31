import vcf
import psycopg2
from glob import glob
import uuid

def add_file_to_db(cursor,filename):
    child_name = filename[filename.rfind("/") + 1:filename.find("_trimmed")]
    reader = vcf.Reader(open(filename, "r"))
    for row in reader:
        vals = (child_name, row.CHROM, row.POS, row.REF, str(row.ALT[0]))
        cursor.execute("INSERT INTO public.mutations(child,chromosome,position,parentallele,childallele) VALUES (%s, %s, %s, %s, %s);", vals)
        if "ANN" in row.INFO:
             for annotation in row.INFO["ANN"]:
                 annotation_list = [None if a == '' else a for a in annotation.split("|")]

                 vals=tuple(annotation_list + [child_name, row.CHROM, row.POS])

                 cursor.execute(
                      "INSERT INTO public.annotations(allele, annotation, annotationimpact, genename, geneid, featuretype,"
                      "featureid, transcriptbiotype, rank, hgvsc, hgvsp, dnapos, cdspos, aapos, distance, notes,"
                      "child, chromosome, position) VALUES "
                      "(%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s);",
                      vals)


def main():
    try:
         conn = psycopg2.connect("dbname='candidadb' host='localhost' user='cwf08523' password='Xrd7756Y'")
    except:
         print "I am unable to connect to the database"
    c = conn.cursor()

    files = glob('/Users/cwf08523/PycharmProjects/CompareVCFs/snpeff_mutect/*.vcf')

    for f in files:
        add_file_to_db(c,f)

    conn.commit()

if __name__ == "__main__":
    main()