import csv
import psycopg2
import config
import codecs

from astropy.io import ascii

def main():

    #try:
    #    conn = psycopg2.connect(
    #        "dbname='postgres' host='localhost' user='" + config.username + "' password='" + config.password + "'")
    #except:
    #    print "I am unable to connect to the database"

    #c = conn.cursor()
    ascii.open("/Users/curtisfocht/Downloads/AllMyFiles/PycharmProjects/CompareVCFs/test_data" /
            "/C_albicans_SC5314_version_A21-s02-m09-r08_features_with_chromosome_sequences.gff")

    #with open(
    #        "/Users/curtisfocht/Downloads/AllMyFiles/PycharmProjects/CompareVCFs/test_data"
    #        "/C_albicans_SC5314_version_A21-s02-m09-r08_features_with_chromosome_sequences.gff") as infile:

    #    myReader = UnicodeReader((line for line in infile if not line.startswith('#')), delimiter='\t', encoding='ascii')

        #for row in myReader:
        #    pass
        #    if row[0] == "##FASTA":
        #        break
        #    newrow = [i if i != '.' else None for i in row[:-1]]
        #    desired_attributes = ['ID', 'Name', 'Gene', 'Parent', 'Note']
        #    attributes = row[-1].split(';')
        #    print attributes
        #    for attribute in desired_attributes:
        #        if attribute in attributes:
        #            newrow[newrow.find('='):]
        #        else:
        #            newrow.append(None)
        #    position_range = '[' + newrow[3] + ',' + newrow[4] + ']'
        #    vals = tuple(newrow[:3] + [position_range] + newrow[5:])
        #    c.execute(
        #        "INSERT INTO public.features(chromosome, source, feature, positionrange, score,"
        #        "strand, frame, id, featurename, gene, parent, note) VALUES"
        #        "(%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s);",
        #        vals)

    #conn.commit()


if __name__ == "__main__":
    main()