import config
import json
import config
import psycopg2
import pandas as pd
from collections import deque
import bisect


class Position:
    def __init__(self, pos, gap):
        self.pos = pos
        self.gap = gap

class LOHTracks(object):
    def __init__(self, parent, child, chromosome, parent_positions, child_positions, conn):
        # type: (str, str, str, list, list) -> object
        self.parent = parent
        self.child = child
        self.chromosome = chromosome
        self.parent_positions = deque(sorted(parent_positions))
        self.child_positions = deque(sorted(child_positions))
        self.conn=conn
        self.cursor = conn.cursor()
        self.compute_tracks()

    def gap_jump(self, x):
        'Find leftmost item greater than or equal to x'
        i = bisect.bisect_left(self.parent_positions, x)
        if i != len(self.parent_positions):
            m = len(self.parent_positions)-i+1
            self.parent_positions = deque(self.parent_positions, maxlen=m)
            return self.parent_positions.popleft()
        raise ValueError

    def walk(self):
        parent_pos = self.gap_jump(self.parent_positions.popleft())
        child_pos = self.child_positions.popleft()
        gap = 0
        while len(self.child_positions) > 1 and len(self.parent_positions) > 1:
            if gap > 10:
                try:
                    child_pos = self.child_positions.popleft()
                    parent_pos = self.gap_jump(child_pos)
                    gap = 0
                except ValueError:
                    pass
            if parent_pos == child_pos:
                yield child_pos
                child_pos = self.child_positions.popleft()
                parent_pos = self.parent_positions.popleft()
                gap = 0
            else:
                parent_pos = self.parent_positions.popleft()
                gap += 1

    def compute_tracks(self):
        first_pos = None
        last_pos = None
        for pos in self.walk():
            if first_pos is None:
                first_pos = pos
                last_pos = pos
            elif (pos - last_pos) < 10000:
                last_pos = pos
            else:
                if last_pos - first_pos > 0:
                    self.to_db(first_pos, last_pos)
                first_pos = None
                last_pos = None

        self.cursor.close()
        self.conn.commit()


    def to_db(self, start, stop):
        position_range = '[' + str(start) + ',' + str(stop) + ']'
        vals = (self.child,self.chromosome,position_range)
        self.cursor.execute(
            "INSERT INTO public.loh_tracks(child, chromosome, track) VALUES (%s, %s, %s);",
            vals)


def main():
    temp_pos = json.load(open('/Users/curtisfocht/AleezaAnalysis/het_pos.json', 'r'))

    het_pos = {}
    for key, values in temp_pos.items():
        start_index = key.rfind('/') + 1
        end_index = key.find('.')
        newkey = key[start_index:end_index]
        het_pos[newkey] = values
    try:
        conn = psycopg2.connect(
            "dbname='postgres' host='"+config.host+"' user='" + config.username + "' password='" + config.password + "'")
    except psycopg2.OperationalError:
        print "Unable to connect to the database"

    df = pd.read_sql_query("SELECT * FROM varscan_strains WHERE somatic_status = 'LOH'", con=conn)
    grouped = df.groupby(['parent', 'child', 'chromosome'])
    for group, data in grouped:
        parent = group[0]
        child = group[1]
        chromosome = group[2]
        parent_positions=het_pos[parent][chromosome]
        print 'Getting tracks for:', child, chromosome
        tracks = LOHTracks(parent, child, chromosome,parent_positions, data.position,conn)


    conn.close()

if __name__ == "__main__":
    main()