'''
    Merge two files by column name
'''

import argparse
import pandas as pd

def Merge2File(file1, file2):
    df1 = pd.read_csv(file1, delimiter = "\t")
    df2 = pd.read_csv(file2, delimiter = "\t")
    df = pd.concat([df1, df2])
    column_intersect = [ x for x in list(df2.columns.values) if x in list(df1.columns.values) ]
    column_intersect.reverse()
    df.to_csv("merged_file.csv", index=False, columns = column_intersect, sep = "\t")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Build tree for vis.')
    parser.add_argument('--file1')
    parser.add_argument('--file2')
    arg = parser.parse_args()

    Merge2File(arg.file1, arg.file2)
