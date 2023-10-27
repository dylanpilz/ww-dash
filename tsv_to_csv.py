import csv

tsv_file = open("app/aggregate_variants.tsv", 'r')
csv_file = open("app/aggregate_variants.csv", 'w')

tsv_reader = csv.reader(tsv_file, delimiter='\t')
csv_writer = csv.writer(csv_file)

for row in tsv_reader:
    csv_writer.writerow(row)

tsv_file.close()
csv_file.close()