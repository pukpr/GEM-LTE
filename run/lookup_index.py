import sys
import csv

def lookup_site_name(csv_file, search_index):
    with open(csv_file, newline='', encoding='utf-8') as f:
        reader = csv.reader(f)
        for row in reader:
            if row and str(row[0]).strip() == str(search_index):
                if int(search_index) > 0:
                    print('"'+row[1]+'"')
                else:
                    print(row[1])
                return
        print(f"No match found for index: {search_index}")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python lookup_site_name.py <index>")
        sys.exit(1)
    search_index = sys.argv[1]
    lookup_site_name("sorted_sites_pmsl.csv", search_index)
