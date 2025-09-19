#!/usr/bin/env python2
import os
import sys
import json

def extract_info(json_file):
    try:
        with open(json_file, 'r') as f:
            data = json.load(f)
            pval = data["test results"]["p-value"]
            kval = data["test results"]["relaxation or intensification parameter"]
            ogg_id = os.path.basename(json_file).split(".")[0]
            return (ogg_id, pval, kval)
    except:
        return None

def main(folder_path):
    print("ogg_id,pval,kval")
    for filename in os.listdir(folder_path):
        if filename.endswith(".json"):
            full_path = os.path.join(folder_path, filename)
            result = extract_info(full_path)
            if result:
                print("%s,%s,%s" % result)

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python get_relax_pvalue_kvalue.py folder/")
    else:
        main(sys.argv[1])
