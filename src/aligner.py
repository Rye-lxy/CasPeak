#! /usr/bin/env python3
# Use LAST for sequence alignment

import os
import subprocess

from .fileReader import *

def aligner(args):
    # args: read, ref, insert, thread
    os.makedirs("lastdb", exist_ok=True)
    try:
        subprocess.check_call(f"lastdb -P{args.thread} -uRY4 lastdb/ref {args.ref}", shell=True)
        subprocess.check_call(f"lastdb -P{args.thread} -uRY4 lastdb/insert {args.insert}", shell=True)
        os.makedirs("lastal", exist_ok=True)
        with open("lastal/read_to_ref.maf", "w") as f:
            subprocess.check_call(f"last-train -P{args.thread} -Q0 lastdb/ref {args.read} | lastal -P{args.thread} --split -p - lastdb/ref {args.read}", 
                                stdout=f, shell=True)
        with open("lastal/read_to_insert.maf", "w") as f:
            subprocess.check_call(f"last-train -P{args.thread} -Q0 lastdb/insert {args.read} | lastal -P{args.thread} --split -p - lastdb/insert {args.read}", 
                                stdout=f, shell=True)
    except subprocess.CalledProcessError:
        print("Error in alignment", file=sys.stderr)
        exit(1)