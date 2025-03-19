#! /usr/bin/env python3
# Use LAST for sequence alignment

import os
import subprocess

from .fileReader import *

from .logConfigure import getLogger
logger = getLogger(__name__)

def aligner(args):
    # args: read, ref, insert, thread
    os.makedirs("lastdb", exist_ok=True)
    os.makedirs("lastal", exist_ok=True)
    try:
        logger.info("Start aligning to reference genome...")
        subprocess.check_call(f"lastdb -P{args.thread} -uRY4 lastdb/ref {args.ref}", shell=True)
        with open("lastal/read_to_ref.maf", "w") as f:
            subprocess.check_call(f"last-train -P{args.thread} -Q0 lastdb/ref {args.read} | lastal -P{args.thread} --split -p - lastdb/ref {args.read}", 
                                stdout=f, shell=True)
        logger.info("Start aligning to ME consensus sequence...")
        subprocess.check_call(f"lastdb -P{args.thread} -uRY4 lastdb/insert {args.insert}", shell=True)
        with open("lastal/read_to_insert.maf", "w") as f:
            subprocess.check_call(f"last-train -P{args.thread} -Q0 lastdb/insert {args.read} | lastal -P{args.thread} --split -p - lastdb/insert {args.read}", 
                                stdout=f, shell=True)
        logger.info("All alignments are done.")
    except subprocess.CalledProcessError:
        print("Error in alignment", file=sys.stderr)
        exit(1)