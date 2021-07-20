import argparse
import scanpy as sc
from diopy import output
import re

parser = argparse.ArgumentParser(description='R read the h5ad ', epilog="diorpy -f file.h5ad -t seurat")
parser.add_argument('-f', '--file', dest = 'file', help = 'The h5ad file', type = str, required = True)
parser.add_argument('-a', '--assay_name', dest = 'assay_name', help = 'The data type assay', type = str, required = True)
args = parser.parse_args()

adata = sc.read(args.file)
output.write_h5(adata = adata,file=re.sub('.h5ad', '_tmp.h5', args.file), assay_name=args.assay_name)

