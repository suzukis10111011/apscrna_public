import os
#!/usr/bin/env python
"""Purpose: Compare gene-set coverage before and after symbol correction across validation cells
Input: frozen/net/collectri_regulons.csv, data/processed/all_cells_validation.h5ad
Output: results/tables/20b_geneset_coverage_before_after.csv

20b: symbol 修正の前後で、遺伝子セットの取りこぼし率を並べる。
"""
import sys, os, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _guard import require_heartbeat, W
require_heartbeat("20b")
import pandas as pd, scanpy as sc
AXJ = json.load(open(f"{W}/frozen/genesets/stage1_geneset_definitions.json"))
NET = pd.read_csv(f"{W}/frozen/net/collectri_regulons.csv")
nms = [l.strip() for l in open(f"{W}/data/counts/CAP1/counts_unfiltered/cells_x_genes.genes.names.txt")]
before = set(n for n in nms if n and not n.startswith("ENSG"))
V = sc.read_h5ad(f"{W}/data/processed/all_cells_validation.h5ad", backed="r")
after = set(str(x) for x in V.var_names); V.file.close()
MARK = ["COL1A1","COL1A2","DCN","LUM","PDGFRB","PECAM1","VWF","CDH5","RGS5","ACTA2","MYH11",
        "PTPRC","LYZ","CD68","CD3E","MS4A1","MZB1","DERL3","TNFRSF17","KRT14","KRT5","S100B","HBB"]
REC = ["TNFRSF1A","TNFRSF1B","IL1R1","IL1RAP","IL6R","IL6ST","OSMR","LIFR","CCR2","CXCR4"]
sets = [("compartment markers", MARK), ("receptor 10", REC),
        ("RELA target", sorted(set(NET.loc[NET.source=="RELA","target"])))] + \
       [(f"axis {k}", sorted(set(v))) for k, v in AXJ["axes"].items()]
rows = []
for lab, gs in sets:
    b = sum(g in before for g in gs); a = sum(g in after for g in gs)
    rows.append(dict(gene_set=lab, n_def=len(gs), found_before=b, found_after=a,
                     pct_before=round(100*b/len(gs),1), pct_after=round(100*a/len(gs),1),
                     recovered=a-b))
C = pd.DataFrame(rows)
TAB = os.environ.get("APSC_TABDIR", f"{W}/results/tables")
C.to_csv(f"{TAB}/20b_geneset_coverage_before_after.csv", index=False)
print(C.to_string(index=False))
print("\n個別確認:", {g: (g in before, g in after) for g in ["PTPRC","MYH11","DERL3"]},
      "（before, after）")
print("検証の遺伝子空間:", len(before), "->", len(after))
require_heartbeat("20b 完了時")
