#!/usr/bin/env python3
"""Aggregate native and cross-task chromosome-learning metrics."""
import argparse, csv, json
from pathlib import Path
from prepare_GRCh38_chromosome_learning_results import CHROMS, binary_metrics, read_scores

def table(path):
    if not path.is_file() or not path.stat().st_size: raise FileNotFoundError(path)
    with path.open(newline="") as h: return list(csv.DictReader(h,delimiter="\t"))

def main():
    p=argparse.ArgumentParser(); p.add_argument("--repository",type=Path,default=Path(__file__).resolve().parents[1]); p.add_argument("--background-type",choices=("real","synthetic"),required=True); p.add_argument("--models",nargs="+",choices=("lsgkm","alphagenome"),default=("lsgkm","alphagenome")); p.add_argument("--native-only",action="store_true"); p.add_argument("--output",type=Path); a=p.parse_args()
    if not a.native_only and a.background_type != "real": p.error("cross-task metrics are currently restricted to real backgrounds")
    root=a.repository.resolve(); default_output=(root/f"results/comparison/GRCh38/{a.background_type}/native_metrics.tsv" if a.native_only else root/f"results/cross_task/GRCh38/{a.background_type}/native_vs_cross_metrics.tsv"); output=a.output or default_output; rows=[]
    for model,folder in (("lsgkm","ls-gkm"),("alphagenome","alphagenome")):
      if model not in a.models: continue
      for mt,et in (("evn","svn"),("svn","evn")):
        native=table(root/f"results/{mt}/GRCh38/{folder}/{a.background_type}/chromosome_learning_curves/metrics.tsv")
        lookup={(r["checkpoint"],r["split"]):r for r in native}
        for idx,chrom in enumerate(CHROMS,1):
          name=f"through_chr{chrom}"
          for split in ("train_100","validation","test"):
            n=lookup[(name,split)]
            common=dict(model=model,model_task=mt,negative_kind=a.background_type,checkpoint=name,checkpoint_index=idx,endpoint_chromosome=f"chr{chrom}",n_training_total=int(n["n_training_total"]),split=split)
            rows.append(dict(**common,evaluation_task=mt,evaluation_mode="native",n_evaluation_total=int(n["n_evaluation_total"]),**{k:float(n[k]) for k in ("auroc","auprc","accuracy","precision","recall","f1")}))
            if a.native_only: continue
            base=root/f"predictions/cross_task/GRCh38/{folder}/{a.background_type}/{mt}_model_on_{et}/chromosome_learning_curves/{name}"
            if model=="lsgkm":
              pos=read_scores(base/f"{split}_pos_scores.txt"); neg=read_scores(base/f"{split}_neg_scores.txt"); m=binary_metrics(pos,neg); count=len(pos)+len(neg)
            else:
              path=base/f"{split}_metrics.json"
              if not path.is_file() or not path.stat().st_size: raise FileNotFoundError(path)
              with path.open() as h: payload=json.load(h)
              m=payload["metrics"]; count=int(payload["n_positive"])+int(payload["n_negative"])
            rows.append(dict(**common,evaluation_task=et,evaluation_mode="cross",n_evaluation_total=count,**{k:float(m[k]) for k in ("auroc","auprc","accuracy","precision","recall","f1")}))
    output.parent.mkdir(parents=True,exist_ok=True); tmp=Path(str(output)+".tmp")
    with tmp.open("w",newline="") as h: w=csv.DictWriter(h,fieldnames=list(rows[0]),delimiter="\t",lineterminator="\n"); w.writeheader(); w.writerows(rows)
    tmp.replace(output); print(f"Wrote {len(rows)} rows to {output}")
if __name__=="__main__": main()
