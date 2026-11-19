import sys, os
from pathlib import Path
sys.path.insert(0, os.getcwd())
import pandas as pd
import matplotlib; matplotlib.use("Agg")
from shared import data_access
from shared.dataset_selection import filtered_masked_distance_path, _merge_companion_extra
from shared import nav15_figure5

root = Path(".")
D = Path(os.environ["VGIC_DATA_ROOT"])
specs = {
 "WT | vanilla":  (D/"nav15/dataDistances/26-07-27_Nav15_wt_vanillaAF2_distances_all_ok_rmsd_3A.csv", "WT"),
 "WT | masked":   (data_access.resolve(filtered_masked_distance_path(root,"nav15","WT","v1")), "WT"),
 "QQQ | vanilla": (D/"nav15/dataDistances/26-07-27_Nav15_qqq_vanilla_AF2_distances_all_ok_rmsd_3A.csv", "QQQ"),
 "QQQ | masked":  (data_access.resolve(filtered_masked_distance_path(root,"nav15","QQQ","v1")), "QQQ"),
}
def read(p):
    p=Path(p); f=pd.read_csv(p, low_memory=False)
    return _merge_companion_extra(p, f)
datasets={k:(read(v[0]), v[1]) for k,v in specs.items()}
for k,(f,_) in datasets.items(): print(f"  {k:<16} rows={len(f)}")

out = root/"docs/figures/nav15_figure6"; out.mkdir(parents=True, exist_ok=True)

s1, fig1 = nav15_figure5.plot_gate_engagement_coupling(datasets)
fig1.savefig(out/"Nav15_Figure5_gate_engagement_coupling.png", dpi=300, bbox_inches="tight")
fig1.savefig(out/"Nav15_Figure5_gate_engagement_coupling.pdf", bbox_inches="tight")
print("\ngate-engagement coupling summary:"); print(s1.round(3).to_string(index=False))

exp = None
try:
    from shared.nav15_figure5 import experimental_gate_shape
    exp,_ = experimental_gate_shape(root/"nav15"/"experimental")
except Exception as e:
    print("experimental_gate_shape unavailable:", type(e).__name__, str(e)[:90])
if exp is not None:
    s2, fig2 = nav15_figure5.plot_orthogonal_gate_geometry(datasets, exp)
    fig2.savefig(out/"Nav15_Figure5_orthogonal_gate_geometry.png", dpi=300, bbox_inches="tight")
    fig2.savefig(out/"Nav15_Figure5_orthogonal_gate_geometry.pdf", bbox_inches="tight")
    print("\northogonal gate summary:"); print(s2.round(3).to_string(index=False))
