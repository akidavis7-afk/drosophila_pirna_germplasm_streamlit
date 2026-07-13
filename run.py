from pathlib import Path
import argparse
from pirna_germplasm_qc.pipeline import run_from_paths
from pirna_germplasm_qc.synthetic import generate_demo_data

def build_parser():
    parser=argparse.ArgumentParser(description="Integrate Drosophila piRNA, Aub imaging, and germ-cell phenotype tables.")
    sub=parser.add_subparsers(dest="command",required=True)
    g=sub.add_parser("generate-demo"); g.add_argument("--output-dir",type=Path,default=Path("data"))
    a=sub.add_parser("analyze")
    a.add_argument("--metadata",type=Path,required=True); a.add_argument("--pirna",type=Path,required=True)
    a.add_argument("--transposons",type=Path,required=True); a.add_argument("--imaging",type=Path,required=True)
    a.add_argument("--germ-cells",type=Path,required=True); a.add_argument("--config",type=Path,required=True)
    a.add_argument("--output-dir",type=Path,required=True)
    d=sub.add_parser("demo"); d.add_argument("--data-dir",type=Path,default=Path("data")); d.add_argument("--config",type=Path,default=Path("configs/demo.yaml")); d.add_argument("--output-dir",type=Path,default=Path("results/demo"))
    return parser

def main():
    args=build_parser().parse_args()
    if args.command=="generate-demo":
        generate_demo_data(args.output_dir); print(f"Synthetic demonstration written to {args.output_dir}"); return
    if args.command=="analyze":
        run_from_paths(args.metadata,args.pirna,args.transposons,args.imaging,args.germ_cells,args.config,args.output_dir); print(f"Analysis complete: {args.output_dir}"); return
    if args.command=="demo":
        generate_demo_data(args.data_dir)
        run_from_paths(args.data_dir/"sample_metadata.csv",args.data_dir/"pirna_counts.csv",args.data_dir/"transposon_expression.csv",args.data_dir/"imaging_summary.csv",args.data_dir/"germ_cell_counts.csv",args.config,args.output_dir)
        print(f"Demo complete: {args.output_dir}")
if __name__=="__main__": main()
