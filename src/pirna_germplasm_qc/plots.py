from pathlib import Path
import matplotlib.pyplot as plt
import pandas as pd

def plot_length(length_df, path: Path):
    fig,ax=plt.subplots(figsize=(8,5))
    for genotype,g in length_df.groupby("genotype"):
        g=g.sort_values("length")
        ax.plot(g["length"],g["fraction"],marker="o",label=genotype)
    ax.set_xlabel("piRNA length (nt)"); ax.set_ylabel("Fraction of reads"); ax.set_title("piRNA length distribution")
    ax.legend(); ax.grid(alpha=.25); fig.tight_layout(); fig.savefig(path,dpi=180); plt.close(fig)

def plot_family_heatmap(family_df, path: Path):
    matrix=family_df.pivot_table(index="transposon_family",columns="genotype",values="mean_cpm",aggfunc="mean").fillna(0)
    fig,ax=plt.subplots(figsize=(8,5))
    im=ax.imshow(matrix.to_numpy(),aspect="auto")
    ax.set_xticks(range(len(matrix.columns))); ax.set_xticklabels(matrix.columns,rotation=30,ha="right")
    ax.set_yticks(range(len(matrix.index))); ax.set_yticklabels(matrix.index)
    ax.set_title("Mean piRNA CPM by transposon family")
    fig.colorbar(im,ax=ax,label="Mean CPM"); fig.tight_layout(); fig.savefig(path,dpi=180); plt.close(fig)

def plot_bar(df, x, y, title, ylabel, path: Path):
    fig,ax=plt.subplots(figsize=(7,5))
    ax.bar(df[x],df[y])
    ax.set_ylabel(ylabel); ax.set_title(title); fig.tight_layout(); fig.savefig(path,dpi=180); plt.close(fig)

def plot_integrated(df, path: Path):
    fig,ax=plt.subplots(figsize=(8,6))
    sizes=35+8*df["mean_germ_cells"].fillna(0)
    ax.scatter(df["mean_pingpong_score"],df["mean_aub_posterior_enrichment"],s=sizes,alpha=.75)
    for _,r in df.iterrows():
        ax.annotate(r["genotype"],(r["mean_pingpong_score"],r["mean_aub_posterior_enrichment"]),xytext=(5,5),textcoords="offset points")
    ax.set_xlabel("Mean synthetic ping-pong score"); ax.set_ylabel("Mean Aub posterior enrichment")
    ax.set_title("Integrated piRNA → Aub localization → germ-cell summary")
    ax.grid(alpha=.25); fig.tight_layout(); fig.savefig(path,dpi=180); plt.close(fig)
