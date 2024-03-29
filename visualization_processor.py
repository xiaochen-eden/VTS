# visualization_processor.py
import tkinter as tk
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import textwrap
import tkinter.filedialog as fd
import PyWGCNA
from tkinter import messagebox
import numpy as np


def save_plot(fig):
    filetypes = (
        ('PNG files', '*.png'),
        ('PDF files', '*.pdf'),
        ('SVG files', '*.svg'),
        ('All files', '*.*')
    )

    filepath = fd.asksaveasfilename(
        title='Save as...',
        initialdir='/',
        filetypes=filetypes,
        defaultextension='.png'
    )

    if not filepath:
        return
    fig.savefig(filepath)

def visualize_plot(plot_type, gene_expression_data, plot_window,plot_period,plot_gene,plot_func):

    if plot_type == "Bar Plot":
        fig_size = (10,6)
    elif plot_type == "Box Plot":
        fig_size = (15,6)
    elif plot_type == "Heatmap":
        fig_size = (15,10)

    tableau_colors = list(mcolors.TABLEAU_COLORS.values())
    colors = tableau_colors[:20]
    fig = Figure(figsize=fig_size, dpi=100)
    ax = fig.add_subplot(111)

    if plot_type == "Bar Plot":
        visualize_bar_plot(gene_expression_data,ax,plot_period,plot_gene,plot_func,colors)
    elif plot_type == "Box Plot":
        visualize_box_plot(gene_expression_data,ax,plot_gene)
    elif plot_type == "Heatmap":
        visualize_heatmap(gene_expression_data,ax,colors)

    canvas = FigureCanvasTkAgg(fig, master=plot_window)
    canvas.get_tk_widget().pack()
    canvas.draw()
    fig.tight_layout()

    save_button = tk.Button(plot_window, text="Save Plot", command=lambda: save_plot(fig))
    save_button.pack(side=tk.BOTTOM,anchor=tk.CENTER,expand=True)


def auto_set_metadata_color(df):
    color_pool = ['green', 'yellow', 'darkviolet', 'deeppink', 'thistle', 'plum', 'violet', 'purple', 'red', 'blue']
    metadata_color_mapping = {}

    for column in df.columns:
        unique_values = df[column].unique()
        value_color_map = {}
        for i, value in enumerate(unique_values):
            color = color_pool[i % len(color_pool)]
            value_color_map[value] = color
        metadata_color_mapping[column] = value_color_map
    return metadata_color_mapping
def WGCNA_plot(WGCNA_adata,usr_par):
    messagebox.showinfo("Upload","Data processing in progress, please wait ......")
    try:
        TPM_cutoff_value = int(usr_par.get("TPM cutoff","1"))
    except ValueError:
        TPM_cutoff_value = 1
    TOMType = usr_par.get("TOM Type", 'signed')
    if not TOMType or TOMType not in ['signed', 'unsigned']:
        TOMType = 'signed'
    print(usr_par)
    pyWGCNA_gene_exp = PyWGCNA.WGCNA(
        species=usr_par.get("Species", 'sugarcane'),
        outputPath=usr_par.get("Output Path", './WGCNA_out/'),
        anndata=WGCNA_adata,
        TPMcutoff= TPM_cutoff_value,
        TOMType=TOMType,
        networkType=usr_par.get("Network Type", 'signed hybrid'),
        save=True,
        figureType='pdf')
    pyWGCNA_gene_exp.geneExpr.to_df().head(5)
    pyWGCNA_gene_exp.preprocess()
    pyWGCNA_gene_exp.findModules()
    metadata_color_mapping = auto_set_metadata_color(WGCNA_adata.obs)
    for column, mapping in metadata_color_mapping.items():
        pyWGCNA_gene_exp.setMetadataColor(column, mapping)
    pyWGCNA_gene_exp.analyseWGCNA()
    pyWGCNA_gene_exp.saveWGCNA()
    messagebox.showinfo("DONE","Data analysis completed!")

def visualize_bar_plot(gene_expression_data, ax,plot_period,plot_gene,plot_func,colors):
    bar_data = gene_expression_data[
        (gene_expression_data['Period'] == plot_period) &
        (gene_expression_data['Gene'] == plot_gene) &
        (gene_expression_data['Functions'] == plot_func)
        ]
    Organ = bar_data['Organ']
    Expression = bar_data['Expression']
    ax.bar(Organ,Expression,color=colors)
    ax.set_title('Bar Plot of {} Expression Across {} Tissues'.format(plot_gene, plot_period))
    ax.set_xlabel('Organ')
    ax.set_ylabel('Expression Level')

def visualize_box_plot(gene_expression_data, ax,plot_gene):
    box_data = gene_expression_data[
        (gene_expression_data['Gene'] == plot_gene)
        ]
    wrapped_labels = ['\n'.join(textwrap.wrap(label, 22)) for label in box_data['Period'].unique()]
    grouped_data = [group["Expression"].values for name, group in box_data.groupby("Period")]
    ax.boxplot(grouped_data, patch_artist=True)
    #violin_data.boxplot(column='Expression',by='Period',ax=ax,patch_artist=True)
    ax.set_title('Box Plot of {} Expression Across all Tissues'.format(plot_gene))
    ax.grid(False)
    ax.set_xticklabels(wrapped_labels)
    ax.set_ylabel('Expression Level')

def visualize_heatmap(gene_expression_data, ax,colors):
    gene_expression_sum = gene_expression_data.groupby('Gene')['Expression'].sum()
    top20_genes = gene_expression_sum.nlargest(20).index
    top20_data = gene_expression_data[gene_expression_data['Gene'].isin(top20_genes)]
    pivot_table = top20_data.pivot_table(index='Gene', columns='Organ', values='Expression', fill_value=0)

    genes = pivot_table.index.tolist()
    organs = pivot_table.columns.tolist()

    c = ax.imshow(pivot_table, cmap='coolwarm', aspect='auto')
    plt.colorbar(c, ax=ax)
    ax.set_xticks(np.arange(len(organs)))
    ax.set_yticks(np.arange(len(genes)))
    ax.set_xticklabels(organs)
    ax.set_yticklabels([textwrap.fill(gene, 22) for gene in genes])
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")

    ax.set_title('Heatmap of Top 20 Genes by Total Expression Across All Organs')
    ax.set_xlabel('Organ')
    ax.set_ylabel('Gene')
