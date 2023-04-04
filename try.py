from flask import Flask, request, jsonify, current_app, render_template, Response, session, redirect, url_for
import scanpy as sc
from werkzeug.utils import secure_filename
import asyncio
import matplotlib.pyplot as plt
import scanpy as sc
import pandas as pd
import numpy as np
import os
import time
import random
import tempfile
import copy
import threading
import matplotlib
matplotlib.use('Agg')
import seaborn as sns


app = Flask(__name__)
app.config['UPLOAD_FOLDER'] = 'static/'
app.config['ALLOWED_EXTENSIONS'] = {'h5', 'h5ad'}
app.config['SECRET_KEY'] = '123zxc'

adata = None
copy_adata = None
genes = []

def allowed_file(filename):
    return '.' in filename and filename.rsplit('.', 1)[1].lower() in app.config['ALLOWED_EXTENSIONS']

async def load_async(file):
    # Your asynchronous code here
    loop = asyncio.get_event_loop()
    # adata = await loop.run_in_executor(None, sc.read_10x_h5, file)
    adata = await loop.run_in_executor(None, sc.read_10x_h5, file)
    os.remove(file)
    # Just an example of an async task that takes 1 second
    return adata



@app.route("/")
def index():
    return render_template("index2.html")

@app.route('/getdata', methods=['POST'])
async def success():
    if request.method == 'POST':
        f = request.files['file']
        f.save(f.filename)
        global adata
        adata = await load_async(f.filename)
        dat = adata.var.head
        fig, ax = plt.subplots(figsize=(4, 3))

        # call sc.pl.highest_expr_genes() to generate the plot and assign the return value to a variable
        highest_expr_genes_plot = sc.pl.highest_expr_genes(adata, n_top=20, show=False, ax=ax)

        # save the figure as a PNG file with a defined size
        fig.savefig('static/highest_expr_genes.png', dpi=100, bbox_inches='tight')
        plt.close()
        img = os.path.join(app.config['UPLOAD_FOLDER'], 'highest_expr_genes.png')

        n_cells, n_genes = adata.shape

        return render_template("upload.html", data = dat, plot=img, det_a = adata, cell = n_cells, genes = n_genes)

@app.route('/prep', methods=['POST'])
def preprocess():
    global adata
    adata.var_names_make_unique()
    gene = request.form['mn_gene']
    cell = request.form['cell']
    mt = request.form['mt']
    num_highly_variable = None
    num_not_highly_variable = None



    #Preprocessing
    adata.var_names_make_unique()
    sc.pp.filter_cells(adata, min_genes=int(gene))
    sc.pp.filter_genes(adata, min_cells=int(cell))
    print("ok done1")


    adata.var['mt'] = adata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'

    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

    # Prep 2
    if 'n_genes' in request.form and 'pct_mt' in request.form:
        n_genes = request.form['n_genes']
        pct_mt = request.form['pct_mt']
        adata = adata[adata.obs.n_genes_by_counts < int(n_genes), :]
        print(cell)
        adata = adata[adata.obs.pct_counts_mt < int(pct_mt), :]
        print("ok done 2")





    fig2, axs = plt.subplots(ncols=3, figsize=(10, 4), sharey=True)

    sns.violinplot(data=adata.obs, y='n_genes_by_counts', ax=axs[0], color='blue')
    sns.violinplot(data=adata.obs, y='total_counts', ax=axs[1], color='green')
    sns.violinplot(data=adata.obs, y='pct_counts_mt', ax=axs[2], color='orange')
    plt.tight_layout()
    plt.close()

    fig2.savefig('static/qc.png', dpi=150)
    img = os.path.join(app.config['UPLOAD_FOLDER'], 'qc.png')

    pltn = plots(adata)
    n_cells, n_genes = adata.shape
    global genes
    genes = adata.var_names.tolist()
    return render_template("prep.html", det_a=adata, plot=img, plt1=pltn, cell=n_cells, genes=n_genes)

@app.route('/analysis', methods=['POST'])
def analysis():
    global adata
    if 'norm' in request.form:
        norm = request.form['norm']
        sc.pp.normalize_total(adata, target_sum=int(norm))

    if 'log' in request.form:
        sc.pp.log1p(adata)


    min_mean = request.form['min_mean']
    max_mean = request.form['max_mean']
    min_disp = request.form['min_disp']
    sc.pp.highly_variable_genes(adata, min_mean=float(min_mean), max_mean=float(max_mean), min_disp=float(min_disp))
    value_counts = adata.var['highly_variable'].value_counts()
    num_highly_variable = value_counts[True]
    num_not_highly_variable = value_counts[False]

    n_cells, n_genes = adata.shape

    return render_template("analysis.html", num_highly_variable = num_highly_variable, cell=n_cells, genes=n_genes)

@app.route('/analysis2', methods=['POST','GET'])
def analysis2():
    global adata
    genes = adata.var_names.tolist()


    img2 = None
    img = None

    sc.tl.pca(adata)
    var_ratio = adata.uns['pca']['variance_ratio']
    n_components = len(var_ratio)

    fig, ax = plt.subplots(figsize=(4, 3))
    ax.plot(range(1, n_components + 1), var_ratio, 'o-')
    ax.set_xlabel('PC')
    ax.set_ylabel('Variance ratio')
    ax.set_xticks(range(1, n_components + 1))

    fig.savefig('static/pca_plot2.png', dpi=150)
    img2 = os.path.join(app.config['UPLOAD_FOLDER'], 'pca_plot2.png')

    selected_gene = None
    selected_gene = request.form.get('gene')

    if 'gene' in request.form:
        # Perform the PCA plot using the selected gene
        sc.pp.pca(adata)  # Perform PCA to reduce the number of dimensions
        fig1, ax = plt.subplots(figsize=(4, 3))
        pca_plot = sc.pl.pca(adata, color=selected_gene ,show=False, ax=ax)
        fig1.savefig('static/pca_plot.png', dpi=150)
        img = os.path.join(app.config['UPLOAD_FOLDER'], 'pca_plot.png')

        return render_template("analysis2.html", genes=genes, img=img, img2=img2)

    # if 'gene' in request.form:
    #     # Perform the PCA plot using the selected gene
    #     sc.pp.pca(adata)  # Perform PCA to reduce the number of dimensions
    #     pca_df = adata.obsm['X_pca'][:, :3]
    #     pca_df = pd.DataFrame(data=pca_df, columns=['PC1', 'PC2', 'PC3'])
    #     pca_df[selected_gene] = adata.obs[selected_gene].values
    #
    #     import plotly.express as px
    #
    #     # Assuming your dataframe is called `df`, you can extract the first three columns like this:
    #     # df_subset = df_original.iloc[:, :3]
    #
    #
    #
    #     fig = px.scatter_3d(pca_df, x='PC1', y='PC2', z='PC3', color=selected_gene)
    #     fig.update_traces(marker_size=5)
    #     plot_div = fig.to_html(full_html=False)
    #
    #     # Render the template with the interactive plot
    #     return render_template("analysis2.html", genes=genes, plot_div=plot_div, img2=img2)


        # sc.pl.pca_variance_ratio(adata, log=True)
        #
        #


    # adata = preprocess_data(adata)
    global copy_adata
    copy_adata = adata
    return render_template("analysis2.html", genes=genes, img = img, img2 = img2)

@app.route("/prep2", methods=['POST','GET'])
def preprocess_data():
    # Perform preprocessing steps on adata

    global adata
    n_pcgs = None
    n_pcgs = request.form.get('n_pcgs')

    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=int(n_pcgs))
    sc.tl.pca(adata)
    # sc.pp.neighbors(adata)
    sc.tl.leiden(adata)
    sc.tl.paga(adata)
    sc.tl.umap(adata)
    sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')
    sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)
    return redirect('/analysis3')


@app.route('/analysis3', methods=['POST','GET'])
def analysis3():
    global adata
    global copy_adata
    global genes
    genes1 = copy_adata.var_names.tolist()
    print(len(genes1))


    sc.pl.paga(adata, show=False)

    gene1 = None
    if 'gene1' in request.form:
        gene1 = request.form['gene1']


    if gene1 != None:
        plt.subplots(figsize=(8, 6))
        sc.pl.umap(adata, color=gene1, use_raw=False)
        plt.savefig('static/umap.png', dpi=300, bbox_inches='tight')
        plt.close()


    else:
        plt.subplots(figsize=(8, 6))
        sc.pl.umap(adata, color=['leiden'], use_raw=False)
        plt.savefig('static/umap.png', dpi=300, bbox_inches='tight')
        plt.close()

    img = os.path.join(app.config['UPLOAD_FOLDER'], 'umap.png')

    if 'file_name' in request.form:
        result_file = request.form['file_name']
        adata.write(result_file)
        return redirect('analysis3')

    return render_template("analysis3.html", genes=genes1, img = img)




 # if 'nbr' in request.form:
    #     sc.pp.neighbors(adata)  # Calculate the neighborhood graph of cells based on PCA components
    # if 'umap' in request.form:
    #     sc.tl.umap(adata)
    # if 'loalgo' in request.form:
    #     sc.tl.louvain(adata)


# def pca_plot():


  # Perform UMAP for visualization
# def preprocess_data(adata):
#     # Perform preprocessing steps on adata
#     sc.tl.pca(adata)
#     sc.pp.neighbors(adata)
#     sc.tl.leiden(adata)
#     sc.tl.paga(adata)
#     sc.tl.umap(adata)
#     return adata

def plots(adata):
    fig, axs = plt.subplots(ncols=2, figsize=(10, 4), sharey=True)

    # plot scatterplot on the left subplot
    sns.scatterplot(data=adata.obs, x='total_counts', y='n_genes_by_counts', color='blue', ax=axs[0])
    axs[0].set_xlabel('Total Counts')
    axs[0].set_ylabel('Number of Genes by Counts')

    # plot scatterplot on the right subplot
    sns.scatterplot(data=adata.obs, x='total_counts', y='pct_counts_mt', color='red', ax=axs[1])
    axs[1].set_xlabel('Total Counts')
    axs[1].set_ylabel('Percentage of Counts MT', labelpad=10)
    axs[1].yaxis.set_label_coords(-0.1, 0.5)

    fig.tight_layout()  # adjust spacing between subplots to fit labels
    fig.savefig('static/tot_vs_n_genes.png', dpi=150)
    img = os.path.join(app.config['UPLOAD_FOLDER'], 'tot_vs_n_genes.png')
    plt.close()
    return img

# def success_backup():
#     if request.method == 'POST':
#         f = request.files['file']
#         f.save(f.filename)
#         adata = sc.read_10x_h5(f.filename)
#         dat = adata.var.head
#         command = f.filename
#         os.remove(command)
#         return render_template("upload.html", data = dat)

@app.route('/about')
def about():
    return render_template('about.html')



@app.after_request
def add_header(response):
    # response.cache_control.no_store = True
    response.headers['Cache-Control'] = 'no-store, no-cache, must-revalidate, post-check=0, pre-check=0, max-age=0'
    response.headers['Pragma'] = 'no-cache'
    response.headers['Expires'] = '-1'
    return response


if __name__ == '__main__':
    app.run(debug=True)
