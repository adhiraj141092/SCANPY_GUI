# SCANPY_GUI
A GUI based frameware for single cell RNA-Seq data analysis using scanpy. This tool runs on Python >= 3.7<br>
It takes only 10x Genomics h5 files for the analysis. An example of the file is given in data folder which contains single cell RNA-Seq data of mouse model as an example sequence. 

<h3>Installation:</h3>

<ol>1. Download and extract the content of the folder.</ol>
<ol>2. Using terminal, navigate to the extracted folder. Eg. cd /path/to/folder. </ol>
<ol>3. Install the necessary python packages by instaling the packages in requirements.txt using pip:</ol>
<pre>pip install -r requirements.txt</pre>

<h3>To run the tool</h3>
<ol>Open terminal and execute try.py by entering the following command:</ol>
<pre>python3 try.py</pre>
<ol>Open a browser and enter the following in the address bar:</ol>
<pre>127.0.0.1:5000</pre> or 
<pre>localhost</pre>
<h4>Follow the steps given below to perform the analysis</h4>
<br>Step 1: Upload the file in .h5 format<br>

![upload1](https://user-images.githubusercontent.com/102503979/229825745-30f9847c-ac6c-4494-81f1-144d4427927c.png)

<br>Step 2: Pre-process the data. Select the minimum number of cells a gene is expressed and minimum genes to be expressed in a cell. Enter the mitochondrial gene annotation symbol<br>

![upload2](https://user-images.githubusercontent.com/102503979/229825757-ffa12510-6285-4c99-ae13-6f35424cd982.png)
<br>Step 3: Visualize the data and further preprocess the data by entering the maximum transcripts and maximum pct mitochondrial genes. Next:
<li>Normalise the data. </li>
<li>Perform Logarithmic Transformation</li>
<li>Calculate the highly variable genes</li>
<br>

![upload3](https://user-images.githubusercontent.com/102503979/229825766-3412b985-e399-429e-b26c-125e00802e5f.png)

<br>Step 4: Principal Component Analysis. Select the optimal number of principal components from the scree plot. Visualize the PC1 vs PC2 plot for the genes.<br>
<br>![upload4](https://user-images.githubusercontent.com/102503979/229825770-5b70e081-e5fe-4771-bfc3-7f9d681dd2a7.png)

<br>Step 5: Uniform Manifold Approximation and Projection (UMAP). A default clustering is done using leiden algorithm. Visualize the UMAP1 vs UMAP2 pot for the genes. 
<br>![upload5](https://user-images.githubusercontent.com/102503979/229825776-97a6c791-4028-4e0d-8efe-f3cdd72d28ec.png)


<br>
Created By:<br>
<a href = "https://iitg.ac.in/stud/adhiraj/">Adhiraj Nath</a><br>
PhD Candidate, IIT Guwahati<br>
