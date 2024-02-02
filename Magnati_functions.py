import subprocess
import importlib

packages_error = []

def check_installation(package):
    try:
        importlib.import_module(package)
        return True
    except ImportError:
        print(f"The package '{package}' is not installed.")
        install = input(f"Do you want to install the '{package}' package? (yes/no): ").lower()
        if install == 'yes' or install == 'y':
            try:
                subprocess.check_call(['pip', 'install', package])
                print(f"The package '{package}' has been successfully installed.")
                return True
            except subprocess.CalledProcessError as e:
                print(f"Error occurred while installing '{package}': {e}")
                packages_error.append(package)
                return False
        else:
            print(f"'{package}' package has not been installed.")
            return False


# Check installation of each required package
packages = [
    'pandas', 'numpy', 'joblib', 'matplotlib', 'scipy', 'statsmodels', 'itertools', 'warnings','pingouin',
    'bootstrapped', 'plotly', 'seaborn', 'scikit-learn', 'kaleido', 'dash-bio'
]

all_packages_installed = all(check_installation(package) for package in packages)

try:
    import os
    import pandas as pd
    import numpy as np
    import joblib
    import warnings
    warnings.filterwarnings("ignore")

    # Stats
    import scipy.stats as stats
    import statsmodels.stats.multitest as smt
    from statsmodels.stats.outliers_influence import variance_inflation_factor
    import math
    import statsmodels.api as sm
    from statsmodels.stats.diagnostic import het_breuschpagan
    from scipy.stats import levene, kruskal, mannwhitneyu, shapiro, zscore, norm, ttest_ind
    import pingouin as pg #Post Hoc Games-Howell
    import bootstrapped.bootstrap as bs
    import bootstrapped.compare_functions as bs_compare
    import bootstrapped.stats_functions as bs_stats
    from statsmodels.tsa.seasonal import seasonal_decompose
    from scipy.stats import randint, uniform
    from statsmodels.tools.tools import add_constant
    from statsmodels.stats.multitest import multipletests
    from itertools import combinations

    # Data Visaulisation
    import matplotlib.pyplot as plt
    import seaborn as sns
    import scipy
    import plotly.express as px
    import plotly.graph_objects as go
    import plotly.io as pio
    import kaleido
    
    from sklearn.preprocessing import StandardScaler, MinMaxScaler
    from sklearn.cluster import KMeans, SpectralClustering, MiniBatchKMeans,DBSCAN 
    from sklearn.manifold import TSNE
    from sklearn.utils.multiclass import unique_labels
    from sklearn.decomposition import PCA
    from sklearn.metrics import silhouette_samples, silhouette_score
    
    #dash bio
    import dash_bio


except ImportError as e:
    print(f"Error importing package: {e}")
    all_packages_installed = False

if all_packages_installed:
    print("All required packages have been imported successfully. You can proceed with your code.")
    
else:
    if packages_error:
        print("\nThe following packages have not been installed or require an upgrade:")
        print(', '.join(packages_error))
        print("Please install them and try again. You can use 'pip install <package_name>'. It might require an internet connection.")
        if 'dash-bio' in packages_error:
            print("\nIf the error pertains to the Dash Bio package, it could indicate a distinct issue. Please refer to the documentation at https://dash.plotly.com/dash-bio and attempt to install it directly on the kernel. Dash Bio is utilized for generating specific graphs, such as volcano plots. If installation proves unsuccessful, you can still proceed with the remainder of the code.")



        
        
def plot_distribution_qqplot(df, test, save_fig=False):
    # Seleziona solo le colonne numeriche
    numerical = df.select_dtypes(include='number')
    normal_genes = []

    for feature in numerical:
        data = df[feature]
        p_value = perform_normal_test(data, test)
            
        if p_value > 0.05:
            normal_genes.append(feature)
            # Plot of data's distribution
            plt.figure(figsize=(12, 4))
            plt.subplot(1, 3, 1)
            sns.histplot(data, kde=True)
            plt.xlabel(feature)
            plt.title("Distribution of Data")

            if save_fig:
                plt.savefig(f'{feature}_Distribution of Data.png')

            # QQ plot
            plt.subplot(1, 3, 2)
            stats.probplot(data, dist="norm", plot=plt)
            plt.xlabel("Theoretical Quantiles")
            plt.title("QQ plot for " + feature)

            plt.tight_layout()

            if save_fig:
                plt.savefig(f'{feature}_qqplot.png')
                
            plt.subplot(1, 3, 3)
            sns.boxplot(data)
            plt.xlabel("Box Plot")
            plt.title("Boxplot for " + feature)
            plt.tight_layout()
            
            if save_fig:
                plt.savefig(f'{feature}_boxplot.png')


    if len(normal_genes) == 0:
        print("From Shapiro Test: I'm confident to 95% to reject the null hypothesis: No feature follows a normal distribution")

    return normal_genes




def grubbs_test_indices(data):
    z_scores = zscore(data)
    n = len(data)
    critical_value = norm.ppf(1 - 0.05 / (2 * n))
    test_statistic = max(abs(z_scores))
    return [i for i, z in enumerate(z_scores) if abs(z) > critical_value]

def find_outliers(df):
    df = df.select_dtypes(include='number')
    outliers_indices = {}
    for col in df.columns:
        outliers = grubbs_test_indices(df[col])
        if outliers:
            outliers_indices[col] = outliers
            print(f"Outliers in column {col}: {outliers}")
    return outliers_indices

def corr_matrix_(df, method_, title, threshold, show_plot, save_fig=False, annot=False, df_target_genes=None):
    numerical = df.select_dtypes(include='number')
    
    corr = numerical.corr(method=method_)  # Calcola la matrice di correlazione
    
    if df_target_genes is not None:
        col_target = list(col for col in df_target_genes)

        significant_genes = []
        for col in col_target:
            significant_genes.extend(corr[corr[col].abs() > threshold].index)

        significant_genes = list(set(significant_genes) - set(col_target))
        corr = corr.loc[significant_genes, col_target]
    
    # Crea e restituisci un DataFrame con i valori di correlazione significativi
    df_corr_values = pd.DataFrame(corr[(corr.abs() > threshold)].stack(), columns=['Correlation'])
    
    if show_plot:
     
        fig, ax = plt.subplots(figsize=(12, 6))
        if annot == True:
            sns.heatmap(corr, annot=True, cmap='mako_r', annot_kws={'size': 6, 'color': 'white'})
        else:
            sns.heatmap(corr, annot=False, cmap='mako_r')

        plt.title(title, y=1.0)
        if save_fig:
            plt.savefig(f'{title}.png')  # Salva la figura prima di mostrarla
        plt.show()  # Visualizza la heatmap per il gruppo corrente

    
    return df_corr_values


def perform_VIF(df):
    numerical = df.select_dtypes(include='number')
    # Calcola il VIF per ogni variabile indipendente (escludendo l'intercetta)
    vif = pd.DataFrame()
    vif["Variable"] = numerical.columns
    vif["VIF"] = [variance_inflation_factor(numerical.values, i) for i in range(numerical.shape[1])]
    vif = vif.sort_values(by='VIF', ascending=False)
    # Visualizza il risultato
    print(vif)



def perform_levene(df_long, diff_col):

    # Otteniamo l'elenco dei geni unici nel DataFrame
    gene_columns = df_long['gene'].unique()

    results_list = []  # Lista per memorizzare i risultati

    # Iteriamo su ciascun gene per eseguire il test di Levene
    for gene in gene_columns:
        # Selezioniamo solo i dati del gene corrente
        gene_data = df_long[df_long['gene'] == gene]['expr']

        # Raggruppiamo i dati per fasce di età e otteniamo le liste di dati per ciascun gruppo
        grouped_data = [group for _, group in gene_data.groupby(df_long[diff_col])]

        # Applichiamo il test di Levene ai dati raggruppati
        statistic, p_value = levene(*grouped_data)

        homoschedasticity = 'No' if p_value < 0.05 else 'Yes'

        result = {'Gene': gene,
                  'Statistic': round(statistic, 2),
                  'P-value': round(p_value, 4),
                  'Homoschedasticity': homoschedasticity
                  }

        results_list.append(result)

    # Creiamo un DataFrame dai risultati
    results_df = pd.DataFrame(results_list)


    return results_df




def perform_normal_test(data, test):
    if test == 'Shapiro':
        _, shapiro_p_value = shapiro(data)
        return shapiro_p_value
    
    elif test == "D'Agostino":
        _, p_value = stats.normaltest(data)
        return p_value


def perform_ANOVA(data, gene_, between_):
    anova_result = pg.anova(data=data[data['gene'] == gene_], dv='expr', between=between_, detailed=True)
    return anova_result

def perform_Ttest(df_long, gene, between_, between_group, result_list, calculate_CI):
    group_data = np.array(df_long[(df_long['gene'] == gene) & (df_long[between_] == between_group)]['expr'])
    gene_data = np.array(df_long[df_long['gene'] == gene]['expr'])

    t_statistic, t_p_value = ttest_ind(group_data, gene_data)
    result_dict = None

    # Aggiungi i risultati alla lista solo se il p-value è inferiore a 0.05
    if t_p_value < 0.05:
        result_dict = {
            'Gene': gene,
            'Contrast': f"{between_group}_vs_All",
            'Test': 'T-test',
            'Statistic': round(t_statistic, 2),
            'P-value': round(t_p_value, 4)
        }

        # Aggiungi intervalli di confidenza solo se richiesto
        if calculate_CI:
            confidence_interval = bs.bootstrap_ab(
                group_data,
                gene_data,
                stat_func=bs_stats.mean,
                compare_func=bs_compare.difference,
                alpha=0.05,
                num_iterations=10000
            )

            conf_int_low = round(confidence_interval.lower_bound, 2)
            conf_int_high = round(confidence_interval.upper_bound, 2)

            result_dict.update({
                'CI_low': conf_int_low,
                'CI_high': conf_int_high
            })


    return result_dict



def perform_kruskal(data, gene_, between_):
    gene_data = data[data['gene'] == gene_]
    kruskal_result = kruskal(*[group['expr'].values for name, group in gene_data.groupby(between_)])
    return kruskal_result

def perform_gameshowell(data, gene, between_, games_howell_dict, calculate_CI):
    dati_gene = data[data['gene'] == gene]
    games_howell_result = pg.pairwise_gameshowell(dati_gene, dv='expr', between=between_)
    games_howell_dict[gene] = games_howell_result
    result_dict = None

    for index, row in games_howell_result.iterrows():
        contrast = f"{row['A']}_vs_{row['B']}"
        statistic = round(row['T'], 2)
        p_val = round(row['pval'], 4)
        if p_val < 0.05:
            
            result_dict = {
                'Gene': gene,
                'Contrast': contrast,
                'Test': 'Games-Howell',
                'Statistic': round(statistic, 2),
                'P-value': round(p_val, 4)
                }
            
            
            if calculate_CI:
                group1_data = np.array(data[(data[between_] == row['A'])]['expr'])
                group2_data = np.array(data[(data[between_] == row['B'])]['expr'])
                confidence_interval = bs.bootstrap_ab(
                    group1_data,
                    group2_data,
                    stat_func=bs_stats.mean,
                    compare_func=bs_compare.difference,
                    alpha=0.05,
                    num_iterations=10000
                )

                conf_int_low = round(confidence_interval.lower_bound, 2)
                conf_int_high = round(confidence_interval.upper_bound, 2)
                
                result_dict.update({
                'CI_low': conf_int_low,
                'CI_high': conf_int_high
                })

           

    return result_dict

                
           
def perform_stats_analysis(df, df_long, between_, test, numerical, calculate_CI):
    gene_columns = df_long['gene'].unique()
    results_list = []
    games_howell_test_results = {}
    significant_genes_kruskal = []
    
    normal_genes = []
    
    for feature in numerical:
        data = df[feature]
        p_value = perform_normal_test(data, test)
        if p_value > 0.05:
            normal_genes.append(feature)

    for gene in gene_columns:
        if gene in normal_genes:
            anova_result = perform_ANOVA(df_long, gene, between_)

            if 'p-unc' in anova_result.columns and anova_result['p-unc'][0] < 0.05:
                for between_group in df_long[between_].unique():
                    test_results = perform_Ttest(df_long, gene, between_, between_group, results_list, calculate_CI)
                    if test_results is not None:
                        results_list.append(test_results)
        else:
            kruskal_result = perform_kruskal(df_long, gene, between_)  
            if hasattr(kruskal_result, 'pvalue') and kruskal_result.pvalue < 0.05:
                significant_genes_kruskal.append(gene)

    if len(significant_genes_kruskal) == 0 and len(results_list) == 0:
        print('No significant difference')
        return None

    elif len(significant_genes_kruskal) != 0:
        for gene in significant_genes_kruskal:
            results_GH = perform_gameshowell(df_long, gene, between_, games_howell_test_results, calculate_CI)
            if results_GH is not None:
                results_list.append(results_GH)

    results_df = pd.DataFrame(results_list)

    if 'P-value' not in results_df.columns:
        print("No significant difference")
        return None

    try:
        results_df['Benjamini-Hochberg correction'] = smt.multipletests(results_df['P-value'], method='fdr_bh')[1]
        results_df_filtered = results_df[results_df['Benjamini-Hochberg correction'] < 0.05]
    
    except KeyError as e:
        print(f"KeyError: {e}. Check your datas")
        return None

    

    results_df_filtered['P-value'] = results_df_filtered['P-value'].apply(lambda x: f'{x:.2e}')
    results_df_filtered['Benjamini-Hochberg correction'] = results_df_filtered['Benjamini-Hochberg correction'].apply(lambda x: f'{x:.2e}')
    return results_df_filtered
    
    



def one_feature_analysis(df_original, between_, test, limit, calculate_CI = False):
    print(df_original[between_].value_counts())
    
    df_cat = df_original.loc[:,between_]
    numerical = df_original.select_dtypes(include='number')
    df = pd.concat([df_cat, numerical], axis=1)

    for variable, count in df[between_].value_counts().items():
        if count < limit:
            print(f'Variable {variable} has frequence < {limit}, it was deleted')
            df.drop(df[df[between_] == variable].index, inplace=True)


    
    if len(df[between_].unique()) == 1:
        print(f'There is just one class, impossible plays the test')
        return None
    
    
    df_long = pd.melt(df, id_vars=[between_],
                        var_name='gene', value_name='expr')

    df_results = perform_stats_analysis(df, df_long, between_, test, numerical, calculate_CI)

        
    return df_results



def create_boxplot(df, between, gene_list):
    between_groups = df[between].unique()
    data_to_plot = df[df[between].isin(between_groups)]
    data_to_plot = data_to_plot.melt(id_vars=between, value_vars=gene_list, var_name='Gene', value_name='Expression')
    
    plt.figure(figsize=(5, 4))
    sns.boxplot(data=data_to_plot, x=between, y='Expression', hue='Gene', orient='v')
    plt.title(f'Boxplot of Gene Expression by {between}')
    plt.xlabel(between)
    plt.ylabel('Gene Expression')
    plt.legend(title='Gene')
    plt.xticks(rotation=90)
    plt.show()


    
def create_map(data, genes, between, title, map_, save_fig):

    # Creazione del dataframe per l'età media
    means = data.groupby(between)[genes].mean()

    # Creazione del clustermap
    plt.figure(figsize=(10, 12))
    if map_ == 'heatmap':
        sns.heatmap(means.T, cmap='mako_r', yticklabels=True, annot=False, fmt=".2f", annot_kws={"size": 8})
        plt.xticks(rotation=30)
    elif map_ == 'clustermap': 
        clustermap = sns.clustermap(means.T, cmap='mako_r', method='average', metric='euclidean', 
                       col_cluster=True, row_cluster = True, cbar_pos=(0.02, 0.8, 0.05, 0.15), 
                                    yticklabels=True, annot=False, fmt=".2f", annot_kws={"size": 8})
        
        # Regola la dimensione del font dell'asse y
        #clustermap.ax_heatmap.yaxis.tick_left()
        clustermap.ax_heatmap.yaxis.set_tick_params(labelsize=6)
        clustermap.ax_heatmap.xaxis.set_tick_params(rotation = 45, labelsize=6)

    # Aggiungi titoli e salva il plot
    plt.title(title)
    plt.xlabel(between)
    plt.subplots_adjust(left=0.1, bottom=0.1, right=0.9, top=0.9)
    if save_fig:
        plt.savefig(f"{title}.png", dpi = 300)
    # Visualizza la figura
    plt.show()
    return


def df_scaler(df):
    numerical = df.select_dtypes(include='number')
    scaler = MinMaxScaler()
    scaled_features = scaler.fit_transform(numerical)
    return scaled_features
    



def elbow_method(scaled_features):
    k_values = range(1,11)

    # Lista per memorizzare l'indice di varianza totale
    total_variance = []
    # silhouette score
    score_list = []
    

    # Iterazione su tutti i valori di k sfruttando metodo K-means
    for k in k_values:
        kmeans = KMeans(n_clusters=k, random_state=42)
        kmeans.fit(scaled_features)
        total_variance.append(kmeans.inertia_)
        
        if k!=1:
            cluster_labels = kmeans.labels_
            silhouette_avg = silhouette_score(scaled_features, cluster_labels)
            score_list.append(silhouette_avg)

    # Plot dell'indice di varianza totale
    plt.figure(figsize=(6,4))
    plt.plot(k_values, total_variance, 'bx-')
    plt.xlabel('Numero di Cluster (k)')
    plt.ylabel('Varianza totale')
    plt.title('Metodo Elbow')
    plt.show()
    
    i=2
    for score in score_list:
        print("For n_clusters={0}, the silhouette score is {1}".format(i, score))
        i+=1
        
        
def pca_explained_components(scaled_features):
    pca_var = PCA()
    pca_var.fit(scaled_features)

    # Plot
    plt.figure(figsize=(9,3))
    xi = np.arange(1, 1+scaled_features.shape[1], step=1)
    yi = np.cumsum(pca_var.explained_variance_ratio_)
    plt.plot(xi, yi, marker='o', linestyle='--', color='b')

    # Aesthetics
    plt.ylim(0.0,1.1)
    plt.xlabel('Number of Components')
    plt.xticks(np.arange(1, 1+scaled_features.shape[1], step=1))
    plt.ylabel('Cumulative variance (%)')
    plt.title('Explained variance by each component')
    plt.axhline(y=1, color='r', linestyle='-')
    plt.gca().xaxis.grid(False)



def cluster_analysis(cluster_reduction, scaled_features, df, variable, color_mapping, subplot_index, save_fig=False):
    if cluster_reduction == 'PCA':
        cluster_reduction_ = PCA(n_components=2, random_state=42)
        reduced_features = cluster_reduction_.fit_transform(scaled_features)

        # Gestione del cambiamento di nome del metodo
        if hasattr(cluster_reduction_, 'explained_variance_ratio_'):
            explained_variance_ratio = cluster_reduction_.explained_variance_ratio_
        elif hasattr(cluster_reduction_, 'explained_variance_ratio'):
            explained_variance_ratio = cluster_reduction_.explained_variance_ratio
        else:
            explained_variance_ratio = None
    elif cluster_reduction == 'T-SNE':
        cluster_reduction_ = TSNE(n_components=2, random_state=42)
        reduced_features = cluster_reduction_.fit_transform(scaled_features)
        explained_variance_ratio = None
        
    else:
        print('cluster_reduction must be "PCA" or "T-SNE"')
        return

    groups = df[variable].values
    colors = np.array([color_mapping[group] for group in groups])

    # Utilizza plt.subplot per creare subplot
    plt.subplot(1, 3, subplot_index)  # Aggiungi questo per i subplot
    plt.scatter(reduced_features[:, 0], reduced_features[:, 1], c=colors, s=50)

    if cluster_reduction == 'PCA' and explained_variance_ratio is not None:
        plt.xlabel('PC1 ({:.2%} varianza spiegata)'.format(explained_variance_ratio[0]))
        plt.ylabel('PC2 ({:.2%} varianza spiegata)'.format(explained_variance_ratio[1]))

    else:
        plt.xlabel('PC1')
        plt.ylabel('PC2')
              
    plt.title(f'{cluster_reduction} components - Colored by {variable}')

    legend_labels = list(color_mapping.keys())
    legend_handles = [plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=color_mapping[label], markersize=10) for label in legend_labels]
    plt.legend(legend_handles, legend_labels, title=variable, bbox_to_anchor=(1.05, 1), loc='upper left')
    
    if save_fig:
        title = input('Choose your title: ')
        # Regola il layout per evitare sovrapposizioni
        plt.tight_layout()
        plt.savefig(f'{title}.png')
        print(f'Figure saved as {title}.png')

    return 


def tsne_3D(df, scaled_features, variable, save_plot=False):
    tsne3 = TSNE(n_components=3)
    components_tsne3 = tsne3.fit_transform(scaled_features)
    
    # Create a DF wth reduced datas to 3 dimensions and use 'variable' as colours
    df_3d_tsne = pd.DataFrame(data=components_tsne3, columns=['PC1', 'PC2', 'PC3'])
    df_3d_tsne['color'] = df[variable].values
    
    # Create interactive 3D graph with Plotly Express library!
    fig_3d = px.scatter_3d(df_3d_tsne, x='PC1', y='PC2', z='PC3', color='color',
                       title='T-sne plot in 3D - Colored by Age Group',
                       labels={'PC1': 'PC 1', 'PC2': 'PC 2', 'PC3': 'PC 3'},
                       width=800, height=600)
    
    return fig_3d.show()


def create_df_volcano(df,group1,group2,target_col):
    print(df[target_col].value_counts())
    
    # Transform the df in a long format
    df_long = pd.melt(df, id_vars=[target_col],
                    var_name='gene', value_name='expr')
    
    results = {}

    # Use the groupby to create am iterable Pandas object 
    grouped = df_long.groupby('gene')

    # Perform Mann-Whitney U for every gene
    for gene, group_data in grouped:
        group1_data = group_data[group_data[target_col] == group1]['expr']
        group2_data = group_data[group_data[target_col] == group2]['expr']

        statistic, p_value = mannwhitneyu(group1_data, group2_data)

        results[gene] = p_value

    # Create a df of the results
    result_df = pd.DataFrame(list(results.items()), columns=['gene', 'p-value'])

    # Use the Benjamin-Hochberg correction for p-value
    result_df['adjusted_p-value'] = multipletests(result_df['p-value'], method='fdr_bh')[1]

    # list of gene names
    gene_columns = [gene for gene in result_df.gene]
    
    # Create a empty list in which we will put the mean value of every group for every gene
    df_mean = []

    for gene in gene_columns:
        # Flter datas for the current gene
        df_gene = df[[target_col, gene]]

        # Calculate the average values
        average_group1 = df_gene[df_gene[target_col] == group1][gene].mean()
        average_group2 = df_gene[df_gene[target_col] == group2][gene].mean()

        # Add the values to the list
        df_mean.append([gene, average_group1, average_group2])

    # Create a df of the 'mean list'
    df_mean = pd.DataFrame(df_mean, columns=['Gene', group1, group2])
    
    # Add the column of the adjusted p value calculated in result_df in the new df_mean
    df_mean['adjusted_p-value'] = result_df['adjusted_p-value']
    
    # Change the name (is not necessary)
    df_volcano=df_mean
    
    # Create new column of the log2 Fold Change
    # print to explain
    print(f'Note: our value are already in log2 format, so I will use the difference among the {group2} - {group1}' )
    df_volcano['log2_FC'] = df_volcano[group2] - df_volcano[group1]
   
    # finally create the last columns of the negative log10 of adjusted p-value
    df_volcano['-log10_pval'] = np.log10(df_volcano['adjusted_p-value']) * (-1)
    

    return df_volcano



def go_library_volcano_plot(df_volcano, title1, title2, save_fig = False):
    print('This is the first and simplest graph in which you can appreciate only the scatter')
    print('I used the *plotly.graph_objects* library')
    fig = go.Figure()

    fig1 = go.Scatter(
     x=df_volcano['log2_FC'],
     y=df_volcano['-log10_pval'],
     mode='markers',
     name='',
     hovertext=list(df_volcano.index)
    )

    fig.add_trace(fig1)

    fig.update_layout(title=title1)
    fig.show()

    
    print('Second plot whith gene name')
    fig = px.scatter(df_volcano, x='log2_FC', y='-log10_pval', text=df_volcano.Gene)
    fig.update_traces(textposition='top center')
    fig.update_layout(
     title_text=title2
    )
    
    if save_fig:
        fig.write_image(f"{title2}.png", engine="kaleido") 
        
    fig.show()
    
    return
def significant_volcano_plot(df_volcano, FC, title, alpha, save_fig = False):
    plt.scatter(x=df_volcano['log2_FC'],y=df_volcano['-log10_pval'],s=1,label="Not significant")

    # highlight down- or up- regulated genes
    down = df_volcano[(df_volcano['log2_FC']<=-FC)&(df_volcano['adjusted_p-value']<=alpha)]
    up = df_volcano[(df_volcano['log2_FC']>=FC)&(df_volcano['adjusted_p-value']<=alpha)]

    plt.scatter(x=down['log2_FC'],y=down['adjusted_p-value'].apply(lambda x:-np.log10(x)),s=3,label="Down-regulated",color="blue")
    plt.scatter(x=up['log2_FC'],y=up['adjusted_p-value'].apply(lambda x:-np.log10(x)),s=3,label="Up-regulated",color="red")
    for i,r in up.iterrows():
        plt.text(x=r['log2_FC'],y=-np.log10(r['adjusted_p-value']),s=df_volcano['Gene'][i])

    for i,r in down.iterrows():
        plt.text(x=r['log2_FC'],y=-np.log10(r['adjusted_p-value']),s=df_volcano['Gene'][i])

    plt.xlabel("log2_FC")
    plt.ylabel("-log10_pval")
    plt.title(title)
    plt.axvline(-FC,color="grey",linestyle="--")
    plt.axvline(FC,color="grey",linestyle="--")
    plt.axhline(1.3,color="red",linestyle="--")
    plt.legend()
    if save_fig:
        plt.savefig(f"{title}.png", dpi = 300)
    plt.show()
    return

def dash_volcano_plot(df_volcano, title, interactive=False, save_fig = False):
    print('I used the *dash_bio* library to create this graph')
    dash_bio_fig = dash_bio.VolcanoPlot(
        dataframe=df_volcano,
        effect_size='log2_FC',
        p='adjusted_p-value',
        snp=None,
        gene='Gene'
    )
    if save_fig:
        fig.write_image(f"{title}.png", engine="kaleido") 
    
    return dash_bio_fig