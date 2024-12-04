import numpy as np
import pandas as pd
import abagen
from abagen import images
import nibabel as nib
import os, glob
from nilearn import image
from scipy import stats
import statsmodels.stats.multitest as multi
from nilearn.maskers import NiftiMasker

def validate_path(filepath, file_type="file"):
    """
    Validates if a given file or directory path exists.

    Parameters:
    filepath (str): Path to validate.
    file_type (str): Type of path - "file" or "directory".

    Raises:
    FileNotFoundError: If the path does not exist.
    """
    if file_type == "file" and not os.path.isfile(filepath):
        raise FileNotFoundError(f"File not found: {filepath}")
    elif file_type == "directory" and not os.path.isdir(filepath):
        raise FileNotFoundError(f"Directory not found: {filepath}")
    
def get_atlas_expr(atlas_path, atlas_info):

    """
    Loads gene expression data using abagen.

    Parameters:
    atlas_path (str): File path to the atlas NIfTI file.
    atlas_info (str): File path to the atlas CSV file.

    Returns:
    tuple:
        - pd.DataFrame: Gene expression data.
        - dict: Abagen report.
    """

    validate_path(atlas_path, "file")
    validate_path(atlas_info, "file")

    atlas_info = pd.read_csv(atlas_info)
    atlas = images.check_atlas(atlas_path,atlas_info=atlas_info) 

    try: 
        expressions = abagen.get_expression_data(
            atlas, 
            atlas_info=atlas_info,
            ibf_threshold=0.5, 
            probe_selection='diff_stability',
            donor_probes='aggregate',
            sim_threshold=None, 
            lr_mirror='bidirectional', 
            missing=None,
            tolerance=2,
            sample_norm='zscore',
            gene_norm='zscore', 
            norm_matched=True,
            norm_structures=False,
            region_agg='donors',
            agg_metric='mean',
            corrected_mni=False, 
            reannotated=True, 
            return_counts=False,
            return_donors=False,
            return_report=True,
            donors='all',
            data_dir='/path/to/microarray/directory',
            verbose=1,
            n_proc=1
        )

    except Exception as e:
        raise RuntimeError(f"Error processing atlas expressions: {e}")
    
    expression = expressions[0]
    report = expressions[1]
    expression.index = expression.index.astype(str)
    expression = expression.dropna(axis=0)
    return expression, report

def get_available_genes(expression, csv_path): 
    """
    Filters genes from the expression data based on a given list.

    Parameters:
    expression (pd.DataFrame): Gene expression data.
    csv_path (str): Path to CSV containing genes of interest.

    Returns:
    pd.DataFrame: Filtered gene expression data.
    """
    validate_path(csv_path, "file")
    avail_genes = expression.columns
    gene_list = pd.read_csv(csv_path)

    if 'name' not in gene_list.columns:
        raise ValueError("CSV file must contain a 'name' column with gene names.")

    goi = gene_list['name'].to_list()
    goi_ls = []
    for g in goi:
        if g in avail_genes:
            goi_ls.append(g)

    g_df = expression[goi_ls]

    return g_df

def load_parcel_fc(mask,path):
    """
    Loads functional connectivity maps masked with a brain atlas.

    Parameters:
    mask (str): Filepath to the brain mask NIfTI image.
    path (list): List of file paths pointing to functional connectivity (Fz) maps.

    Returns:
    tuple: 
        - NiftiMasker: A Nilearn NiftiMasker object fitted with the brain mask.
        - np.ndarray: Masked functional connectivity maps as an array.
    """
    validate_path(mask, "file")
    for p in path:
        validate_path(p, "file")

    nifti_masker = NiftiMasker(mask_img=mask)

    try:
        fmri_masked = nifti_masker.fit_transform(path)
    except Exception as e:
        raise RuntimeError(f"Error loading parcel FC maps: {e}")

    return nifti_masker, fmri_masked

def corr2_coeff(A, B):
    """
    Computes the correlation coefficients and p-values for the relationship 
    between two datasets (A and B) using a beta distribution approximation.

    Parameters:
    A (np.ndarray): A 2D array where each row represents a variable, and each column is an observation.
    B (np.ndarray): A 2D array with the same structure as A.

    Returns:
    tuple:
        - corrs (np.ndarray): Flattened array of correlation coefficients.
        - ps (np.ndarray): Flattened array of p-values corresponding to the correlations.
    
    Raises:
    ValueError: If A or B is not a 2D array or if their shapes are incompatible.
    """
    if not (isinstance(A, np.ndarray) and isinstance(B, np.ndarray)):
        raise TypeError("Inputs A and B must be numpy arrays.")
    
    if A.ndim != 2 or B.ndim != 2:
        raise ValueError("Inputs A and B must be 2D arrays.")
    
    if A.shape[1] != B.shape[1]:
        raise ValueError("Inputs A and B must have the same number of columns (observations).")

    try:
        n = A.shape[1]

        if n <= 2:
            raise ValueError("Each array must have at least 3 observations (columns) for correlation calculation.")
        
        dist = stats.beta(n/2 - 1, n/2 - 1, loc=-1, scale=2)
        
        A_mA = A - A.mean(1)[:, None]
        B_mB = B - B.mean(1)[:, None]

        ssA = (A_mA**2).sum(1)
        ssB = (B_mB**2).sum(1)

        corrs = (np.dot(A_mA, B_mB.T) / np.sqrt(np.dot(ssA[:, None],ssB[None]))).ravel()
        ps = 2*dist.cdf(-abs(corrs))
        ps = np.around(ps, 200)

        return corrs, ps
    
    except Exception as e:
        raise RuntimeError(f"Error in correlation coefficient calculation: {e}")

def save_nifti_images(R_map, p_map, r_method, nifti_masker, expr_col, corr_method, output_dir):
    """
    Saves NIfTI images of the results.

    Parameters:
    expr_col (str): Column name of the expression data.
    R_thr (ndarray): Thresholded correlation coefficients.
    p_corrected (ndarray): Corrected p-values.
    nifti_masker (NiftiMasker): A Nilearn NiftiMasker object fitted with the brain mask.
    output_dir (str): Directory where the Nifti images will be saved.
    corr_method (str): Correlation method used.
    """
    try:
        new_R_image = nifti_masker.inverse_transform(R_map)
        new_p_image = nifti_masker.inverse_transform(p_map)

        if r_method == 'R2':
            nib.save(new_R_image, os.path.join(output_dir, f'R2_{expr_col}_{corr_method}.nii.gz'))
            nib.save(new_p_image, os.path.join(output_dir, f'p_{expr_col}_{corr_method}.nii.gz'))
        else:
            nib.save(new_R_image, os.path.join(output_dir, f'R_{expr_col}_{corr_method}.nii.gz'))
            nib.save(new_p_image, os.path.join(output_dir, f'p_{expr_col}_{corr_method}.nii.gz'))
    except Exception as e:
        raise RuntimeError(f"Error saving NIfTI images: {e}")

def generate_model(expr_data, fz_mx, corr_method):
    """
    Computes association between gene expression data and functional connectivity matrix.

    Parameters:
    expr_data (array-like): Gene expression data.
    fz_mx (array-like): Functional connectivity matrix.
    corr_method (str): Correlation method ('pearson' or 'spearman').

    Returns:
    tuple: Tuple containing:
        - R (ndarray): Correlation coefficients.
        - p (ndarray): P-values for the correlation.
    """
    if corr_method not in ['pearson', 'spearman']:
        raise ValueError("Invalid correlation method. Choose 'pearson' or 'spearman'.")
    
    try:
        if corr_method == 'pearson':
            R, p = corr2_coeff(fz_mx.T, np.reshape(np.array(expr_data), (-1, 1)).T)
        elif corr_method == 'spearman':
            fz_rank = stats.rankdata(fz_mx, axis=1)
            expr_rank = np.array(expr_data).reshape(1, -1)
            expr_rank = stats.rankdata(expr_rank, axis=1)
            R, p = corr2_coeff(fz_rank.T, expr_rank)
    except Exception as e:
        raise RuntimeError(f"Error generating correlation model: {e}")

    return R, p

def apply_multiple_comparisons_correction(R, p, method, threshold):
    """
    Applies multiple comparisons correction.

    Parameters:
    R (ndarray): Correlation coefficients.
    p (ndarray): P-values for the correlation.
    method (str): Correction method ('FDR', 'Bonferroni', or None).
    threshold (float): Significance threshold.

    Returns:
    tuple: Tuple containing:
        - R_thr (ndarray): Thresholded correlation coefficients.
        - p_corrected (ndarray): Corrected p-values.
    """

    if method not in ['FDR', 'Bonferroni', None]:
        raise ValueError("Invalid correction method. Choose 'FDR', 'Bonferroni', or None.")

    try:
        R[np.isnan(R)] = 0
        p[np.isnan(p)] = 0
        mask = R.copy()
        mask[mask != 0] = 1

        if method == 'FDR':
            p[mask==1] = multi.fdrcorrection(p[mask == 1], alpha=threshold)[0]
        elif method == 'Bonferroni':
            p[mask==1] = multi.multipletests(p[mask == 1], alpha=threshold, method='bonferroni')[0]
        else:
            p[p >= threshold] = 0
            p[p != 0] = 1

        R_thr = R.copy()
        R_thr[p != 1.] = 0

        return R_thr, p
    
    except Exception as e:
        raise RuntimeError(f"Error during multiple comparisons correction: {e}")

def simple_r_map(data, fz_mx, r_method, multcomp, nifti_masker, save_first_level, output_dir=None, p_thr=0.05, corr='pearson'):
    """
    Generates first-level R and p maps from expression data and functional connectivity.

    Parameters:
    data (pd.DataFrame): Gene expression data for regions of interest.
    fz_mx (np.ndarray): Functional connectivity maps.
    r_method (str): Type of correlation ('R2' or 'R').
    multcomp (str): Multiple comparisons correction method ('FDR' or 'Bonferroni').
    nifti_masker (NiftiMasker): A Nilearn NiftiMasker object fitted with the brain mask.
    save_first_level (bool): Whether to save first-level maps as NIfTI images.
    output_dir (str, optional): Directory to save outputs. Defaults to None.
    p_thr (float, optional): P-value threshold for significance. Defaults to 0.05.
    corr (str, optional): Correlation method ('pearson' or 'spearman').

    Returns:
    tuple:
        - list: List of thresholded R maps as NIfTI images.
        - list: List of corrected p maps as NIfTI images.
    """
    if r_method not in ['R2', 'R']:
        raise ValueError("Invalid r_method. Choose 'R2' or 'R'.")
    if multcomp not in ['FDR', 'Bonferroni']:
        raise ValueError("Invalid multcomp method. Choose 'FDR' or 'Bonferroni'.")
    if not isinstance(data, pd.DataFrame):
        raise TypeError("data must be a pandas DataFrame.")
    if not isinstance(fz_mx, np.ndarray):
        raise TypeError("fz_mx must be a numpy ndarray.")
    
    final_R_maps = []
    final_p_maps = []

    for expr_col in data:
        try:
            R, p = generate_model(data[expr_col], fz_mx, corr)
        except Exception as e:
            raise RuntimeError(f"Error generating model for {expr_col}: {e}")
        
        if output_dir:
            try:
                R_thr, p_corrected = apply_multiple_comparisons_correction(R, p, multcomp, p_thr)
            except Exception as e:
                raise RuntimeError(f"Error applying multiple comparisons correction: {e}")

            if r_method == 'R2':
                R_thr = R_thr**2

            try:
                final_R_maps.append(nifti_masker.inverse_transform(R_thr))
                final_p_maps.append(nifti_masker.inverse_transform(p_corrected))

                if save_first_level==True: 
                    save_nifti_images(R_thr, p_corrected, r_method, nifti_masker, expr_col, corr, output_dir)
            except Exception as e:
                raise RuntimeError(f"Error processing NIfTI images for {expr_col}: {e}")

    return final_R_maps, final_p_maps

def first_level(expression, fz, r_method, multcomp, masker, output_dir, p_thr, save_first_level):
    """
    Generates gene-network maps based on gene expression and functional connectivity data.

    Parameters:
    expression (pd.DataFrame): Gene expression data for regions of interest.
    fz (np.ndarray): Functional connectivity matrix (Fz maps).
    r_method (str): Correlation metric ('R2' or 'R').
    multcomp (str): Method for multiple comparisons correction ('FDR' or 'Bonferroni').
    masker (NiftiMasker): A Nilearn NiftiMasker object fitted with the brain mask.
    output_dir (str): Directory to save output files.
    p_thr (float): P-value threshold for significance.
    save_first_level (bool): Whether to save first-level maps as NIfTI files.

    Returns:
    tuple:
        - list: List of thresholded R maps as NIfTI images.
        - list: List of corrected p maps as NIfTI images.
    """
    if r_method not in ['R2', 'R']:
        raise ValueError("Invalid r_method. Choose 'R2' or 'R'.")
    if multcomp not in ['FDR', 'Bonferroni']:
        raise ValueError("Invalid multcomp method. Choose 'FDR' or 'Bonferroni'.")

    try:
        r_map, p_map = simple_r_map(
            data=expression, 
            fz_mx=fz, 
            r_method=r_method, 
            multcomp=multcomp, 
            nifti_masker=masker, 
            output_dir=output_dir, 
            p_thr=p_thr, 
            corr='pearson', 
            save_first_level=save_first_level
        )

    except Exception as e:
        raise RuntimeError(f"Error in first-level analysis: {e}")

    return r_map, p_map

def generate_n_map(r_map, n_thr, expression):
    """
    Generates an 'n-map', indicating the number of significant correlations across gene network maps.

    Parameters:
    r_map (list of Nifti1Image): List of first-level R maps.
    n_thr (float): Fractional threshold for the number of significant correlations.
    expression (pd.DataFrame): Gene expression data used to determine the total number of genes.

    Returns:
    Nifti1Image: Thresholded n-map indicating significant correlation counts.
    """
    if not isinstance(r_map, list) or not all(isinstance(i, nib.Nifti1Image) for i in r_map):
        raise TypeError("r_map must be a list of Nifti1Image objects.")
    if not isinstance(expression, pd.DataFrame):
        raise TypeError("expression must be a pandas DataFrame.")

    try:
        n_map = np.array([i.get_fdata() for i in r_map])
        n_map[n_map != 0] = 1
        n_map = n_map.sum(axis=0)
        n_thr2 = round(expression.shape[1] * n_thr)
        n_map = nib.Nifti1Image(n_map, affine=r_map[0].affine)
        thr_n_map = image.threshold_img(n_map, threshold=n_thr2)
    except Exception as e:
        raise RuntimeError(f"Error generating n-map: {e}")

    return thr_n_map

def generate_LNM(r_map, n_thr, expression):
    """
    Uses the Lesion Network Mapping (LNM) approach, which combines positive and negative correlations across all R maps.

    Parameters:
    r_map (list of Nifti1Image): List of first-level R maps.
    n_thr (float): Fractional threshold for the number of significant correlations.
    expression (pd.DataFrame): Gene expression data used to determine the total number of genes.

    Returns:
    Nifti1Image: Thresholded Local Network Map (LNM).
    """
    if not isinstance(r_map, list) or not all(isinstance(i, nib.Nifti1Image) for i in r_map):
        raise TypeError("r_map must be a list of Nifti1Image objects.")
    if not isinstance(expression, pd.DataFrame):
        raise TypeError("expression must be a pandas DataFrame.")

    try:
        lnm = np.array([i.get_fdata() for i in r_map])
        lnm[lnm>0] = 1
        lnm[lnm<0] = -1
        lnm = lnm.sum(axis=0)
        n_thr2 = round(expression.shape[1] * n_thr)
        lnm = nib.Nifti1Image(lnm, affine=r_map[0].affine)
        thr_lnm = image.threshold_img(lnm, threshold=n_thr2)
    except Exception as e:
        raise RuntimeError(f"Error generating Local Network Map (LNM): {e}")

    return thr_lnm

def second_level(data, expression, sl_method, r_method, n_thr, output_dir, system_of_interest, multcomp, p_thr, masker):
    """
    Generates second-level correlation maps using a specified method (e.g., n-map or LNM).

    Parameters:
    data (list of Nifti1Image): List of first-level R maps.
    expression (pd.DataFrame): Gene expression data used to determine thresholds.
    sl_method (str): Second-level mapping method ('n-map' or 'LNM').
    r_method (str): Correlation metric ('R2'or 'R').
    n_thr (float): Fractional threshold for the number of significant correlations.
    output_dir (str): Directory to save output files.
    system_of_interest (str): Identifier for the system of interest.
    multcomp (str): Identifier for the multiple comparisons correction method.
    p_thr (float): Identifier for the P-value threshold for significance.
    masker (NiftiMasker): A Nilearn NiftiMasker object fitted with the brain mask.

    Returns:
    Nifti1Image: Second-level correlation map saved as a NIfTI file.
    """
    if sl_method not in ['n-map', 'LNM']:
        raise ValueError("Invalid second-level method. Choose 'n-map' or 'LNM'.")

    try:
        if sl_method == 'n-map':
            second_level_map = generate_n_map(data, n_thr, expression)
        elif sl_method == 'LNM':
            second_level_map = generate_LNM(data, n_thr, expression)

        second_level_map.to_filename(f'{output_dir}{system_of_interest}_{sl_method}_p{p_thr}_multcomp{multcomp}_n_{n_thr}_{r_method}.nii.gz')
        
        return second_level_map
    
    except Exception as e:
        raise RuntimeError(f"Error generating second-level map: {e}")

#Define variables
output_dir = '/path/to/output/directory/'
#Define atlas parcellations to be used with abagen
atlas = '/path/to/Translating_The_Transcriptome/Gene_Network_Mapping/Atlas_parcellation.nii.gz'
atlas_csv = '/path/to/Translating_The_Transcriptome/Gene_Network_Mapping/Atlas_parcellation_annotation_table.csv'
#Define brain mask
masks = '/path/to/Translating_The_Transcriptome/Gene_Network_Mapping/Brain_mask.nii.gz'
#Define list containing systems of interest
system_of_interest = '/path/to/Translating_The_Transcriptome/Gene_Network_Mapping/Gene_lists/Parkinsonism_genes.csv'
#Set path to fc-maps
atlas_file_suffix = '*.nii.gz'
path_fc = 'path/to/fingerprint/maps/'
#Get expressions
gene_expression, report = get_atlas_expr(atlas,atlas_csv)
#Load fc-maps
fc_path = [path_fc+str(i)+atlas_file_suffix for i in gene_expression.index]
masker, fz = load_parcel_fc(masks,fc_path)
#generate gene- and disease network maps
expression = get_available_genes(gene_expression, system_of_interest)
r_map, _ = first_level(expression=expression, fz=fz, r_method='R2', multcomp='FDR', masker=masker, output_dir=output_dir, p_thr=0.05, save_first_level=True)
processed_image = second_level(r_map, expression, sl_method='LNM', r_method='R2', n_thr=0.0, output_dir=output_dir, system_of_interest=system_of_interest, multcomp='FDR', p_thr=0.05, masker=masker)