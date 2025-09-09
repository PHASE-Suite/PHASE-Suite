#!/usr/bin/env python3
# -*- coding: utf-8 -*-

################################################################################
# Script Name: Interactive ecNGS Data Visualization App
# Author: Jaime Miranda
# Date: September 2, 2025
# Version: 5.0
#
# Description:
# This Streamlit web application serves as the interactive visualization component
# for the PHASE (PacBio HiFi Analysis of Somatic Events) Suite. It transforms
# the output from the PHASE bioinformatics pipeline (.bases and .muts files)
# into a user-friendly graphical interface, enabling researchers to generate
# publication-quality plots and data tables without command-line interaction.
#
# The application is designed for maximum flexibility, handling various file
# naming conventions and allowing for user-defined experimental controls.
#
# Key Features:
# 1.  Flexible File Handling: Automatically parses filenames with an
#     'Identifier_Condition.bcBarcode' format.
# 2.  Smart Control Detection: Can perform standard statistical analysis against
#     a simple control condition OR a time-matched, shared control identifier.
# 3.  Interactive Visualizations: Utilizes the Plotly library to create dynamic
#     charts where users can hover for details, zoom, and pan.
# 4.  Rigorous Statistical Analysis: Implements a One-Way ANOVA followed by a
#     Dunnett's post-hoc test for simple designs, and time-matched t-tests for
#     complex, two-factor designs.
# 5.  Advanced Plot Customization: Offers granular control over plot titles,
#     axis labels, and Y-axis synchronization for publication-ready figures.
# 6.  Batch Export: Allows for one-click download of all generated plots and
#     data tables in a single, organized .zip archive.
#
# Input Requirements:
# -   An 'indexing' folder containing two subdirectories: 'bases/' and 'mutations/'.
# -   '.bases' files: Text files containing the total number of bases analyzed.
# -   '.muts' files: Tab-separated files detailing identified mutations.
# -   Reference Genome: A FASTA file (.fa) for trinucleotide context extraction.
#
# Python Dependencies:
# - streamlit, pandas, numpy, matplotlib, scipy, pyfaidx, scikit-posthocs,
#   plotly, kaleido
#
################################################################################


import os
import glob
import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import sem, f_oneway, ttest_ind
from pyfaidx import Fasta
import scikit_posthocs as sp
import matplotlib.colors as mcolors
import streamlit as st
from io import BytesIO
import zipfile
import plotly.graph_objects as go
from plotly.subplots import make_subplots

# --- Configuration and Constants ---

# Maps base changes to the 6 standard substitution groups.
MUTATION_GROUPS = {
    ("G", "T"): "G:C>T:A", ("C", "A"): "G:C>T:A",
    ("G", "C"): "G:C>C:G", ("C", "G"): "G:C>C:G",
    ("G", "A"): "G:C>A:T", ("C", "T"): "G:C>A:T",
    ("A", "T"): "A:T>T:A", ("T", "A"): "A:T>T:A",
    ("A", "G"): "A:T>G:C", ("T", "C"): "A:T>G:C",
    ("A", "C"): "A:T>C:G", ("T", "G"): "A:T>C:G"
}
# Defines the standard order for plotting the 6 mutation groups.
DESIRED_GROUP_ORDER = [
    "G:C>T:A", "G:C>C:G", "G:C>A:T",
    "A:T>T:A", "A:T>G:C", "A:T>C:G",
]
# Defines the color for each of the 6 mutation groups in summary plots.
GROUP_COLOR_MAP = {
    "G:C>T:A": "blue", "G:C>C:G": "black", "G:C>A:T": "red",
    "A:T>T:A": "grey", "A:T>G:C": "green", "A:T>C:G": "pink",
}

# Defines the 6 canonical substitution classes for 96-channel signatures.
SUBSTITUTION_CLASSES = ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"]
# Defines the standard COSMIC colors for each substitution class.
SUBSTITUTION_COLORS = {
    "C>A": "#03BCE3", "C>G": "#000000", "C>T": "#E32221",
    "T>A": "#A0A0A0", "T>C": "#008101", "T>G": "#E4A9C8"
}
# Defines the complementary base for DNA strands.
COMPLEMENT_MAP = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}


# --- Helper and Utility Functions ---

def to_csv(df):
    """Converts a pandas DataFrame to a UTF-8 encoded CSV string for downloading."""
    return df.to_csv(index=False).encode('utf-8')

def find_data_files(directory, extension):
    """Discovers all files in a directory that have a given extension."""
    if not os.path.exists(directory) or not os.path.isdir(directory):
        return []
    return glob.glob(os.path.join(directory, f"*.{extension}"))

def generate_96_channel_order():
    """Generates the standard 96-channel trinucleotide context order for plotting."""
    order = []
    for sub in SUBSTITUTION_CLASSES:
        ref, alt = sub[0], sub[-1]
        for left in "ACGT":
            for right in "ACGT":
                context_label = f"{left}[{ref}>{alt}]{right}"
                order.append((sub, context_label))
    return order

STANDARD_96_ORDER = generate_96_channel_order()

def flexible_sort_key(condition_str, control_condition=None):
    if condition_str == control_condition:
        return (0, 0)
    numbers = re.findall(r'(\d+\.?\d*)', condition_str)
    if numbers:
        return (1, float(numbers[0]))
    return (2, condition_str)

def canonical_96_context(context_str):
    try:
        parts = context_str.split('[')
        left_part = parts[0]
        parts2 = parts[1].split(']')
        mid_part = parts2[0]
        right_part = parts2[1]
        left_base = left_part[-1] if left_part else "N"
        ref_base, alt_base = mid_part.split('>')
        right_base = right_part[0] if right_part else "N"
    except Exception:
        return None, None
    if ref_base in ["A", "G"]:
        left_base, ref_base, alt_base, right_base = (
            COMPLEMENT_MAP[right_base], COMPLEMENT_MAP[ref_base],
            COMPLEMENT_MAP[alt_base], COMPLEMENT_MAP[left_base]
        )
    sub_class = f"{ref_base}>{alt_base}"
    if sub_class not in SUBSTITUTION_CLASSES:
        return None, None
    final_label = f"{left_base}[{ref_base}>{alt_base}]{right_base}"
    return sub_class, final_label

def get_baseline_context_counts(reference):
    first_contig = next(iter(reference.keys()))
    baseline_seq = str(reference[first_contig][:300000]).upper()
    baseline_counts = {}
    total_baseline = 0
    for i in range(1, len(baseline_seq) - 1):
        trinuc = baseline_seq[i-1:i+2]
        left, center, right = trinuc
        if center in ["C", "T"]:
            context_label = f"{left}[{center}]{right}"
        else:
            new_left, new_center, new_right = COMPLEMENT_MAP[right], COMPLEMENT_MAP[center], COMPLEMENT_MAP[left]
            context_label = f"{new_left}[{new_center}]{new_right}"
        if context_label:
            baseline_counts[context_label] = baseline_counts.get(context_label, 0) + 1
            total_baseline += 1
    return baseline_counts, total_baseline

def save_plotly_fig_to_buffer(fig, format='png'):
    """Saves a Plotly figure to an in-memory BytesIO buffer for downloading."""
    buf = BytesIO()
    fig.write_image(buf, format=format)
    buf.seek(0)
    return buf

# --- Data Processing Core Logic ---

def process_files(selected_bases_paths, mutations_dir, reference):
    samples_data = []
    pattern = re.compile(r'^(?P<identifier>[^_]+)_(?P<condition>.+?)\.bc(?P<barcode>\d+)')
    for base_file_path in selected_bases_paths:
        base_filename = os.path.basename(base_file_path)
        mut_filename = base_filename.replace('.bases', '.muts')
        corresponding_mut_file = os.path.join(mutations_dir, mut_filename)
        if not os.path.exists(corresponding_mut_file):
            st.warning(f"Skipping `{base_filename}`: Corresponding mutation file `{mut_filename}` not found.")
            continue
        match = pattern.match(base_filename)
        if not match:
            st.warning(f"Skipping `{base_filename}`: Filename does not match expected format.")
            continue
        identifier, condition, barcode = match.groups()
        try:
            with open(base_file_path, 'r') as f: bases_analyzed = float(f.read().strip())
            with open(corresponding_mut_file, 'r') as f: mutation_lines = [line for line in f.readlines() if line.strip()]
        except Exception as e:
            st.error(f"Error reading files for `{base_filename}`: {e}")
            continue
        num_mutations = len(mutation_lines)
        mutation_group_counts, mutation_context_96_counts = {}, {}
        for line in mutation_lines:
            cols = line.strip().split('\t')
            if len(cols) < 4: continue
            chrom, pos_str, ref_base, alt_base = cols[0], cols[1], cols[2].upper(), cols[3].upper()
            group = MUTATION_GROUPS.get((ref_base, alt_base))
            if group: mutation_group_counts[group] = mutation_group_counts.get(group, 0) + 1
            try:
                pos = int(pos_str)
                left_base = str(reference[chrom][pos-2]).upper() if pos > 1 else "N"
                right_base = str(reference[chrom][pos]).upper() if pos < len(reference[chrom]) else "N"
                raw_context = f"{left_base}[{ref_base}>{alt_base}]{right_base}"
            except Exception:
                raw_context = "N[N>N]N"
            sub_class, final_label = canonical_96_context(raw_context)
            if sub_class and final_label:
                key = (sub_class, final_label)
                mutation_context_96_counts[key] = mutation_context_96_counts.get(key, 0) + 1
        mutation_freq = num_mutations / bases_analyzed if bases_analyzed > 0 else 0
        samples_data.append({
            "identifier": identifier, "condition": condition, "barcode": barcode, "bases": bases_analyzed,
            "num_mutations": num_mutations, "mutation_freq": mutation_freq,
            "mutation_group_counts": mutation_group_counts,
            "mutation_context_96_counts": mutation_context_96_counts
        })
    return pd.DataFrame(samples_data) if samples_data else pd.DataFrame()

# --- Interactive Plotting Functions (Plotly) ---

def plot_dose_response_plotly(df, overlay=False, control_selection=None, y_max=None, custom_title=None, xaxis_title=None, yaxis_title=None):
    if df.empty or 'mutation_freq' not in df.columns: return None, None
    agg_df = df.groupby(['identifier', 'condition'])['mutation_freq'].agg(['mean', sem, 'count']).reset_index()
    
    # --- Smart Statistics Block ---
    if control_selection:
        # Check if we are doing a time-matched comparison
        is_time_matched = (control_selection == "Identifier: Control (Time-Matched)")
        stats_results = []
        
        # --- LOGIC FOR TIME-MATCHED T-TESTS ---
        if is_time_matched:
            control_df = df[df['identifier'] == 'Control']
            treatment_identifiers = [i for i in df['identifier'].unique() if i != 'Control']
            
            for identifier in treatment_identifiers:
                identifier_df = df[df['identifier'] == identifier]
                for condition in identifier_df['condition'].unique():
                    treatment_freqs = identifier_df[identifier_df['condition'] == condition]['mutation_freq']
                    control_freqs = control_df[control_df['condition'] == condition]['mutation_freq']
                    
                    if len(treatment_freqs) >= 2 and len(control_freqs) >= 2:
                        _, p_val = ttest_ind(treatment_freqs, control_freqs, equal_var=False, nan_policy='omit')
                        stats_results.append({'identifier': identifier, 'condition': condition, 'p_value_vs_control': p_val})
                    else:
                        stats_results.append({'identifier': identifier, 'condition': condition, 'p_value_vs_control': np.nan})

        # --- LOGIC FOR ANOVA/DUNNETT'S TEST ---
        else:
            control_condition = control_selection
            for identifier in df['identifier'].unique():
                identifier_df = df[df['identifier'] == identifier]
                groups = [group['mutation_freq'].values for name, group in identifier_df.groupby('condition')]
                
                if len(groups) > 1 and all(len(g) >= 2 for g in groups):
                    _, anova_p = f_oneway(*groups)
                    dunnett_p_values = sp.posthoc_dunnett(identifier_df, val_col='mutation_freq', group_col='condition', control=control_condition)
                    for cond in identifier_df['condition'].unique():
                        p_val = dunnett_p_values.loc[cond, control_condition] if cond in dunnett_p_values.index else np.nan
                        stats_results.append({'identifier': identifier, 'condition': cond, 'ANOVA_p_value': anova_p, 'p_value_vs_control': p_val})
                else:
                    for cond in identifier_df['condition'].unique():
                        stats_results.append({'identifier': identifier, 'condition': cond, 'ANOVA_p_value': np.nan, 'p_value_vs_control': np.nan})

        # Merge results back into the aggregated dataframe
        if stats_results:
            stats_df = pd.DataFrame(stats_results)
            agg_df = pd.merge(agg_df, stats_df, on=['identifier', 'condition'], how='left')
            def get_significance_stars(p):
                if pd.isna(p): return 'n/a'
                if p < 0.001: return '***'
                if p < 0.01: return '**'
                if p < 0.05: return '*'
                return 'ns'
            agg_df['significance'] = agg_df['p_value_vs_control'].apply(get_significance_stars)

    # Use the simple string for sorting if not time-matched
    sort_control = control_selection if not is_time_matched else None
    unique_conditions = sorted(df['condition'].unique(), key=lambda c: flexible_sort_key(c, sort_control))
    figs = {}

    # (Plotting logic uses agg_df)
    if overlay:
        fig = go.Figure()
        # In overlay mode, we plot all identifiers, including the 'Control' identifier if it exists
        plot_identifiers = df['identifier'].unique()
        for i, identifier in enumerate(plot_identifiers):
            cell_df = agg_df[agg_df['identifier'] == identifier].set_index('condition').reindex(unique_conditions).reset_index()
            fig.add_trace(go.Scatter(x=cell_df['condition'], y=cell_df['mean'], error_y=dict(type='data', array=cell_df['sem'], visible=True), mode='lines+markers', name=identifier))
        title = custom_title if custom_title else "Comparative Response by Condition"
        x_ax_title = xaxis_title if xaxis_title else "Condition"
        y_ax_title = yaxis_title if yaxis_title else "Mutation Frequency"
        fig.update_layout(title=title, xaxis_title=x_ax_title, yaxis_title=y_ax_title, yaxis_tickformat=".1e")
        if y_max: fig.update_yaxes(range=[0, y_max])
        figs["combined_plot"] = fig
    else:
        plot_identifiers = [i for i in df['identifier'].unique() if not (is_time_matched and i == 'Control')]
        for identifier in plot_identifiers:
            fig = go.Figure()
            cell_df = agg_df[agg_df['identifier'] == identifier].set_index('condition').reindex(unique_conditions).dropna(subset=['mean']).reset_index()
            n_conditions = len(cell_df['condition'])
            cmap = plt.get_cmap("Pastel1", n_conditions)
            bar_colors = []
            color_idx = 0
            for cond in cell_df['condition']:
                # The concept of a single 'control' color is less relevant in time-matched mode
                bar_colors.append(mcolors.to_hex(cmap(color_idx)))
                color_idx += 1
            fig.add_trace(go.Bar(x=cell_df['condition'], y=cell_df['mean'], error_y=dict(type='data', array=cell_df['sem'], visible=True), marker=dict(color=bar_colors, line=dict(color='black', width=1)), name=identifier))
            if control_selection and 'significance' in cell_df.columns:
                for i, row in cell_df.iterrows():
                    sig = row['significance']
                    if sig not in ['ns', 'n/a']:
                        y_pos = row['mean'] + row['sem'] if pd.notna(row['sem']) else row['mean']
                        fig.add_annotation(x=row['condition'], y=y_pos, text=sig, showarrow=False, yshift=10)
            title = custom_title if custom_title else f"Response by Condition for {identifier}"
            x_ax_title = xaxis_title if xaxis_title else "Condition"
            y_ax_title = yaxis_title if yaxis_title else "Mutation Frequency"
            fig.update_layout(title=title, xaxis_title=x_ax_title, yaxis_title=y_ax_title, yaxis_tickformat=".1e", showlegend=False)
            if y_max: fig.update_yaxes(range=[0, y_max])
            figs[identifier] = fig
    return figs, agg_df

def plot_mutation_groups_plotly(df, control_condition=None, y_max=None, custom_title=None, yaxis_title=None):
    # This function is not affected by the new stats logic and remains the same
    if df.empty: return None, None
    figs, mutation_group_list = {}, []
    for identifier in df['identifier'].unique():
        cell_df, conditions = df[df['identifier'] == identifier], sorted(df[df['identifier'] == identifier]['condition'].unique(), key=lambda c: flexible_sort_key(c, control_condition))
        if not conditions: continue
        fig = make_subplots(rows=1, cols=len(conditions), shared_yaxes=True, subplot_titles=conditions)
        for i, condition in enumerate(conditions):
            condition_df, total_bases, group_counts = cell_df[cell_df['condition'] == condition], cell_df[cell_df['condition'] == condition]['bases'].sum(), {}
            for _, row in condition_df.iterrows():
                for grp, count in row['mutation_group_counts'].items(): group_counts[grp] = group_counts.get(grp, 0) + count
            row_data, freqs = {"identifier": identifier, "condition": condition, "total_bases": total_bases}, []
            for grp in DESIRED_GROUP_ORDER:
                freq = group_counts.get(grp, 0) / total_bases if total_bases > 0 else 0
                freqs.append(freq); row_data[f"{grp}_freq"] = freq
            mutation_group_list.append(row_data)
            fig.add_trace(go.Bar(x=DESIRED_GROUP_ORDER, y=freqs, marker_color=[GROUP_COLOR_MAP[g] for g in DESIRED_GROUP_ORDER], name=condition), row=1, col=i+1)
        title = custom_title if custom_title else f"Mutation Frequency by Group for {identifier}"
        y_ax_title = yaxis_title if yaxis_title else "Mutation Frequency"
        fig.update_layout(title_text=title, showlegend=False, yaxis_tickformat=".1e")
        fig.update_yaxes(title_text=y_ax_title, row=1, col=1)
        if y_max: fig.update_yaxes(range=[0, y_max], row=1, col=1)
        figs[identifier] = fig
    return figs, pd.DataFrame(mutation_group_list)

def plot_96_channel_signature_plotly(df, baseline_counts, total_baseline, mode='traditional', control_selection=None, custom_title=None, xaxis_title=None, yaxis_title=None):
    if df.empty: return None, None
    plot_data, figs = [], {}
    is_time_matched = (control_selection == "Identifier: Control (Time-Matched)")
    
    for identifier in df['identifier'].unique():
        if is_time_matched and identifier == 'Control': continue # Don't plot the control itself in this mode
        cell_df = df[df['identifier'] == identifier]
        conditions = sorted(cell_df['condition'].unique(), key=lambda c: flexible_sort_key(c, None if is_time_matched else control_selection))
        
        for condition in conditions:
            # --- Smart control selection for subtraction ---
            control_norm_values = {}
            if mode == 'subtracted':
                control_df = pd.DataFrame()
                if is_time_matched:
                    # Find the control data that matches both the identifier 'Control' and the current condition
                    control_df = df[(df['identifier'] == 'Control') & (df['condition'] == condition)]
                elif control_selection:
                    # Original logic: find control data within the same identifier
                    control_df = cell_df[cell_df['condition'] == control_selection]

                if not control_df.empty:
                    total_bases_control = control_df['bases'].sum()
                    control_counts_96 = {}
                    for _, row in control_df.iterrows():
                        for key, count in row['mutation_context_96_counts'].items():
                            control_counts_96[key] = control_counts_96.get(key, 0) + count
                    for sub_class, context_label in STANDARD_96_ORDER:
                        observed_count = control_counts_96.get((sub_class, context_label), 0)
                        baseline_label = f"{context_label[0]}[{context_label[2]}]{context_label[-1]}"
                        baseline_freq = baseline_counts.get(baseline_label, 0) / total_baseline if total_baseline > 0 else 0
                        expected_count = baseline_freq * total_bases_control
                        control_norm_values[(sub_class, context_label)] = observed_count / expected_count if expected_count > 0 else 0
                else:
                    st.warning(f"Control data for '{identifier} - {condition}' not found. Cannot generate subtracted plot.")
                    continue
            
            # (Rest of the plotting logic remains the same)
            condition_df = cell_df[cell_df['condition'] == condition]
            total_bases_condition, counts_96 = condition_df['bases'].sum(), {}
            for _, row in condition_df.iterrows():
                for key, count in row['mutation_context_96_counts'].items(): counts_96[key] = counts_96.get(key, 0) + count
            total_mutations = sum(counts_96.values())
            y_values, x_labels, bar_colors = [], [], []
            for sub_class, context_label in STANDARD_96_ORDER:
                observed_count = counts_96.get((sub_class, context_label), 0)
                baseline_label = f"{context_label[0]}[{context_label[2]}]{context_label[-1]}"
                baseline_freq = baseline_counts.get(baseline_label, 0) / total_baseline if total_baseline > 0 else 0
                row_data = {"identifier": identifier, "condition": condition, "sub_class": sub_class, "context_label": context_label}
                value = 0
                if mode == 'traditional':
                    observed_freq = observed_count / total_mutations if total_mutations > 0 else 0
                    value = observed_freq / baseline_freq if baseline_freq > 0 else 0
                    row_data.update({"observed_freq": observed_freq, "baseline_freq": baseline_freq, "adjusted_value": value})
                elif mode == 'new':
                    expected_count = baseline_freq * total_bases_condition
                    value = observed_count / expected_count if expected_count > 0 else 0
                    row_data.update({"observed_count": observed_count, "expected_count": expected_count, "new_norm_value": value})
                elif mode == 'subtracted':
                    expected_count = baseline_freq * total_bases_condition
                    condition_norm_value = observed_count / expected_count if expected_count > 0 else 0
                    control_value = control_norm_values.get((sub_class, context_label), 0)
                    value = max(0, condition_norm_value - control_value)
                    row_data.update({"condition_norm_value": condition_norm_value, "control_norm_value": control_value, "diff_value": value})
                y_values.append(value); plot_data.append(row_data); x_labels.append(context_label); bar_colors.append(SUBSTITUTION_COLORS[sub_class])
            fig = go.Figure(data=[go.Bar(x=x_labels, y=y_values, marker_color=bar_colors)])
            title_map = {'traditional': "Traditional", 'new': "New", 'subtracted': "New (Control Subtracted)"}
            ylabel_map = {'traditional': "Enrichment", 'new': "Normalized Value", 'subtracted': "Difference"}
            default_title = f"{title_map[mode]} Normalization: {identifier} - {condition}"
            x_ax_title = xaxis_title if xaxis_title else "Trinucleotide Context"
            y_ax_title = yaxis_title if yaxis_title else ylabel_map[mode]
            fig.update_layout(title=custom_title if custom_title else default_title, xaxis_title=x_ax_title, yaxis_title=y_ax_title, xaxis_tickangle=-90)
            if mode in ['new', 'subtracted']: fig.update_layout(yaxis_tickformat=".1e")
            figs[f"{identifier}_{condition}"] = fig
    return figs, pd.DataFrame(plot_data)

# --- Streamlit UI ---

def main():
    st.set_page_config(layout="wide", page_title="PHASE Explorer App")
    st.title("🔬 Interactive PHASE Data Visualization App")
    st.markdown("Explore ecNGS data with interactive plots. Hover to see values, zoom into signatures, and customize your analysis.")
    if 'analysis_complete' not in st.session_state:
        st.session_state.analysis_complete = False
        st.session_state.results = {}
        st.session_state.custom_titles = {}
    with st.sidebar:
        st.header("⚙️ Inputs & Settings")
        st.subheader("1. File & Folder Paths")
        genome_path = st.text_input("Path to Reference Genome", placeholder="/path/to/your/genome.fa")
        indexing_folder = st.text_input("Path to Indexing Folder", placeholder="/path/to/your/data_folder")
        st.subheader("2. Data File Selection")
        base_files, all_conditions, all_identifiers = [], [], []
        mutations_dir = ""
        if indexing_folder and os.path.isdir(indexing_folder):
            bases_dir = os.path.join(indexing_folder, "bases")
            mutations_dir = os.path.join(indexing_folder, "mutations")
            base_files = find_data_files(bases_dir, "bases")
            id_pattern, cond_pattern = re.compile(r'^([^_]+)'), re.compile(r'[^_]+_(.+?)\.bc\d+')
            conditions_found, identifiers_found = set(), set()
            for f in base_files:
                basename = os.path.basename(f)
                id_match, cond_match = id_pattern.match(basename), cond_pattern.match(basename)
                if id_match: identifiers_found.add(id_match.group(1))
                if cond_match: conditions_found.add(cond_match.group(1))
            all_identifiers = sorted(list(identifiers_found))
            all_conditions = sorted(list(conditions_found), key=flexible_sort_key)
            if not base_files: st.warning("No `.bases` files found.")
            if not os.path.exists(mutations_dir): st.warning("`mutations` subdirectory not found.")
        elif indexing_folder: st.error("Indexing Folder path does not exist.")
        filter_identifier = st.selectbox("Filter by Identifier (Optional)", options=["All"] + all_identifiers)
        filtered_base_files = [f for f in base_files if filter_identifier == "All" or os.path.basename(f).startswith(filter_identifier + '_')]
        selected_base_filenames = st.multiselect("Select .bases files", [os.path.basename(f) for f in filtered_base_files], default=[os.path.basename(f) for f in filtered_base_files], disabled=not base_files)
        selected_bases_paths = [os.path.join(indexing_folder, "bases", f) for f in selected_base_filenames]
        st.subheader("3. Analysis Options")
        
        # --- Smart Control Selection UI ---
        control_options = [None]
        if 'Control' in all_identifiers:
            control_options.append("Identifier: Control (Time-Matched)")
        control_options.extend(all_conditions)
        control_selection = st.selectbox("Select Control", options=control_options, format_func=lambda x: 'None (Disable Stats/Subtraction)' if x is None else x)

        overlay_plots = st.checkbox("Overlay identifiers for comparison", value=False)
        st.subheader("4. Plot Customization")
        sync_yaxis = st.checkbox("Synchronize Y-axis for Summary plots", value=True)
        with st.expander("Customize Plot Text"):
            st.session_state.custom_titles['response_title'] = st.text_input("Response Plot: Title", placeholder="Default")
            st.session_state.custom_titles['response_xaxis'] = st.text_input("Response Plot: X-Axis", placeholder="Condition")
            st.session_state.custom_titles['response_yaxis'] = st.text_input("Response Plot: Y-Axis", placeholder="Mutation Frequency")
            st.markdown("---"); st.session_state.custom_titles['groups_title'] = st.text_input("Mutation Group Plot: Title", placeholder="Default"); st.session_state.custom_titles['groups_yaxis'] = st.text_input("Mutation Group Plot: Y-Axis", placeholder="Mutation Frequency")
            st.markdown("---"); st.session_state.custom_titles['sig_title'] = st.text_input("96-Channel Signature: Main Title", placeholder="Default"); st.session_state.custom_titles['sig_xaxis'] = st.text_input("96-Channel Signature: X-Axis", placeholder="Trinucleotide Context"); st.session_state.custom_titles['sig_yaxis'] = st.text_input("96-Channel Signature: Y-Axis", placeholder="Default")
        download_format = st.selectbox("Download Format", ["png", "pdf", "svg"], index=0, key="download_format")
        st.subheader("5. Run Analysis")
        col1, col2 = st.columns(2)
        with col1: analyze_button = st.button("🚀 Generate", type="primary", use_container_width=True)
        with col2: clear_button = st.button("❌ Clear Results", use_container_width=True)
        with st.expander("ℹ️ Help & Information", expanded=False):
            st.markdown("""...""") 
    if clear_button:
        st.session_state.analysis_complete = False; st.session_state.results = {}; st.session_state.custom_titles = {}; st.rerun()
    if analyze_button:
        with st.spinner("Analyzing data and generating plots..."):
            try:
                fai_path = f"{genome_path}.fai"
                if os.path.exists(fai_path): os.remove(fai_path)
                reference = Fasta(genome_path, as_raw=True, sequence_always_upper=True)
                baseline_counts, total_baseline = get_baseline_context_counts(reference)
                df = process_files(selected_bases_paths, mutations_dir, reference)
                if df.empty:
                    st.error("No valid data processed. Check file formats/selections."); st.session_state.analysis_complete = False; return
                y_limits = {}
                if sync_yaxis:
                    agg_df_resp = df.groupby(['identifier', 'condition'])['mutation_freq'].agg(['mean', sem]).reset_index()
                    y_limits['response'] = (agg_df_resp['mean'] + agg_df_resp['sem'].fillna(0)).max() * 1.15
                    max_group_freq = 0
                    grouped = df.groupby(['identifier', 'condition'])
                    for _, group in grouped:
                        total_bases = group['bases'].sum()
                        if total_bases > 0:
                            group_counts = {}
                            for _, row in group.iterrows():
                                for grp, count in row['mutation_group_counts'].items(): group_counts[grp] = group_counts.get(grp, 0) + count
                            for grp in DESIRED_GROUP_ORDER:
                                freq = group_counts.get(grp, 0) / total_bases
                                if freq > max_group_freq: max_group_freq = freq
                    y_limits['groups'] = max_group_freq * 1.15
                titles = st.session_state.custom_titles
                st.session_state.results['dose_response_figs'], st.session_state.results['dose_response_df'] = plot_dose_response_plotly(df, overlay=overlay_plots, control_selection=control_selection, y_max=y_limits.get('response'), custom_title=titles.get('response_title'), xaxis_title=titles.get('response_xaxis'), yaxis_title=titles.get('response_yaxis'))
                st.session_state.results['mutation_group_figs'], st.session_state.results['mutation_group_df'] = plot_mutation_groups_plotly(df, control_condition=control_selection, y_max=y_limits.get('groups'), custom_title=titles.get('groups_title'), yaxis_title=titles.get('groups_yaxis'))
                st.session_state.results['traditional_figs'], st.session_state.results['traditional_norm_df'] = plot_96_channel_signature_plotly(df, baseline_counts, total_baseline, mode='traditional', control_selection=control_selection, custom_title=titles.get('sig_title'), xaxis_title=titles.get('sig_xaxis'), yaxis_title=titles.get('sig_yaxis'))
                st.session_state.results['new_norm_figs'], st.session_state.results['new_norm_df'] = plot_96_channel_signature_plotly(df, baseline_counts, total_baseline, mode='new', control_selection=control_selection, custom_title=titles.get('sig_title'), xaxis_title=titles.get('sig_xaxis'), yaxis_title=titles.get('sig_yaxis'))
                st.session_state.results['subtracted_figs'], st.session_state.results['new_norm_diff_df'] = plot_96_channel_signature_plotly(df, baseline_counts, total_baseline, mode='subtracted', control_selection=control_selection, custom_title=titles.get('sig_title'), xaxis_title=titles.get('sig_xaxis'), yaxis_title=titles.get('sig_yaxis'))
                st.session_state.results['raw_df'] = df
                st.session_state.analysis_complete = True
            except Exception as e:
                st.error(f"An unexpected error occurred during analysis: {e}")
                st.session_state.analysis_complete = False
    if st.session_state.analysis_complete:
        st.success("✅ Analysis complete!")
        download_format_selected = st.session_state.get('download_format', 'png')
        zip_buffer = BytesIO()
        with zipfile.ZipFile(zip_buffer, 'a', zipfile.ZIP_DEFLATED, False) as zip_file:
            results_to_zip = st.session_state.results
            for title, data in results_to_zip.items():
                if title.endswith('_figs') and data:
                    folder_name = title.replace('_figs', '')
                    for key, fig in data.items():
                        plot_buffer = save_plotly_fig_to_buffer(fig, download_format_selected)
                        zip_file.writestr(f"plots/{folder_name}/{key}.{download_format_selected}", plot_buffer.getvalue())
                elif title.endswith('_df') and data is not None and not data.empty:
                    csv_bytes = to_csv(data)
                    zip_file.writestr(f"data/{title.replace('_df', '_data')}.csv", csv_bytes)
        st.download_button("📥 Download All Results (.zip)", zip_buffer.getvalue(), "ecNGS_analysis_results.zip", "application/zip", use_container_width=True)
        st.markdown("---")
        tab1, tab2, tab3, tab4 = st.tabs(["📊 Summary Plots", "🔬 96-Channel Signatures", "📄 Raw Data", "📈 Aggregated Data"])
        with tab1:
            st.header("Response by Condition")
            dose_figs = st.session_state.results.get('dose_response_figs')
            if dose_figs:
                for key, fig in dose_figs.items():
                    st.plotly_chart(fig, use_container_width=True)
                    st.download_button(f"Download Plot as .{download_format_selected.upper()}", save_plotly_fig_to_buffer(fig, download_format_selected), f"response_{key}.{download_format_selected}", f"image/{download_format_selected}", key=f"dl_dose_{key}")
            st.header("Mutation Frequency by Group")
            group_figs = st.session_state.results.get('mutation_group_figs')
            if group_figs:
                for key, fig in group_figs.items():
                    st.plotly_chart(fig, use_container_width=True)
                    st.download_button(f"Download Plot as .{download_format_selected.upper()}", save_plotly_fig_to_buffer(fig, download_format_selected), f"mutation_group_{key}.{download_format_selected}", f"image/{download_format_selected}", key=f"dl_group_{key}")
        with tab2:
            st.header("Traditional Normalization")
            traditional_figs = st.session_state.results.get('traditional_figs')
            if traditional_figs:
                for key, fig in traditional_figs.items():
                    st.plotly_chart(fig, use_container_width=True)
                    st.download_button(f"Download Plot as .{download_format_selected.upper()}", save_plotly_fig_to_buffer(fig, download_format_selected), f"traditional_{key}.{download_format_selected}", f"image/{download_format_selected}", key=f"dl_trad_{key}")
            st.header("New Normalization")
            new_figs = st.session_state.results.get('new_norm_figs')
            if new_figs:
                for key, fig in new_figs.items():
                    st.plotly_chart(fig, use_container_width=True)
                    st.download_button(f"Download Plot as .{download_format_selected.upper()}", save_plotly_fig_to_buffer(fig, download_format_selected), f"new_norm_{key}.{download_format_selected}", f"image/{download_format_selected}", key=f"dl_new_{key}")
            st.header("New Normalization (Control Subtracted)")
            sub_figs = st.session_state.results.get('subtracted_figs')
            if sub_figs:
                for key, fig in sub_figs.items():
                    st.plotly_chart(fig, use_container_width=True)
                    st.download_button(f"Download Plot as .{download_format_selected.upper()}", save_plotly_fig_to_buffer(fig, download_format_selected), f"subtracted_{key}.{download_format_selected}", f"image/{download_format_selected}", key=f"dl_sub_{key}")
        with tab3:
            st.header("Raw Processed Data")
            raw_df = st.session_state.results.get('raw_df')
            if raw_df is not None and not raw_df.empty:
                st.dataframe(raw_df.drop(columns=['mutation_group_counts', 'mutation_context_96_counts'], errors='ignore'))
                st.download_button("Download Raw Data as CSV", to_csv(raw_df), "raw_data.csv", 'text/csv', key="dl_raw_csv")
        with tab4:
            st.header("Aggregated Data for Plots")
            df_ui_map = {
                "Response Data": (st.session_state.results.get('dose_response_df'), "response_data.csv"),
                "Mutation Group Data": (st.session_state.results.get('mutation_group_df'), "mutation_group_data.csv"),
                "Traditional Normalization Data": (st.session_state.results.get('traditional_norm_df'), "traditional_norm_data.csv"),
                "New Normalization Data": (st.session_state.results.get('new_norm_df'), "new_norm_data.csv"),
                "Control-Subtracted Normalization Data": (st.session_state.results.get('new_norm_diff_df'), "subtracted_norm_data.csv"),
            }
            for title, (data, filename) in df_ui_map.items():
                with st.expander(f"View and Download: {title}"):
                    if data is not None and not data.empty:
                        st.dataframe(data)
                        st.download_button(f"Download {title} as CSV", to_csv(data), filename, 'text/csv', key=f"csv_dl_{filename}")
if __name__ == "__main__":
    main()

