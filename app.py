import streamlit as st
import pandas as pd
import plotly.express as px

st.set_page_config(page_title="Metabolomics Data Hub", layout="wide")
st.title("LC-HRMS Data Processing & Visualization")

st.sidebar.header("1. Upload Data Tables")
quant_file = st.sidebar.file_uploader("Upload Quant Table (setac0305...quant.csv)", type=['csv'])
cmmc_file = st.sidebar.file_uploader("Upload CMMC Annotations (tsv)", type=['tsv', 'csv'])
cluster_file = st.sidebar.file_uploader("Upload Cluster Summary (tsv)", type=['tsv', 'csv'])
meta_file = st.sidebar.file_uploader("Upload Metadata (csv/xlsx)", type=['csv', 'xlsx'])

if quant_file and cmmc_file and cluster_file and meta_file:
    
    # Load data
    df_quant = pd.read_csv(quant_file)
    df_cmmc = pd.read_csv(cmmc_file, sep='\t') 
    df_cluster = pd.read_csv(cluster_file, sep='\t') 
    
    # Handle metadata as string to preserve formatting (e.g., sample codes with leading zeros)
    if meta_file.name.endswith('.csv'):
        df_meta = pd.read_csv(meta_file, dtype=str)
    else:
        df_meta = pd.read_excel(meta_file, dtype=str)

    st.header("Step 1 & 2: Merging and Annotation Filtering")
    
    # Outer Merging and Key Coalescing
    # 1st Outer Merge: Quant + CMMC
    df_merged = pd.merge(df_quant, df_cmmc, left_on='row ID', right_on='query_scan', how='outer')
    df_merged['row ID'] = df_merged['row ID'].combine_first(df_merged['query_scan'])
    
    # 2nd Outer Merge: Result + Cluster Summary
    df_merged = pd.merge(df_merged, df_cluster, left_on='row ID', right_on='cluster index', how='outer')
    df_merged['row ID'] = df_merged['row ID'].combine_first(df_merged['cluster index'])
    
    st.success(f"Tables successfully merged (Outer Mode)! Initial feature count: {len(df_merged)}")

    # ---------------------------------------------------------
    # UPDATE: Filter by Annotation using Compound_Name
    # ---------------------------------------------------------
    # Ensure columns exist just to avoid Streamlit crashing if a bad file is uploaded
    if 'Compound_Name' not in df_merged.columns:
        df_merged['Compound_Name'] = None
    if 'input_name' not in df_merged.columns:
        df_merged['input_name'] = None

    # Keeps rows where either 'Compound_Name' is not null OR 'input_name' is not null
    df_annotated = df_merged[df_merged['Compound_Name'].notna() | df_merged['input_name'].notna()].copy()
    st.write(f"Features remaining after retaining only annotated rows: {len(df_annotated)}")

    st.header("Step 3 & 4: Intensity and Blank Correction")
    
    # Isolate sample and blank columns using metadata
    sample_files = df_meta[df_meta['ATTRIBUTE_SampleType'] == 'Sample']['filename'].tolist()
    blank_files = df_meta[df_meta['ATTRIBUTE_SampleType'].astype(str).str.contains('Blank', na=False)]['filename'].tolist()
    
    sample_cols = [f"{f} Peak area" for f in sample_files if f"{f} Peak area" in df_annotated.columns]
    blank_cols = [f"{f} Peak area" for f in blank_files if f"{f} Peak area" in df_annotated.columns]

    # Calculate Max Intensities row by row
    df_annotated['Max_Sample_Intensity'] = df_annotated[sample_cols].max(axis=1)
    df_annotated['Max_Blank_Intensity'] = df_annotated[blank_cols].max(axis=1).fillna(0)

    # UI Inputs
    col1, col2 = st.columns(2)
    with col1:
        X = st.number_input("Blank Correction Factor (X)", min_value=0.0, value=3.0, step=0.5)
    with col2:
        Y = st.number_input("Minimum Sample Intensity (Y)", min_value=0.0, value=10000.0, step=1000.0)

    # Apply Filters
    condition_blank = df_annotated['Max_Sample_Intensity'] > (X * df_annotated['Max_Blank_Intensity'])
    condition_intensity = df_annotated['Max_Sample_Intensity'] > Y
    df_final = df_annotated[condition_blank & condition_intensity].copy()
    
    st.success(f"Final feature count after intensity filtering: {len(df_final)}")

    st.header("Step 5: Sample Visualization")
    
    if not df_final.empty:
        # Define the unified label for the plot dropdown and title
        def define_annotation(row):
            gnps = row.get('Compound_Name', None)
            cmmc = row.get('input_name', None)
            
            labels = []
            if pd.notna(gnps) and str(gnps).strip() != "":
                labels.append(f"GNPS: {gnps}")
            if pd.notna(cmmc) and str(cmmc).strip() != "":
                labels.append(f"CMMC: {cmmc}")
                
            if labels:
                return f"Feature {row['row ID']} | " + " | ".join(labels)
            return f"Feature {row['row ID']} | Unknown"
            
        df_final['Plot_Name'] = df_final.apply(define_annotation, axis=1)

        # Melt data for Plotly
        df_melted = df_final.melt(
            id_vars=['row ID', 'Plot_Name'], 
            value_vars=sample_cols, 
            var_name='Column_Name', 
            value_name='Intensity'
        )
        
        df_melted['filename'] = df_melted['Column_Name'].str.replace(' Peak area', '', regex=False)
        df_plot = df_melted.merge(df_meta, on='filename', how='left')

        # Dropdowns for user selection
        plot_col1, plot_col2 = st.columns(2)
        with plot_col1:
            selected_feature = st.selectbox("Select a Feature to View:", options=df_final['Plot_Name'].unique())
        with plot_col2:
            meta_attributes = [col for col in df_meta.columns if col.startswith('ATTRIBUTE_')]
            grouping_var = st.selectbox("Group Boxplots By:", options=meta_attributes)

        # Plotting
        data_to_plot = df_plot[df_plot['Plot_Name'] == selected_feature]
        
        fig = px.box(
            data_to_plot, 
            x=grouping_var, 
            y='Intensity', 
            color=grouping_var,
            points="all", 
            title=f"Distribution of {selected_feature}",
            labels={grouping_var: grouping_var.replace('ATTRIBUTE_', '')}
        )
        
        st.plotly_chart(fig, use_container_width=True)
        
        # Download Cleaned Data
        st.download_button(
            label="Download Cleaned Data Table",
            data=df_final.to_csv(index=False).encode('utf-8'),
            file_name='cleaned_metabolomics_data.csv',
            mime='text/csv'
        )
    else:
        st.warning("No features passed the blank and intensity thresholds. Try lowering X or Y.")
else:
    st.info("Please upload all four files in the sidebar to begin processing.")
