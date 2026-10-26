from pathlib import Path
import numpy as np
import pandas as pd
import plotly.express as px
import streamlit as st

st.set_page_config(page_title='Genomics Research & NGS Quality Intelligence', page_icon='🧬', layout='wide')
st.markdown('''<style>
.block-container{padding-top:2rem}h1,h2,h3{letter-spacing:-.03em}
div[data-testid="stMetric"]{background:#1e293b!important;border:1px solid #334155;border-radius:14px;padding:18px}
div[data-testid="stMetricLabel"] *,div[data-testid="stMetricValue"] *{color:#f8fafc!important}
</style>''',unsafe_allow_html=True)
DATA=Path(__file__).parent/'data'/'GSE232336_RNA-seq_TPM.txt.gz'
@st.cache_data
def load_expression():
    return pd.read_csv(DATA,sep='\t',compression='gzip')
@st.cache_data
def metadata():
    return pd.DataFrame([{'sample':s,'cell_line':'HaCaT' if '_HaCaT_' in s else 'TE-7','condition':'IL-13' if '_IL-13_' in s else 'Unstimulated'} for s in load_expression().columns[2:]])

def qc_data():
    rng=np.random.default_rng(43)
    n=24
    return pd.DataFrame({'sample_id':[f'DEMO-{i:03d}' for i in range(1,n+1)],'q30_pct':rng.normal(88,7,n).clip(60,99),'alignment_pct':rng.normal(92,9,n).clip(55,99),'duplication_pct':rng.normal(24,12,n).clip(1,70),'reads_millions':rng.normal(32,9,n).clip(5,65)})

st.sidebar.title('🧬 Genomics Intelligence')
page=st.sidebar.radio('Explore', ['Overview','Research expression explorer','Advanced RNA-seq analysis','RNA-seq expression QC','NGS QC demonstration','NGS pipeline workflow','Variant review demonstration','Methods & provenance'])
st.sidebar.caption('Independent portfolio project • Not affiliated with or endorsed by Cincinnati Children’s')
if page=='Overview':
    st.title('Genomics Research & NGS Quality Intelligence')
    st.caption('Real public RNA-seq expression data + clearly labeled synthetic clinical workflow demonstrations')
    df=load_expression(); m=metadata()
    a,b,c,d=st.columns(4)
    a.metric('Genes',f'{len(df):,}');b.metric('Research samples',len(m));c.metric('Cell lines',m.cell_line.nunique());d.metric('Missing values',int(df.isna().sum().sum()))
    st.info('Research data: NCBI GEO GSE232336, processed TPM values. Clinical QC and variant review are synthetic educational demonstrations—not clinical data or validated software.')
    st.subheader('Study design')
    st.dataframe(m.groupby(['cell_line','condition']).size().reset_index(name='replicates'),hide_index=True,use_container_width=True)
    st.subheader('Research highlight: CCL26 expression')
    g=df[df.gene_name=='CCL26']
    if not g.empty:
        long=g.melt(id_vars=['gene_id','gene_name'],var_name='sample',value_name='TPM').merge(m,on='sample')
        fig=px.box(long,x='cell_line',y='TPM',color='condition',points='all',title='CCL26 by cell line and stimulation',labels={'TPM':'Expression (TPM)'})
        st.plotly_chart(fig,use_container_width=True)
elif page=='Research expression explorer':
    st.title('Public RNA-seq research explorer')
    st.caption('NCBI GEO GSE232336 • processed gene expression in TPM • 10 samples')
    df=load_expression();m=metadata()
    st.write('**Study groups**');st.dataframe(m,hide_index=True,use_container_width=True)
    options=df.gene_name.dropna().astype(str).unique().tolist()
    defaults=[x for x in ['CCL26','TSPAN6','DPM1','IL13RA1'] if x in options]
    genes=st.multiselect('Choose genes',options,default=defaults,max_selections=30)
    if genes:
        chosen=df[df.gene_name.isin(genes)].copy()
        long=chosen.melt(id_vars=['gene_id','gene_name'],value_vars=m['sample'].tolist(),var_name='sample',value_name='TPM').merge(m,on='sample')
        st.subheader('Expression across study conditions')
        gene=st.selectbox('Gene for group comparison',genes)
        fig=px.box(long[long.gene_name==gene],x='cell_line',y='TPM',color='condition',points='all',labels={'TPM':'Expression (TPM)'})
        st.plotly_chart(fig,use_container_width=True)
        st.subheader('Replicate-level expression profiles')
        profile=long.groupby(['gene_name','cell_line','condition','sample'],as_index=False)['TPM'].mean()
        profile['group']=profile['cell_line']+' · '+profile['condition']
        fig=px.strip(profile,x='gene_name',y='TPM',color='group',hover_data=['sample'],
                     title='Selected genes across biological groups (individual replicates)',
                     labels={'gene_name':'Gene','TPM':'Expression (TPM)'})
        fig.update_layout(xaxis_tickangle=-35)
        st.plotly_chart(fig,use_container_width=True)
        st.subheader('Expression heatmap (log₂[TPM + 1])')
        matrix=chosen.drop_duplicates('gene_name').set_index('gene_name')[m['sample'].tolist()]
        fig=px.imshow(np.log2(matrix+1),aspect='auto',labels={'x':'Sample','y':'Gene','color':'log₂(TPM+1)'},color_continuous_scale='Blues')
        st.plotly_chart(fig,use_container_width=True)
        st.subheader('Exploratory condition contrasts')
        rows=[]
        for cell in ['HaCaT','TE-7']:
            control=m[(m.cell_line==cell)&(m.condition=='Unstimulated')]['sample'].tolist()
            treated=m[(m.cell_line==cell)&(m.condition=='IL-13')]['sample'].tolist()
            sub=chosen[['gene_name']].copy();sub['cell_line']=cell
            sub['unstimulated_mean_TPM']=chosen[control].mean(axis=1)
            sub['IL13_mean_TPM']=chosen[treated].mean(axis=1)
            sub['exploratory_log2_ratio']=np.log2((sub['IL13_mean_TPM']+1)/(sub['unstimulated_mean_TPM']+1))
            rows.append(sub)
        contrasts=pd.concat(rows,ignore_index=True)
        st.dataframe(contrasts,use_container_width=True,hide_index=True)
        st.subheader('Exploratory IL-13 response')
        fig=px.bar(contrasts,x='gene_name',y='exploratory_log2_ratio',color='cell_line',barmode='group',
                   title='Descriptive log₂ mean TPM ratio (IL-13 / unstimulated; +1 pseudocount)',
                   labels={'gene_name':'Gene','exploratory_log2_ratio':'log₂ ratio'})
        fig.add_hline(y=0,line_dash='dash',line_color='#94a3b8')
        st.plotly_chart(fig,use_container_width=True)
        st.download_button('Download exploratory contrasts CSV',contrasts.to_csv(index=False),'GSE232336_exploratory_contrasts.csv','text/csv')
        st.warning('Ratios use a +1 TPM pseudocount and are descriptive only. No significance testing, differential-expression inference, or clinical interpretation is performed.')
elif page=='Advanced RNA-seq analysis':
    st.title('Advanced RNA-seq analysis')
    st.caption('Real processed TPM data · NCBI GEO GSE232336 · exploratory, descriptive comparisons only')
    df=load_expression()
    m=metadata()
    sample_cols=m['sample'].tolist()
    # Gene symbols can be duplicated across different gene identifiers; aggregate
    # for browsing and clearly disclose this limitation.
    expr=df.copy()
    expr['gene_name']=expr['gene_name'].fillna('').astype(str).str.strip()
    expr=expr[expr['gene_name']!=''].copy()
    expr[sample_cols]=expr[sample_cols].apply(pd.to_numeric,errors='coerce')
    expr=expr.groupby('gene_name',as_index=False)[sample_cols].sum(min_count=1)
    st.info('Gene symbols with multiple source rows are summed for this exploratory browser. For rigorous analyses, retain stable gene IDs and use appropriate count-level workflows.')
    st.subheader('1 · Gene search and replicate-level expression')
    gene_options=expr['gene_name'].tolist()
    default_gene='CCL26' if 'CCL26' in gene_options else gene_options[0]
    selected=st.selectbox('Search a gene symbol',gene_options,index=gene_options.index(default_gene))
    gene_row=expr.loc[expr.gene_name==selected].iloc[0]
    long=m.copy()
    long['TPM']=long['sample'].map(gene_row[sample_cols].to_dict()).astype(float)
    long['group']=long['cell_line']+' · '+long['condition']
    left,right=st.columns([3,2])
    with left:
        fig=px.strip(long,x='group',y='TPM',color='condition',hover_data=['sample','cell_line'],
                     title=f'{selected}: individual sample TPM values',
                     labels={'group':'Cell line · condition','TPM':'TPM'},
                     color_discrete_map={'IL-13':'#38bdf8','Unstimulated':'#f59e0b'})
        fig.update_traces(marker=dict(size=12,line=dict(width=1,color='#e2e8f0')))
        fig.update_layout(xaxis_tickangle=-20)
        st.plotly_chart(fig,use_container_width=True)
    with right:
        summary=long.groupby(['cell_line','condition'],as_index=False).agg(
            n=('TPM','size'),mean_TPM=('TPM','mean'),median_TPM=('TPM','median'),sd_TPM=('TPM','std'))
        st.markdown('**Group summary**')
        st.dataframe(summary.round(3),hide_index=True,use_container_width=True)
        st.caption('Each dot is a provided sample replicate. Replicate type is not independently verified; error estimates are descriptive.')
    st.download_button('Download selected gene replicate data',long.to_csv(index=False),
                       file_name=f'GSE232336_{selected}_replicates.csv',mime='text/csv')
    st.divider()
    st.subheader('2 · Explore genome-wide IL-13 response rankings')
    a,b,c=st.columns([1,1,1])
    cell=a.selectbox('Cell line',['HaCaT','TE-7'],key='advanced_cell')
    rank_mode=b.selectbox('Rank genes by',['Largest increase','Largest decrease','Largest absolute change'],key='advanced_rank')
    min_expr=c.number_input('Minimum mean TPM in either condition',min_value=0.0,max_value=100000.0,value=1.0,step=1.0)
    control=m.loc[(m.cell_line==cell)&(m.condition=='Unstimulated'),'sample'].tolist()
    treated=m.loc[(m.cell_line==cell)&(m.condition=='IL-13'),'sample'].tolist()
    ranking=expr[['gene_name']].copy()
    ranking['unstimulated_mean_TPM']=expr[control].mean(axis=1)
    ranking['IL13_mean_TPM']=expr[treated].mean(axis=1)
    ranking['exploratory_log2_ratio']=np.log2((ranking['IL13_mean_TPM']+1)/(ranking['unstimulated_mean_TPM']+1))
    ranking=ranking[ranking[['unstimulated_mean_TPM','IL13_mean_TPM']].max(axis=1)>=min_expr].dropna()
    if rank_mode=='Largest increase':
        ranking=ranking.sort_values('exploratory_log2_ratio',ascending=False)
    elif rank_mode=='Largest decrease':
        ranking=ranking.sort_values('exploratory_log2_ratio',ascending=True)
    else:
        ranking=ranking.assign(abs_change=ranking.exploratory_log2_ratio.abs()).sort_values('abs_change',ascending=False).drop(columns='abs_change')
    top_n=st.slider('Number of genes to show',5,40,15,5)
    top=ranking.head(top_n).copy()
    if top.empty:
        st.warning('No genes match this expression filter. Try lowering the minimum mean TPM.')
    else:
        fig=px.bar(top.iloc[::-1],x='exploratory_log2_ratio',y='gene_name',orientation='h',
                   color='exploratory_log2_ratio',color_continuous_scale='RdBu',color_continuous_midpoint=0,
                   labels={'exploratory_log2_ratio':'Descriptive log₂ ratio (+1 TPM)','gene_name':'Gene'},
                   title=f'{cell}: {rank_mode.lower()} after IL-13 stimulation')
        fig.add_vline(x=0,line_dash='dash',line_color='#94a3b8')
        fig.update_layout(coloraxis_showscale=False,height=max(440,top_n*24))
        st.plotly_chart(fig,use_container_width=True)
        st.dataframe(top.round(4),hide_index=True,use_container_width=True)
    st.download_button('Download full filtered ranking CSV',ranking.to_csv(index=False),
                       file_name=f'GSE232336_{cell}_exploratory_gene_ranking.csv',mime='text/csv')
    st.warning('Descriptive log₂ ratio = log₂((mean IL-13 TPM + 1) / (mean unstimulated TPM + 1)). No raw counts, p-values, FDR, or differential-expression significance are calculated. TPM is not a substitute for raw-count differential-expression modeling.')
elif page=='RNA-seq expression QC':
    st.title('RNA-seq expression QC & sample relationships')
    st.caption('NCBI GEO GSE232336 · exploratory analysis of processed TPM values · not FASTQ-level QC')
    st.info('PCA, correlations, and expression distributions describe processed expression profiles. They do not measure sequencing quality (Q30, adapter content, mapping rate) or establish batch effects or biological significance.')
    m=metadata()
    sample_cols=m['sample'].tolist()
    raw=load_expression()
    expression=raw[sample_cols].apply(pd.to_numeric,errors='coerce').clip(lower=0)
    expression=expression.replace([np.inf,-np.inf],np.nan)
    expression=expression.fillna(0)
    log_expr=np.log2(expression+1)
    st.subheader('1 · Principal component analysis (PCA)')
    min_mean=st.number_input('Minimum mean TPM for PCA/correlation',min_value=0.0,max_value=100000.0,value=1.0,step=1.0,help='Filters low-abundance genes before sample relationship analysis.')
    variable_n=st.slider('Most variable genes to include',100,5000,1000,100)
    filtered=log_expr.loc[expression.mean(axis=1)>=min_mean]
    if len(filtered)<3:
        st.warning('Too few genes remain. Lower the minimum mean TPM threshold.')
    else:
        variability=filtered.var(axis=1)
        selected_matrix=filtered.loc[variability.nlargest(min(variable_n,len(filtered))).index]
        # PCA: samples x genes, center each gene across samples; use SVD (no sklearn dependency).
        X=selected_matrix.to_numpy(dtype=float).T
        X=X-X.mean(axis=0,keepdims=True)
        U,S,Vt=np.linalg.svd(X,full_matrices=False)
        scores=U[:,:2]*S[:2]
        variance=S**2
        explained=variance[:2]/variance.sum()*100 if variance.sum()>0 else np.zeros(2)
        points=m.copy()
        points['PC1']=scores[:,0]
        points['PC2']=scores[:,1]
        fig=px.scatter(points,x='PC1',y='PC2',color='cell_line',symbol='condition',hover_name='sample',
                       hover_data=['cell_line','condition'],
                       labels={'PC1':f'PC1 ({explained[0]:.1f}% variance)','PC2':f'PC2 ({explained[1]:.1f}% variance)'},
                       title=f'PCA using {len(selected_matrix):,} variable genes')
        fig.update_traces(marker={'size':13,'line':{'width':1,'color':'#e2e8f0'}})
        st.plotly_chart(fig,use_container_width=True)
        st.caption('PCA is computed from log₂(TPM + 1), centered by gene across the ten samples, using the most variable genes. The signs of PCA axes are arbitrary. No batch-effect correction is performed.')
        st.download_button('Download PCA sample coordinates',points.to_csv(index=False),'GSE232336_PCA_coordinates.csv','text/csv')
        st.subheader('2 · Sample-to-sample expression correlation')
        corr=selected_matrix.corr(method='pearson')
        short_names=[f"{r.cell_line} | {'IL13' if r.condition=='IL-13' else 'Ctrl'} | {i+1}" for i,r in m.iterrows()]
        fig=px.imshow(corr.to_numpy(),x=short_names,y=short_names,zmin=0,zmax=1,
                      color_continuous_scale='Blues',aspect='auto',
                      labels={'x':'Sample','y':'Sample','color':'Pearson r'},
                      title='Pearson correlation of log₂(TPM + 1) profiles')
        fig.update_layout(height=650)
        st.plotly_chart(fig,use_container_width=True)
        st.caption('Correlation is descriptive and depends on the chosen gene filter; high correlation alone does not establish replicate validity.')
        st.subheader('3 · Within-group sample concordance')
        group_rows=[]
        for (cell,condition),g in m.groupby(['cell_line','condition'],sort=False):
            cols=g['sample'].tolist()
            if len(cols)>1:
                vals=corr.loc[cols,cols].to_numpy()
                upper=vals[np.triu_indices(len(cols),1)]
                for i in range(len(cols)):
                    for j in range(i+1,len(cols)):
                        group_rows.append({'cell_line':cell,'condition':condition,'sample_1':cols[i],
                                           'sample_2':cols[j],'pearson_r':float(vals[i,j])})
        pairs=pd.DataFrame(group_rows)
        cutoff=st.slider('Exploratory low-correlation flag threshold (Pearson r)',0.5,1.0,0.9,0.01)
        pairs['flag']=np.where(pairs['pearson_r']<cutoff,'Review','No flag')
        st.dataframe(pairs.round(4),hide_index=True,use_container_width=True)
        st.caption('Flags are illustrative screening prompts, not validated QC decisions. Biological and technical replicate types are not independently verified.')
        st.download_button('Download within-group correlations',pairs.to_csv(index=False),'GSE232336_within_group_correlations.csv','text/csv')
    st.subheader('4 · Expression distributions across samples')
    st.caption('Per-sample log₂(TPM + 1) distributions for all source gene rows. Zero-expression rows are retained.')
    dist=log_expr.copy()
    # Plotly box traces directly to avoid melting ~580k values into a long dataframe.
    import plotly.graph_objects as go
    fig=go.Figure()
    for _,r in m.iterrows():
        fig.add_trace(go.Box(y=dist[r['sample']].to_numpy(),name=r['sample'].replace('RNA-seq_',''),
                             boxpoints=False))
    fig.update_layout(title='Log-transformed TPM distribution by sample',
                      yaxis_title='log₂(TPM + 1)',xaxis_title='Sample',height=540,
                      xaxis_tickangle=-40,showlegend=False)
    st.plotly_chart(fig,use_container_width=True)
    st.warning('This page uses public processed TPM data, not raw reads or counts. For FASTQ-level QC use FastQC/MultiQC; for rigorous differential expression use an appropriate raw-count statistical workflow. No formal outlier detection or validated assay thresholds are claimed.')
elif page=='NGS QC demonstration':
    st.title('NGS sample quality monitoring')
    st.warning('Synthetic sequencing metrics only. Thresholds are illustrative and not validated for any assay.')
    df=qc_data()
    with st.expander('Adjust demonstration QC thresholds',expanded=True):
        q=st.slider('Minimum Q30 (%)',60,99,80)
        a=st.slider('Minimum alignment (%)',50,99,85)
        dup=st.slider('Maximum duplication (%)',5,80,45)
    df['review_reasons']=df.apply(lambda r:', '.join(
        (['Low Q30'] if r.q30_pct<q else [])+
        (['Low alignment'] if r.alignment_pct<a else [])+
        (['High duplication'] if r.duplication_pct>dup else [])),axis=1)
    # A two-or-more-flag failure is a DEMO policy, not a laboratory rule.
    df['flags']=df['review_reasons'].apply(lambda x:0 if not x else len(x.split(', ')))
    df['status']=df['flags'].map(lambda n:'Pass' if n==0 else ('Review' if n==1 else 'Fail'))
    a1,a2,a3,a4=st.columns(4)
    a1.metric('Samples',len(df));a2.metric('Pass',int((df.status=='Pass').sum()))
    a3.metric('Review',int((df.status=='Review').sum()));a4.metric('Fail',int((df.status=='Fail').sum()))
    st.caption('Illustrative status logic: Pass = 0 flags; Review = 1 flag; Fail = 2+ flags. Not clinical release criteria.')
    status_filter=st.multiselect('Filter by QC status',['Pass','Review','Fail'],default=['Pass','Review','Fail'])
    visible=df[df.status.isin(status_filter)]
    left,right=st.columns([2,1])
    with left:
        fig=px.scatter(visible,x='q30_pct',y='alignment_pct',color='status',size='reads_millions',
            hover_name='sample_id',hover_data=['review_reasons','duplication_pct'],
            color_discrete_map={'Pass':'#22c55e','Review':'#f59e0b','Fail':'#ef4444'},
            labels={'q30_pct':'Q30 (%)','alignment_pct':'Alignment (%)'},title='Sequencing quality by sample')
        fig.add_vline(x=q,line_dash='dash',line_color='#94a3b8')
        fig.add_hline(y=a,line_dash='dash',line_color='#94a3b8')
        st.plotly_chart(fig,use_container_width=True)
    with right:
        counts=df.status.value_counts().reindex(['Pass','Review','Fail'],fill_value=0).reset_index()
        counts.columns=['status','count']
        st.plotly_chart(px.pie(counts,names='status',values='count',hole=.58,
            color='status',color_discrete_map={'Pass':'#22c55e','Review':'#f59e0b','Fail':'#ef4444'},
            title='QC status distribution'),use_container_width=True)
    st.subheader('Sample QC audit table')
    st.dataframe(visible.round(2),hide_index=True,use_container_width=True)
    st.download_button('Download synthetic QC audit CSV',df.to_csv(index=False),'synthetic_ngs_qc_audit.csv','text/csv')
    st.divider()
    st.subheader('🔎 Sample QC investigation panel')
    st.caption('Investigate one synthetic sample against the current demonstration thresholds. Recommendations are educational, not clinical instructions.')
    sample_id=st.selectbox('Choose a sample to investigate',df['sample_id'].tolist(),
        index=next((i for i, sid in enumerate(df['sample_id']) if df.iloc[i]['status']!='Pass'),0))
    sample=df.loc[df['sample_id']==sample_id].iloc[0]
    st.markdown(f"**{sample_id} — {sample['status']}** · {int(sample['flags'])} QC flag(s)")
    c1,c2,c3,c4=st.columns(4)
    c1.metric('Q30',f"{sample['q30_pct']:.1f}%",f"{sample['q30_pct']-q:+.1f} pp vs min")
    c2.metric('Alignment',f"{sample['alignment_pct']:.1f}%",f"{sample['alignment_pct']-a:+.1f} pp vs min")
    c3.metric('Duplication',f"{sample['duplication_pct']:.1f}%",f"{dup-sample['duplication_pct']:+.1f} pp headroom")
    c4.metric('Reads',f"{sample['reads_millions']:.1f}M")
    checks=[
        {'Metric':'Q30 (%)','Observed':round(float(sample['q30_pct']),2),'Demo threshold':f'≥ {q}%',
         'Result':'Flag' if sample['q30_pct']<q else 'Pass',
         'Investigation idea':'Review FastQC per-base quality, sequencing-cycle trends, and adapter content.'},
        {'Metric':'Alignment (%)','Observed':round(float(sample['alignment_pct']),2),'Demo threshold':f'≥ {a}%',
         'Result':'Flag' if sample['alignment_pct']<a else 'Pass',
         'Investigation idea':'Check reference/aligner settings, contamination indicators, read composition, and mapping summaries.'},
        {'Metric':'Duplication (%)','Observed':round(float(sample['duplication_pct']),2),'Demo threshold':f'≤ {dup}%',
         'Result':'Flag' if sample['duplication_pct']>dup else 'Pass',
         'Investigation idea':'Inspect library complexity and duplication metrics; interpret in assay context.'},
    ]
    st.dataframe(pd.DataFrame(checks),hide_index=True,use_container_width=True)
    flagged=[c for c in checks if c['Result']=='Flag']
    if flagged:
        st.warning('Flagged metrics: '+', '.join(c['Metric'] for c in flagged))
        st.markdown('**Suggested investigation checklist (demonstration)**')
        for item in flagged:
            st.checkbox(item['Investigation idea'],key=f"check_{sample_id}_{item['Metric']}")
    else:
        st.success('No flags under the currently selected demonstration thresholds.')
    st.info('Classification is based on the number of threshold flags, not an assay-validated release decision. A sample may require investigation even if these three metrics pass.')
    investigation=pd.DataFrame([{
        'sample_id':sample_id,'status':sample['status'],'flags':int(sample['flags']),
        'q30_pct':round(float(sample['q30_pct']),2),'minimum_q30_pct':q,
        'alignment_pct':round(float(sample['alignment_pct']),2),'minimum_alignment_pct':a,
        'duplication_pct':round(float(sample['duplication_pct']),2),'maximum_duplication_pct':dup,
        'reads_millions':round(float(sample['reads_millions']),2),
        'flag_reasons':sample['review_reasons'] or 'None',
        'data_type':'Synthetic demonstration; not validated clinical QC'
    }])
    st.download_button('Download selected sample investigation CSV',investigation.to_csv(index=False),
        file_name=f'{sample_id}_qc_investigation_demo.csv',mime='text/csv')
elif page=='NGS pipeline workflow':
    st.title('NGS pipeline execution & traceability')
    st.warning('All run IDs, stages, timings, and statuses are synthetic. This is a workflow demonstration, not an actual Nextflow or SLURM connection.')
    stages=[('FASTQ intake','Validate files and sample sheet','FASTQ'),
            ('FastQC / MultiQC','Read quality and adapter review','QC reports'),
            ('Alignment','STAR for RNA or BWA for DNA (assay dependent)','BAM'),
            ('SAMtools QC','Sorting, indexing, alignment metrics','BAM + metrics'),
            ('Variant calling','Assay-specific calling and filters','VCF'),
            ('Annotation / review','Annotation and documented review','Annotated VCF'),
            ('Report / archive','Controlled output and provenance','Review report')]
    run=st.selectbox('Select demonstration run',['DEMO-RUN-1042','DEMO-RUN-1043','DEMO-RUN-1044'])
    scenario={'DEMO-RUN-1042':6,'DEMO-RUN-1043':3,'DEMO-RUN-1044':1}[run]
    st.progress((scenario+1)/len(stages),text=f'{scenario+1} of {len(stages)} stages reached (synthetic)')
    rows=[]
    for i,(name,description,output) in enumerate(stages):
        status='Complete' if i<scenario else ('Running' if i==scenario else 'Queued')
        if run=='DEMO-RUN-1044' and i==scenario:status='Needs review'
        rows.append({'stage':f'{i+1:02d} · {name}','status':status,'purpose':description,'expected_artifact':output})
    pipeline=pd.DataFrame(rows)
    st.subheader('Pipeline stages')
    for r in rows:
        symbol={'Complete':'✅','Running':'🔵','Queued':'⚪','Needs review':'🟠'}[r['status']]
        with st.expander(f"{symbol} {r['stage']} — {r['status']}",expanded=r['status'] in ['Running','Needs review']):
            st.write(r['purpose']);st.caption('Expected artifact: '+r['expected_artifact'])
    st.subheader('Execution metadata (demonstration)')
    st.dataframe(pd.DataFrame([{'run_id':run,'orchestrator':'Nextflow (illustrative)',
        'scheduler':'SLURM (illustrative)','reference_build':'Not specified — demonstration',
        'sample_sheet':'Synthetic','execution_state':rows[scenario]['status']}]),hide_index=True,use_container_width=True)
    st.download_button('Download synthetic pipeline audit CSV',pipeline.to_csv(index=False),f'{run}_stages.csv','text/csv')
    st.info('Interview talking point: In a real implementation, Nextflow trace/timeline files, FastQC/MultiQC outputs, and SLURM job logs would provide provenance and operational status. This prototype does not execute those tools.')
elif page=='Variant review demonstration':
    st.title('Variant review — synthetic demonstration')
    st.warning('Invented variants for UI demonstration only; no patient, VCF, pathogenicity, or clinical interpretation data.')
    rng=np.random.default_rng(19)
    genes=['TP53','KRAS','BRAF','EGFR','PIK3CA','ALK']
    df=pd.DataFrame({'record_id':[f'SYN-V{i:03d}' for i in range(1,31)],'gene':rng.choice(genes,30),'variant_type':rng.choice(['SNV','Insertion','Deletion'],30),'allele_fraction':rng.uniform(.02,.9,30).round(3),'depth':rng.integers(30,900,30)})
    g=st.multiselect('Filter genes',genes,default=genes)
    sub=df[df.gene.isin(g)]
    st.plotly_chart(px.histogram(sub,x='gene',color='variant_type',barmode='group'),use_container_width=True)
    st.dataframe(sub,hide_index=True,use_container_width=True)
    st.download_button('Download synthetic variant table',sub.to_csv(index=False),'synthetic_variants.csv','text/csv')
else:
    st.title('Methods, provenance & limitations')
    st.markdown('''**Real data:** NCBI Gene Expression Omnibus accession [GSE232336](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE232336), file `GSE232336_RNA-seq_TPM.txt.gz`. The matrix has gene identifiers, gene symbols and TPM expression across HaCaT and TE-7 cell lines with unstimulated and IL-13-treated replicates. Data are a fixed public research snapshot, **not live hospital data**.\n\n**Research processing:** The app reads published processed TPM values. Heatmaps use log₂(TPM+1). Exploratory log₂ ratios use mean TPM with +1 pseudocount. These are not formal differential expression statistics.\n\n**Clinical demonstration:** QC, pipeline runs, and variant tables are seeded synthetic records, generated deterministically for demonstration. Adjustable QC thresholds are not assay-specific validated clinical cutoffs.\n\n**Clinical readiness gap:** Production use would require validated assays, controlled reference datasets, documented sensitivity/specificity and reproducibility, change control, security, audit logs, appropriate regulatory oversight, and laboratory sign-off.\n\n**Project status:** Independent educational prototype; not an official Cincinnati Children’s application or diagnostic tool.''')
