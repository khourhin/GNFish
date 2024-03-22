import streamlit as st

st.title("GNFish")
st.write(
    """A python library and scripts following the pipeline detailed in (DOI of article yet to be published) for:\n
(1) “fishing” specific gene sequences of interest from genomic data available in NCBI databases\n
(2) processing and depuration of retrieved sequences\n
(3) production of a multiple sequence alignment\n
(4) selection of best-fit model of evolution\n
(5) solid reconstruction of a phylogenetic tree."""
)

st.header("Analysis")

st.file_uploader("Upload your input file here")
