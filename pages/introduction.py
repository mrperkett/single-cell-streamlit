import streamlit as st


def run():
    st.markdown("# Single-Cell Analysis Web Tool")
    st.image("resources/single-cell-banner.jpg")
    text = """## Introduction

Welcome to the single-cell analysis web application! 🦠🧬📊

This application allows you to work through standard preprocessing, quality control, and analysis for single-cell RNA-Seq experiments.

To run the analysis, navigate to each of the pages listed on the left side bar from top to bottom starting with *Load Data*. Each page will have instructions on how to complete the analysis for the steps.  Parameters that can be customized by the user are located on the left side bar. You can return to previous steps in the analysis to inspect plots or rerun the analysis with different parameters.  

> The datasets that are registered for analysis can be adjusted in the config file used when deploying this application.  See the [GitHub repository](https://github.com/mrperkett/single-cell-streamlit) for details.

> ⚠️ If you rerun a previously completed analysis step, it will require all downstream steps to also be rerun.  For example, if you rerun the *Normalization* step, you will need to rerun from the *Feature Selection* step on down.

## Resources

- [Scanpy preprocessing and clustering tutorial](https://scanpy.readthedocs.io/en/1.10.x/tutorials/basics/clustering.html)
- [Single-cell best practices](https://www.sc-best-practices.org/preamble.html)
    """
    st.markdown(text)


run()
